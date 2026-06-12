"""
Segmentation Refinement Utilities for Pyable

Generic functions for refining binary ROI masks using various strategies:
  - Watershed segmentation
  - Confidence-connected region growing
  - Geodesic active contour (level-set) boundary smoothing
  - Threshold level-set refinement
  - Laplacian level-set refinement
  - Shape detection level-set
  - Chan-Vese (region-based) segmentation
  - Otsu / multi-Otsu / Li / Yen / Triangle / Huang thresholding
  - Connected threshold region growing
  - Neighbourhood connected region growing
  - Isolated connected region growing
  - Morphological watershed from markers
  - N4 bias field correction (preprocessing)
  - Anisotropic diffusion smoothing (preprocessing)
  - Probability-gated expansion / shrinkage
  - STAPLE multi-segmentation fusion
  - Distance-constrained clipping
  - Mask subtraction / intersection
  - Morphological utilities (fill holes, filter components, distance map, shell)

All functions operate on SimpleITK images and return SimpleITK images.
The corresponding methods on Roiable / LabelMapable / Imaginable are thin wrappers.

Example:
    >>> from pyable import Roiable, Imaginable
    >>> roi = Roiable('mask.nii.gz')
    >>> img = Imaginable('scan.nii.gz')
    >>> roi.refineWatershed(img)
    >>> roi.refineGeodesicActiveContour(img, propagation=0.3, curvature=0.8)
    >>> roi.write('refined_mask.nii.gz')
"""

import numpy as np
import SimpleITK as sitk

from scipy.ndimage import (
    binary_dilation,
    binary_erosion,
    binary_fill_holes,
    distance_transform_edt,
    label as nd_label,
)


# ============================================================================
# HELPERS
# ============================================================================

def _to_sitk(obj):
    """Extract sitk.Image from an Imaginable-like object or return as-is."""
    if hasattr(obj, 'getImage') and callable(obj.getImage):
        return obj.getImage()
    return obj


def _spacing_zyx(image: sitk.Image) -> np.ndarray:
    """Return spacing in numpy (Z,Y,X) order."""
    return np.array(image.GetSpacing()[::-1])


def _ensure_same_grid(moving: sitk.Image, reference: sitk.Image,
                      interpolator=sitk.sitkNearestNeighbor) -> sitk.Image:
    """Resample *moving* onto *reference* grid if sizes differ."""
    if moving.GetSize() == reference.GetSize():
        return moving
    return sitk.Resample(moving, reference, sitk.Transform(),
                         interpolator, 0.0)


def _arr_to_sitk_binary(arr: np.ndarray, reference: sitk.Image) -> sitk.Image:
    """Convert a boolean/uint8 numpy array to a sitk binary image with geometry."""
    result = sitk.GetImageFromArray(arr.astype(np.uint8))
    result.CopyInformation(reference)
    return result


# ============================================================================
# MORPHOLOGICAL UTILITIES
# ============================================================================

def fill_binary_holes(roi: sitk.Image) -> sitk.Image:
    """
    Fill all holes inside a binary ROI using scipy binary_fill_holes.

    Unlike skimage remove_small_holes, this fills *all* enclosed cavities
    regardless of size.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI image.

    Returns
    -------
    sitk.Image
        ROI with all internal holes filled.
    """
    arr = sitk.GetArrayFromImage(roi) > 0
    filled = binary_fill_holes(arr).astype(np.uint8)
    return _arr_to_sitk_binary(filled, roi)


def filter_small_components(roi: sitk.Image, min_voxels: int = 50) -> sitk.Image:
    """
    Remove connected components smaller than *min_voxels*.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI image.
    min_voxels : int
        Minimum component size to keep.

    Returns
    -------
    sitk.Image
        Filtered binary ROI.
    """
    arr = sitk.GetArrayFromImage(roi) > 0
    if not np.any(arr):
        return roi

    labeled, n = nd_label(arr)
    filtered = np.zeros_like(arr)
    for i in range(1, n + 1):
        if np.sum(labeled == i) >= min_voxels:
            filtered |= (labeled == i)

    return _arr_to_sitk_binary(filtered, roi)


def keep_largest_component(roi: sitk.Image) -> sitk.Image:
    """
    Keep only the largest connected component.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI image.

    Returns
    -------
    sitk.Image
        Binary ROI containing only the largest component.
    """
    arr = sitk.GetArrayFromImage(roi) > 0
    if not np.any(arr):
        return roi

    labeled, n = nd_label(arr)
    if n <= 1:
        return roi

    sizes = [(i, np.sum(labeled == i)) for i in range(1, n + 1)]
    largest_label = max(sizes, key=lambda x: x[1])[0]
    result = (labeled == largest_label).astype(np.uint8)
    return _arr_to_sitk_binary(result, roi)


def compute_distance_map(roi: sitk.Image) -> sitk.Image:
    """
    Euclidean distance transform from the ROI surface.

    Zero inside the ROI, positive outside.  Uses image spacing.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI image.

    Returns
    -------
    sitk.Image
        Float32 distance map (mm).
    """
    arr = sitk.GetArrayFromImage(roi) > 0
    spacing = _spacing_zyx(roi)
    dist = distance_transform_edt(~arr, sampling=spacing).astype(np.float32)
    result = sitk.GetImageFromArray(dist)
    result.CopyInformation(roi)
    return result


def compute_shell(roi: sitk.Image, width_mm: float = 5.0) -> sitk.Image:
    """
    Compute a shell (annular band) around the ROI surface.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI image.
    width_mm : float
        Shell width in mm (extends outward from ROI surface).

    Returns
    -------
    sitk.Image
        Binary image of the shell region.
    """
    arr = sitk.GetArrayFromImage(roi) > 0
    spacing = _spacing_zyx(roi)
    mean_sp = float(np.mean(spacing))
    n_dilations = max(1, int(width_mm / mean_sp))

    struct = np.ones((3, 3, 3), dtype=bool)
    dilated = binary_dilation(arr, structure=struct, iterations=n_dilations)
    shell = dilated & ~arr
    return _arr_to_sitk_binary(shell, roi)


def compute_edge_map(image: sitk.Image, sigma: float = 1.0) -> sitk.Image:
    """
    Gradient-magnitude edge potential normalised to [0, 1].

    Parameters
    ----------
    image : sitk.Image
        Scalar image.
    sigma : float
        Gaussian sigma in mm for gradient computation.

    Returns
    -------
    sitk.Image
        Float32 edge map in [0, 1].
    """
    from scipy.ndimage import gaussian_gradient_magnitude

    arr = sitk.GetArrayFromImage(image).astype(np.float32)
    spacing = _spacing_zyx(image)

    # Normalise intensity to [0, 1]
    valid = arr[arr > 0]
    if len(valid) > 0:
        lo, hi = np.percentile(valid, [2, 98])
    else:
        lo, hi = 0.0, 1.0
    arr_norm = np.clip((arr - lo) / (hi - lo + 1e-6), 0, 1)

    sigma_vox = sigma / float(np.mean(spacing))
    grad = gaussian_gradient_magnitude(arr_norm, sigma=sigma_vox)
    if grad.max() > 0:
        grad = grad / grad.max()

    result = sitk.GetImageFromArray(grad.astype(np.float32))
    result.CopyInformation(image)
    return result


# ============================================================================
# WATERSHED REFINEMENT
# ============================================================================

def watershed_refine(
    roi: sitk.Image,
    image: sitk.Image,
    height_map: sitk.Image = None,
    erosion_iters: int = 3,
    dilation_iters: int = 5,
    min_voxels: int = 50,
) -> sitk.Image:
    """
    Refine a binary ROI using marker-based watershed segmentation.

    Markers are created from an eroded interior (foreground) and a dilated
    exterior (background).  The height map controls where boundaries form;
    by default the gradient magnitude of *image* is used.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI (initial segmentation).
    image : sitk.Image
        Reference intensity image (e.g. MRI scan).
    height_map : sitk.Image, optional
        Custom height/cost map for watershed.  If *None*, the gradient
        magnitude of *image* is computed automatically.
    erosion_iters : int
        Erosion iterations for foreground markers.
    dilation_iters : int
        Dilation iterations for background markers.
    min_voxels : int
        Remove components smaller than this after watershed.

    Returns
    -------
    sitk.Image
        Refined binary ROI.

    Example
    -------
    >>> from pyable.segmentation import watershed_refine
    >>> refined = watershed_refine(roi_sitk, scan_sitk)
    """
    from skimage.segmentation import watershed as sk_watershed

    image = _ensure_same_grid(image, roi, sitk.sitkLinear)
    roi_arr = sitk.GetArrayFromImage(roi) > 0

    # Build height map
    if height_map is not None:
        height_map = _ensure_same_grid(height_map, roi, sitk.sitkLinear)
        h_arr = sitk.GetArrayFromImage(height_map).astype(np.float32)
    else:
        h_arr = sitk.GetArrayFromImage(
            compute_edge_map(_ensure_same_grid(image, roi, sitk.sitkLinear))
        ).astype(np.float32)

    # Markers
    struct = np.ones((3, 3, 3), dtype=bool)
    foreground = binary_erosion(roi_arr, structure=struct, iterations=erosion_iters)
    background = ~binary_dilation(roi_arr, structure=struct, iterations=dilation_iters)

    if not np.any(foreground):
        # ROI too small to erode — return as-is
        return roi

    markers = np.zeros(roi_arr.shape, dtype=np.int32)
    markers[foreground] = 1
    markers[background] = 2

    labels = sk_watershed(h_arr, markers, mask=~background)
    refined = (labels == 1).astype(np.uint8)
    refined = binary_fill_holes(refined).astype(np.uint8)

    result = _arr_to_sitk_binary(refined, roi)
    result = filter_small_components(result, min_voxels=min_voxels)
    return result


# ============================================================================
# REGION GROWING REFINEMENT
# ============================================================================

def region_growing_refine(
    roi: sitk.Image,
    image: sitk.Image,
    multiplier: float = 2.5,
    neighborhood_radius: int = 1,
    n_iterations: int = 3,
    max_distance_mm: float = 10.0,
    n_seeds: int = 100,
    prob_map: sitk.Image = None,
    prob_threshold: float = 0.1,
    min_voxels: int = 50,
) -> sitk.Image:
    """
    Refine a binary ROI via confidence-connected region growing.

    Seeds are sampled from the eroded ROI interior.  Growth is constrained
    by an optional maximum distance from the original ROI surface and an
    optional probability map.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI (initial segmentation).
    image : sitk.Image
        Reference intensity image for intensity-based growing.
    multiplier : float
        Number of standard deviations for the intensity acceptance range.
    neighborhood_radius : int
        Radius for local statistics computation.
    n_iterations : int
        Region-growing iterations.
    max_distance_mm : float
        Maximum distance (mm) from the original ROI surface for growth.
        Set to *None* or 0 to disable.
    n_seeds : int
        Maximum number of seed points sampled from eroded interior.
    prob_map : sitk.Image, optional
        Probability map — grown region is intersected with
        ``prob_map > prob_threshold``.
    prob_threshold : float
        Threshold applied to *prob_map* (ignored when *prob_map* is None).
    min_voxels : int
        Remove components smaller than this.

    Returns
    -------
    sitk.Image
        Refined binary ROI.
    """
    image = _ensure_same_grid(image, roi, sitk.sitkLinear)
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    spacing = _spacing_zyx(roi)

    # Build seed points from eroded interior
    struct = np.ones((3, 3, 3), dtype=bool)
    eroded = binary_erosion(roi_arr, structure=struct, iterations=2)
    seed_indices = np.argwhere(eroded)

    if len(seed_indices) == 0:
        return roi

    # Sub-sample seeds
    if len(seed_indices) > n_seeds:
        idx = np.linspace(0, len(seed_indices) - 1, n_seeds, dtype=int)
        seed_indices = seed_indices[idx]

    # Convert (z,y,x) -> (x,y,z) for sitk
    seeds = [tuple(int(v) for v in reversed(s)) for s in seed_indices]

    # Run confidence-connected region growing
    cc = sitk.ConfidenceConnectedImageFilter()
    cc.SetSeedList(seeds)
    cc.SetMultiplier(multiplier)
    cc.SetNumberOfIterations(n_iterations)
    cc.SetInitialNeighborhoodRadius(neighborhood_radius)

    grown = cc.Execute(sitk.Cast(image, sitk.sitkFloat32))
    grown_arr = sitk.GetArrayFromImage(grown) > 0

    # Distance constraint
    if max_distance_mm and max_distance_mm > 0:
        dist_from_roi = distance_transform_edt(~roi_arr, sampling=spacing)
        grown_arr = grown_arr & (dist_from_roi <= max_distance_mm)

    # Probability constraint
    if prob_map is not None:
        prob_map = _ensure_same_grid(prob_map, roi, sitk.sitkLinear)
        prob_arr = sitk.GetArrayFromImage(prob_map)
        grown_arr = grown_arr & (prob_arr > prob_threshold)

    # Preserve eroded core
    grown_arr |= eroded

    grown_arr = binary_fill_holes(grown_arr).astype(np.uint8)
    result = _arr_to_sitk_binary(grown_arr, roi)
    result = filter_small_components(result, min_voxels=min_voxels)
    return result


# ============================================================================
# GEODESIC ACTIVE CONTOUR REFINEMENT
# ============================================================================

def geodesic_active_contour_refine(
    roi: sitk.Image,
    image: sitk.Image,
    propagation: float = 0.5,
    curvature: float = 0.3,
    advection: float = 1.0,
    iterations: int = 50,
    rms_tolerance: float = 0.001,
    allow_shrink: bool = False,
    speed_image: sitk.Image = None,
) -> sitk.Image:
    """
    Refine ROI boundaries using a geodesic active contour level-set.

    The contour evolves under propagation (balloon), curvature (smoothness)
    and advection (edge attraction) forces.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI used as the initial contour.
    image : sitk.Image
        Reference intensity image for edge computation.
    propagation : float
        Balloon force.  Positive → expand, negative → shrink.
    curvature : float
        Smoothing force (higher = smoother boundaries).
    advection : float
        Edge-attraction force.
    iterations : int
        Maximum GAC iterations.
    rms_tolerance : float
        Convergence threshold (lower = more iterations).
    allow_shrink : bool
        If *False* the result is the union of the evolved contour and the
        original seed (expansion-only mode).
    speed_image : sitk.Image, optional
        Pre-computed speed / feature image.  If *None* one is derived
        from *image* automatically.

    Returns
    -------
    sitk.Image
        Refined binary ROI.
    """
    image = _ensure_same_grid(image, roi, sitk.sitkLinear)

    seed = sitk.Cast(roi > 0, sitk.sitkUInt8)

    # Check for empty seed
    stats = sitk.StatisticsImageFilter()
    stats.Execute(seed)
    if stats.GetSum() == 0:
        return seed

    # Signed distance map (positive inside)
    distance_map = sitk.SignedMaurerDistanceMap(
        seed, insideIsPositive=True, squaredDistance=False, useImageSpacing=True
    )

    # Speed image
    if speed_image is not None:
        speed_image = _ensure_same_grid(speed_image, roi, sitk.sitkLinear)
        speed = sitk.Cast(speed_image, sitk.sitkFloat32)
    else:
        smoothed = sitk.CurvatureAnisotropicDiffusion(
            sitk.Cast(image, sitk.sitkFloat32),
            timeStep=0.01, conductanceParameter=3.0, numberOfIterations=5
        )
        grad_mag = sitk.GradientMagnitudeRecursiveGaussian(smoothed, sigma=1.0)
        minmax = sitk.MinimumMaximumImageFilter()
        minmax.Execute(grad_mag)
        if minmax.GetMaximum() > 0:
            grad_norm = sitk.Cast(grad_mag, sitk.sitkFloat32) / minmax.GetMaximum()
        else:
            grad_norm = sitk.Cast(grad_mag, sitk.sitkFloat32)
        speed = sitk.Exp(sitk.Multiply(grad_norm, -5.0))

    try:
        gac = sitk.GeodesicActiveContourLevelSetImageFilter()
        gac.SetPropagationScaling(propagation)
        gac.SetCurvatureScaling(curvature)
        gac.SetAdvectionScaling(advection)
        gac.SetMaximumRMSError(rms_tolerance)
        gac.SetNumberOfIterations(iterations)

        level_set = gac.Execute(
            sitk.Cast(distance_map, sitk.sitkFloat32),
            sitk.Cast(speed, sitk.sitkFloat32),
        )

        refined = sitk.BinaryThreshold(level_set, 0, 1000, 1, 0)
        refined = sitk.Cast(refined, sitk.sitkUInt8)

        if not allow_shrink:
            refined = sitk.Or(refined, seed)

        refined = sitk.BinaryMorphologicalClosing(refined, [1, 1, 1])
        return sitk.Cast(refined, sitk.sitkUInt8)

    except Exception:
        # Fallback: return seed unchanged
        return seed


# ============================================================================
# PROBABILITY-GATED EXPANSION
# ============================================================================

def probability_expansion(
    roi: sitk.Image,
    prob_map: sitk.Image,
    threshold: float = 0.25,
    max_layers: int = 5,
    max_distance_mm: float = None,
    fill_holes: bool = True,
    min_voxels: int = 50,
) -> sitk.Image:
    """
    Expand a binary ROI layer-by-layer, gated by a probability map.

    At each layer the ROI surface is dilated by one voxel and only voxels
    where ``prob_map > threshold`` are accepted.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.
    prob_map : sitk.Image
        Probability map (float [0, 1]).
    threshold : float
        Minimum probability for accepting a new voxel.
    max_layers : int
        Maximum expansion layers.
    max_distance_mm : float, optional
        Hard distance cap from original ROI surface.
    fill_holes : bool
        Fill holes after expansion.
    min_voxels : int
        Remove small components.

    Returns
    -------
    sitk.Image
        Expanded binary ROI.
    """
    prob_map = _ensure_same_grid(prob_map, roi, sitk.sitkLinear)
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    prob_arr = sitk.GetArrayFromImage(prob_map)
    spacing = _spacing_zyx(roi)

    # Distance constraint
    if max_distance_mm and max_distance_mm > 0:
        dist = distance_transform_edt(~roi_arr, sampling=spacing)
        dist_mask = dist <= max_distance_mm
    else:
        dist_mask = np.ones_like(roi_arr, dtype=bool)

    candidates = (prob_arr > threshold) & dist_mask

    struct = np.ones((3, 3, 3), dtype=bool)
    current = roi_arr.copy()

    for _ in range(max_layers):
        dilated = binary_dilation(current, structure=struct, iterations=1)
        layer = dilated & ~current & candidates
        if not np.any(layer):
            break
        current |= layer

    if fill_holes:
        current = binary_fill_holes(current)

    result = _arr_to_sitk_binary(current.astype(np.uint8), roi)
    result = filter_small_components(result, min_voxels=min_voxels)
    return result


# ============================================================================
# PROBABILITY-GATED SHRINKAGE
# ============================================================================

def probability_shrinkage(
    roi: sitk.Image,
    prob_map: sitk.Image,
    threshold: float = 0.15,
    max_layers: int = 5,
    min_preserve_fraction: float = 0.5,
    fill_holes: bool = True,
) -> sitk.Image:
    """
    Shrink a binary ROI by removing low-probability boundary voxels.

    At each layer the surface voxels are identified, and those with
    ``prob_map < threshold`` are removed.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.
    prob_map : sitk.Image
        Probability map (float [0, 1]).
    threshold : float
        Boundary voxels below this probability are removed.
    max_layers : int
        Maximum shrinkage layers.
    min_preserve_fraction : float
        Stop if the ROI would shrink below this fraction of its
        original volume.
    fill_holes : bool
        Fill holes after shrinkage.

    Returns
    -------
    sitk.Image
        Shrunk binary ROI.
    """
    prob_map = _ensure_same_grid(prob_map, roi, sitk.sitkLinear)
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    prob_arr = sitk.GetArrayFromImage(prob_map)

    original_voxels = int(np.sum(roi_arr))
    if original_voxels == 0:
        return roi
    min_voxels = int(original_voxels * min_preserve_fraction)

    struct = np.ones((3, 3, 3), dtype=bool)
    current = roi_arr.copy()

    for _ in range(max_layers):
        n_current = int(np.sum(current))
        if n_current <= min_voxels:
            break

        eroded = binary_erosion(current, structure=struct, iterations=1)
        boundary = current & ~eroded
        low_prob = boundary & (prob_arr < threshold)
        n_remove = int(np.sum(low_prob))

        if n_remove == 0:
            break

        # Respect volume floor
        if n_current - n_remove < min_voxels:
            # Remove only as many as we can, prioritising lowest prob
            n_allowed = n_current - min_voxels
            if n_allowed <= 0:
                break
            probs = prob_arr[low_prob]
            cutoff = np.partition(probs, n_allowed)[n_allowed]
            low_prob = low_prob & (prob_arr < cutoff)

        current = current & ~low_prob

    if fill_holes:
        current = binary_fill_holes(current)

    return _arr_to_sitk_binary(current.astype(np.uint8), roi)


# ============================================================================
# STAPLE FUSION
# ============================================================================

def staple_fusion(
    segmentations,
    reference: sitk.Image = None,
    confidence_threshold: float = 0.5,
    min_voxels: int = 50,
) -> sitk.Image:
    """
    Fuse multiple binary segmentations using the STAPLE algorithm.

    STAPLE (Simultaneous Truth and Performance Level Estimation) computes
    a probabilistic estimate of the true segmentation from a collection of
    imperfect segmentations.

    Parameters
    ----------
    segmentations : list
        List of binary sitk.Image or Imaginable-like objects.
    reference : sitk.Image, optional
        Reference geometry.  Defaults to the first segmentation.
    confidence_threshold : float
        Threshold on STAPLE probability for the final binary output.
    min_voxels : int
        Remove small components from the fused result.

    Returns
    -------
    sitk.Image
        Fused binary ROI.
    """
    segs = [_to_sitk(s) for s in segmentations]
    if len(segs) < 2:
        return segs[0] if segs else None

    if reference is None:
        reference = segs[0]

    resampled = []
    for s in segs:
        s = _ensure_same_grid(s, reference, sitk.sitkNearestNeighbor)
        resampled.append(sitk.Cast(s > 0, sitk.sitkUInt8))

    sf = sitk.STAPLEImageFilter()
    sf.SetForegroundValue(1)
    prob = sf.Execute(resampled)

    fused = sitk.BinaryThreshold(prob, confidence_threshold, 1.0, 1, 0)
    fused = sitk.Cast(fused, sitk.sitkUInt8)

    fused_arr = binary_fill_holes(
        sitk.GetArrayFromImage(fused) > 0
    ).astype(np.uint8)
    result = _arr_to_sitk_binary(fused_arr, reference)
    result = filter_small_components(result, min_voxels=min_voxels)
    return result


# ============================================================================
# LABEL OVERLAP RESOLUTION
# ============================================================================

def resolve_label_overlaps(
    label_rois: dict,
    reference: sitk.Image = None,
) -> dict:
    """
    Resolve overlapping voxels between multiple label ROIs.

    Overlapping voxels are assigned to the label whose centroid is closest
    (Euclidean distance in physical mm).

    Parameters
    ----------
    label_rois : dict
        ``{label_value: sitk.Image}`` binary ROI per label.
    reference : sitk.Image, optional
        Reference geometry (defaults to first ROI).

    Returns
    -------
    dict
        ``{label_value: sitk.Image}`` with overlaps resolved.
    """
    if len(label_rois) < 2:
        return label_rois

    labels = sorted(label_rois.keys())
    ref = reference or label_rois[labels[0]]
    spacing = _spacing_zyx(ref)

    # Gather arrays
    arrays = {}
    for lbl in labels:
        img = _ensure_same_grid(label_rois[lbl], ref, sitk.sitkNearestNeighbor)
        arrays[lbl] = sitk.GetArrayFromImage(img) > 0

    # Compute distance from each label (EDT of complement)
    distances = {}
    for lbl in labels:
        distances[lbl] = distance_transform_edt(~arrays[lbl], sampling=spacing)

    # Find overlapping voxels (belonging to 2+ labels)
    overlap_count = np.zeros_like(arrays[labels[0]], dtype=np.int32)
    for lbl in labels:
        overlap_count += arrays[lbl].astype(np.int32)

    overlap_mask = overlap_count >= 2

    if not np.any(overlap_mask):
        return label_rois

    # Assign overlapping voxels to the nearest label
    for lbl in labels:
        others_closer = np.zeros_like(overlap_mask, dtype=bool)
        for other in labels:
            if other == lbl:
                continue
            others_closer |= (distances[other] < distances[lbl])
        # Remove voxels where another label is closer
        arrays[lbl][overlap_mask & others_closer] = False

    # Rebuild sitk images
    result = {}
    for lbl in labels:
        result[lbl] = _arr_to_sitk_binary(arrays[lbl], ref)
    return result


# ============================================================================
# THRESHOLD LEVEL-SET REFINEMENT
# ============================================================================

def threshold_level_set_refine(
    roi: sitk.Image,
    image: sitk.Image,
    lower_threshold: float = 0.1,
    upper_threshold: float = 0.9,
    propagation: float = 1.0,
    curvature: float = 1.0,
    iterations: int = 100,
    rms_tolerance: float = 0.02,
    allow_shrink: bool = True,
) -> sitk.Image:
    """
    Refine a binary ROI using a threshold-based level-set.

    The contour expands into voxels whose intensity falls within
    [lower_threshold, upper_threshold] (after normalisation to [0, 1]).
    This is useful when the target tissue has a well-defined intensity
    range but irregular boundaries.

    Parameters
    ----------
    roi : sitk.Image
        Binary seed ROI.
    image : sitk.Image
        Scalar intensity image.
    lower_threshold : float
        Lower intensity bound (normalised [0, 1]).
    upper_threshold : float
        Upper intensity bound (normalised [0, 1]).
    propagation : float
        Balloon force.
    curvature : float
        Smoothing force.
    iterations : int
        Maximum iterations.
    rms_tolerance : float
        Convergence threshold.
    allow_shrink : bool
        If False, result is unioned with the seed.

    Returns
    -------
    sitk.Image
        Refined binary ROI (UInt8).
    """
    image = _ensure_same_grid(image, roi, sitk.sitkLinear)
    seed = sitk.Cast(roi > 0, sitk.sitkUInt8)

    stats = sitk.StatisticsImageFilter()
    stats.Execute(seed)
    if stats.GetSum() == 0:
        return seed

    # Normalise image intensity to [0, 1]
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    mmf = sitk.MinimumMaximumImageFilter()
    mmf.Execute(img_f)
    lo, hi = mmf.GetMinimum(), mmf.GetMaximum()
    if hi - lo > 0:
        img_norm = (img_f - lo) / (hi - lo)
    else:
        img_norm = img_f

    # Signed distance initialisation
    dist = sitk.SignedMaurerDistanceMap(
        seed, insideIsPositive=True, squaredDistance=False, useImageSpacing=True,
    )

    try:
        tls = sitk.ThresholdSegmentationLevelSetImageFilter()
        tls.SetLowerThreshold(lower_threshold)
        tls.SetUpperThreshold(upper_threshold)
        tls.SetPropagationScaling(propagation)
        tls.SetCurvatureScaling(curvature)
        tls.SetMaximumRMSError(rms_tolerance)
        tls.SetNumberOfIterations(iterations)

        ls_out = tls.Execute(
            sitk.Cast(dist, sitk.sitkFloat32),
            sitk.Cast(img_norm, sitk.sitkFloat32),
        )

        refined = sitk.BinaryThreshold(ls_out, 0, 1e10, 1, 0)
        refined = sitk.Cast(refined, sitk.sitkUInt8)

        if not allow_shrink:
            refined = sitk.Or(refined, seed)

        return refined
    except Exception:
        return seed


# ============================================================================
# LAPLACIAN LEVEL-SET REFINEMENT
# ============================================================================

def laplacian_level_set_refine(
    roi: sitk.Image,
    image: sitk.Image,
    propagation: float = 1.0,
    curvature: float = 1.0,
    iterations: int = 100,
    rms_tolerance: float = 0.02,
    allow_shrink: bool = True,
) -> sitk.Image:
    """
    Refine a binary ROI using a Laplacian-based level-set.

    The speed function is the Laplacian of the image, which drives
    the contour towards intensity edges (zero-crossings).

    Parameters
    ----------
    roi : sitk.Image
        Binary seed ROI.
    image : sitk.Image
        Scalar intensity image.
    propagation : float
        Balloon force.
    curvature : float
        Smoothing force.
    iterations : int
        Maximum iterations.
    rms_tolerance : float
        Convergence threshold.
    allow_shrink : bool
        If False, result is unioned with the seed.

    Returns
    -------
    sitk.Image
        Refined binary ROI (UInt8).
    """
    image = _ensure_same_grid(image, roi, sitk.sitkLinear)
    seed = sitk.Cast(roi > 0, sitk.sitkUInt8)

    stats = sitk.StatisticsImageFilter()
    stats.Execute(seed)
    if stats.GetSum() == 0:
        return seed

    img_f = sitk.Cast(image, sitk.sitkFloat32)
    dist = sitk.SignedMaurerDistanceMap(
        seed, insideIsPositive=True, squaredDistance=False, useImageSpacing=True,
    )

    try:
        lls = sitk.LaplacianSegmentationLevelSetImageFilter()
        lls.SetPropagationScaling(propagation)
        lls.SetCurvatureScaling(curvature)
        lls.SetMaximumRMSError(rms_tolerance)
        lls.SetNumberOfIterations(iterations)

        ls_out = lls.Execute(
            sitk.Cast(dist, sitk.sitkFloat32),
            sitk.Cast(img_f, sitk.sitkFloat32),
        )

        refined = sitk.BinaryThreshold(ls_out, 0, 1e10, 1, 0)
        refined = sitk.Cast(refined, sitk.sitkUInt8)

        if not allow_shrink:
            refined = sitk.Or(refined, seed)

        return refined
    except Exception:
        return seed


# ============================================================================
# SHAPE DETECTION LEVEL-SET
# ============================================================================

def shape_detection_level_set_refine(
    roi: sitk.Image,
    image: sitk.Image,
    propagation: float = 1.0,
    curvature: float = 0.5,
    iterations: int = 100,
    rms_tolerance: float = 0.02,
    sigma_mm: float = 1.0,
    allow_shrink: bool = True,
) -> sitk.Image:
    """
    Refine ROI using a shape detection level-set.

    Similar to geodesic active contour but without advection.
    The speed is derived from the edge potential of the image.

    Parameters
    ----------
    roi : sitk.Image
        Binary seed ROI.
    image : sitk.Image
        Scalar intensity image.
    propagation : float
        Balloon force.
    curvature : float
        Smoothing force.
    iterations : int
        Maximum iterations.
    rms_tolerance : float
        Convergence threshold.
    sigma_mm : float
        Gaussian sigma for edge computation.
    allow_shrink : bool
        If False, result is unioned with the seed.

    Returns
    -------
    sitk.Image
        Refined binary ROI (UInt8).
    """
    image = _ensure_same_grid(image, roi, sitk.sitkLinear)
    seed = sitk.Cast(roi > 0, sitk.sitkUInt8)

    stats = sitk.StatisticsImageFilter()
    stats.Execute(seed)
    if stats.GetSum() == 0:
        return seed

    img_f = sitk.Cast(image, sitk.sitkFloat32)
    dist = sitk.SignedMaurerDistanceMap(
        seed, insideIsPositive=True, squaredDistance=False, useImageSpacing=True,
    )

    # Edge potential (speed image)
    grad_mag = sitk.GradientMagnitudeRecursiveGaussian(img_f, sigma=sigma_mm)
    speed = sitk.BoundedReciprocal(grad_mag)

    try:
        sdls = sitk.ShapeDetectionLevelSetImageFilter()
        sdls.SetPropagationScaling(propagation)
        sdls.SetCurvatureScaling(curvature)
        sdls.SetMaximumRMSError(rms_tolerance)
        sdls.SetNumberOfIterations(iterations)

        ls_out = sdls.Execute(
            sitk.Cast(dist, sitk.sitkFloat32),
            sitk.Cast(speed, sitk.sitkFloat32),
        )

        refined = sitk.BinaryThreshold(ls_out, 0, 1e10, 1, 0)
        refined = sitk.Cast(refined, sitk.sitkUInt8)

        if not allow_shrink:
            refined = sitk.Or(refined, seed)

        return refined
    except Exception:
        return seed


# ============================================================================
# CHAN-VESE (REGION-BASED) SEGMENTATION
# ============================================================================

def chan_vese_refine(
    roi: sitk.Image,
    image: sitk.Image,
    lambda1: float = 1.0,
    lambda2: float = 1.0,
    curvature_weight: float = 0.0,
    area_weight: float = 0.0,
    volume_weight: float = 0.0,
    iterations: int = 100,
    rms_tolerance: float = 0.02,
    allow_shrink: bool = True,
) -> sitk.Image:
    """
    Refine a binary ROI using the Chan-Vese (region-based) level-set.

    This method does not rely on edge information and works well for
    images with weak or absent edges.  It separates the image into
    foreground and background based on intensity homogeneity.

    Parameters
    ----------
    roi : sitk.Image
        Binary seed ROI.
    image : sitk.Image
        Scalar intensity image.
    lambda1 : float
        Weight for inside-region variance penalty.
    lambda2 : float
        Weight for outside-region variance penalty.
    curvature_weight : float
        Curvature regularisation.
    area_weight : float
        Area penalty weight.
    volume_weight : float
        Volume penalty weight.
    iterations : int
        Maximum iterations.
    rms_tolerance : float
        Convergence threshold.
    allow_shrink : bool
        If False, result is unioned with the seed.

    Returns
    -------
    sitk.Image
        Refined binary ROI (UInt8).
    """
    image = _ensure_same_grid(image, roi, sitk.sitkLinear)
    seed = sitk.Cast(roi > 0, sitk.sitkUInt8)

    stats = sitk.StatisticsImageFilter()
    stats.Execute(seed)
    if stats.GetSum() == 0:
        return seed

    img_f = sitk.Cast(image, sitk.sitkFloat32)

    try:
        cv = sitk.ScalarChanAndVeseDenseLevelSetImageFilter()
        cv.SetLambda1(lambda1)
        cv.SetLambda2(lambda2)
        cv.SetCurvatureWeight(curvature_weight)
        cv.SetAreaWeight(area_weight)
        cv.SetVolumeWeight(volume_weight)
        cv.SetMaximumRMSError(rms_tolerance)
        cv.SetNumberOfIterations(iterations)

        dist = sitk.SignedMaurerDistanceMap(
            seed, insideIsPositive=True, squaredDistance=False, useImageSpacing=True,
        )

        ls_out = cv.Execute(
            sitk.Cast(dist, sitk.sitkFloat32),
            sitk.Cast(img_f, sitk.sitkFloat32),
        )

        refined = sitk.BinaryThreshold(ls_out, 0, 1e10, 1, 0)
        refined = sitk.Cast(refined, sitk.sitkUInt8)

        if not allow_shrink:
            refined = sitk.Or(refined, seed)

        return refined
    except Exception:
        return seed


# ============================================================================
# THRESHOLDING METHODS (Otsu, Multi-Otsu, Li, Yen, Triangle, Huang)
# ============================================================================

def otsu_threshold(image: sitk.Image, n_bins: int = 128) -> sitk.Image:
    """
    Segment an image using Otsu's automatic threshold.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    n_bins : int
        Number of histogram bins for threshold computation.

    Returns
    -------
    sitk.Image
        Binary segmentation (UInt8).
    """
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    filt = sitk.OtsuThresholdImageFilter()
    filt.SetInsideValue(0)
    filt.SetOutsideValue(1)
    filt.SetNumberOfHistogramBins(n_bins)
    return sitk.Cast(filt.Execute(img_f), sitk.sitkUInt8)


def multi_otsu_threshold(
    image: sitk.Image,
    n_thresholds: int = 2,
    n_bins: int = 256,
) -> sitk.Image:
    """
    Segment an image using multi-level Otsu thresholding.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    n_thresholds : int
        Number of thresholds (produces n_thresholds + 1 classes).
    n_bins : int
        Number of histogram bins.

    Returns
    -------
    sitk.Image
        Multi-label segmentation (UInt8).
    """
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    filt = sitk.OtsuMultipleThresholdsImageFilter()
    filt.SetNumberOfThresholds(n_thresholds)
    filt.SetNumberOfHistogramBins(n_bins)
    return sitk.Cast(filt.Execute(img_f), sitk.sitkUInt8)


def li_threshold(image: sitk.Image, n_bins: int = 128) -> sitk.Image:
    """
    Segment an image using Li's iterative minimum cross-entropy threshold.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    n_bins : int
        Number of histogram bins for threshold computation.

    Returns
    -------
    sitk.Image
        Binary segmentation (UInt8).
    """
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    filt = sitk.LiThresholdImageFilter()
    filt.SetInsideValue(0)
    filt.SetOutsideValue(1)
    filt.SetNumberOfHistogramBins(n_bins)
    return sitk.Cast(filt.Execute(img_f), sitk.sitkUInt8)


def yen_threshold(image: sitk.Image, n_bins: int = 128) -> sitk.Image:
    """
    Segment an image using Yen's entropy-based threshold.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    n_bins : int
        Number of histogram bins for threshold computation.

    Returns
    -------
    sitk.Image
        Binary segmentation (UInt8).
    """
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    filt = sitk.YenThresholdImageFilter()
    filt.SetInsideValue(0)
    filt.SetOutsideValue(1)
    filt.SetNumberOfHistogramBins(n_bins)
    return sitk.Cast(filt.Execute(img_f), sitk.sitkUInt8)


def triangle_threshold(image: sitk.Image, n_bins: int = 128) -> sitk.Image:
    """
    Segment an image using the triangle (Zack) threshold method.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    n_bins : int
        Number of histogram bins for threshold computation.

    Returns
    -------
    sitk.Image
        Binary segmentation (UInt8).
    """
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    filt = sitk.TriangleThresholdImageFilter()
    filt.SetInsideValue(0)
    filt.SetOutsideValue(1)
    filt.SetNumberOfHistogramBins(n_bins)
    return sitk.Cast(filt.Execute(img_f), sitk.sitkUInt8)


def huang_threshold(image: sitk.Image, n_bins: int = 128) -> sitk.Image:
    """
    Segment an image using Huang's fuzzy threshold method.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    n_bins : int
        Number of histogram bins for threshold computation.

    Returns
    -------
    sitk.Image
        Binary segmentation (UInt8).
    """
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    filt = sitk.HuangThresholdImageFilter()
    filt.SetInsideValue(0)
    filt.SetOutsideValue(1)
    filt.SetNumberOfHistogramBins(n_bins)
    return sitk.Cast(filt.Execute(img_f), sitk.sitkUInt8)


def manual_threshold(
    image: sitk.Image,
    lower: float = 0.0,
    upper: float = 1.0,
) -> sitk.Image:
    """
    Segment an image by applying a manual intensity threshold.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    lower : float
        Lower intensity bound (inclusive).
    upper : float
        Upper intensity bound (inclusive).

    Returns
    -------
    sitk.Image
        Binary segmentation (UInt8).
    """
    return sitk.Cast(
        sitk.BinaryThreshold(sitk.Cast(image, sitk.sitkFloat32), lower, upper, 1, 0),
        sitk.sitkUInt8,
    )


# ============================================================================
# CONNECTED THRESHOLD REGION GROWING
# ============================================================================

def connected_threshold_grow(
    image: sitk.Image,
    seed_roi: sitk.Image,
    lower: float = None,
    upper: float = None,
    n_seeds: int = 200,
    replace_value: int = 1,
    face_connected: bool = True,
) -> sitk.Image:
    """
    Region growing from seed points with explicit intensity bounds.

    Seeds are sampled from the seed ROI.  Unlike confidence-connected,
    the intensity bounds are specified directly.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    seed_roi : sitk.Image
        Binary seed ROI — points are sampled from its interior.
    lower : float, optional
        Lower intensity bound.  If None, computed as mean - 2*std of
        intensities inside the seed ROI.
    upper : float, optional
        Upper intensity bound.  If None, computed as mean + 2*std.
    n_seeds : int
        Maximum number of seed points.
    replace_value : int
        Value for the grown region.
    face_connected : bool
        If True, use 6-connectivity (face). If False, use 26-connectivity
        (full) which allows diagonal growth.

    Returns
    -------
    sitk.Image
        Binary segmentation (UInt8).
    """
    image = _ensure_same_grid(image, seed_roi, sitk.sitkLinear)
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    roi_arr = sitk.GetArrayFromImage(seed_roi) > 0

    seed_indices = np.argwhere(roi_arr)
    if len(seed_indices) == 0:
        return sitk.Cast(seed_roi, sitk.sitkUInt8)

    # Auto-compute bounds from seed region
    if lower is None or upper is None:
        img_arr = sitk.GetArrayFromImage(img_f)
        vals = img_arr[roi_arr]
        mean_v, std_v = float(vals.mean()), float(vals.std())
        if lower is None:
            lower = mean_v - 2.0 * std_v
        if upper is None:
            upper = mean_v + 2.0 * std_v

    # Sub-sample seeds
    if len(seed_indices) > n_seeds:
        idx = np.linspace(0, len(seed_indices) - 1, n_seeds, dtype=int)
        seed_indices = seed_indices[idx]

    seeds = [tuple(int(v) for v in reversed(s)) for s in seed_indices]

    ct = sitk.ConnectedThresholdImageFilter()
    ct.SetLower(float(lower))
    ct.SetUpper(float(upper))
    ct.SetReplaceValue(replace_value)
    if face_connected:
        ct.SetConnectivity(0)  # Face connectivity (6-connected)
    else:
        ct.SetConnectivity(1)  # Full connectivity (26-connected)
    for s in seeds:
        ct.AddSeed(s)

    result = ct.Execute(img_f)
    return sitk.Cast(result, sitk.sitkUInt8)


# ============================================================================
# NEIGHBOURHOOD CONNECTED REGION GROWING
# ============================================================================

def neighbourhood_connected_grow(
    image: sitk.Image,
    seed_roi: sitk.Image,
    lower: float = None,
    upper: float = None,
    radius: int = 1,
    n_seeds: int = 200,
) -> sitk.Image:
    """
    Region growing with neighbourhood connectivity constraint.

    Similar to connected threshold but also checks the neighbourhood
    of each candidate voxel.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    seed_roi : sitk.Image
        Binary seed ROI.
    lower : float, optional
        Lower intensity bound (auto-computed if None).
    upper : float, optional
        Upper intensity bound (auto-computed if None).
    radius : int
        Neighbourhood radius.
    n_seeds : int
        Maximum number of seed points.

    Returns
    -------
    sitk.Image
        Binary segmentation (UInt8).
    """
    image = _ensure_same_grid(image, seed_roi, sitk.sitkLinear)
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    roi_arr = sitk.GetArrayFromImage(seed_roi) > 0

    seed_indices = np.argwhere(roi_arr)
    if len(seed_indices) == 0:
        return sitk.Cast(seed_roi, sitk.sitkUInt8)

    if lower is None or upper is None:
        img_arr = sitk.GetArrayFromImage(img_f)
        vals = img_arr[roi_arr]
        mean_v, std_v = float(vals.mean()), float(vals.std())
        if lower is None:
            lower = mean_v - 2.0 * std_v
        if upper is None:
            upper = mean_v + 2.0 * std_v

    if len(seed_indices) > n_seeds:
        idx = np.linspace(0, len(seed_indices) - 1, n_seeds, dtype=int)
        seed_indices = seed_indices[idx]

    seeds = [tuple(int(v) for v in reversed(s)) for s in seed_indices]

    nc = sitk.NeighborhoodConnectedImageFilter()
    nc.SetLower(float(lower))
    nc.SetUpper(float(upper))
    nc.SetRadius([radius] * 3)
    for s in seeds:
        nc.AddSeed(s)

    result = nc.Execute(img_f)
    return sitk.Cast(result, sitk.sitkUInt8)


# ============================================================================
# ISOLATED CONNECTED REGION GROWING
# ============================================================================

def isolated_connected_grow(
    image: sitk.Image,
    seed1_roi: sitk.Image,
    seed2_roi: sitk.Image,
    n_seeds: int = 50,
) -> sitk.Image:
    """
    Find the intensity threshold that separates two seed regions.

    Grows from seed1 while staying disconnected from seed2.
    Useful for separating two adjacent structures.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    seed1_roi : sitk.Image
        Binary ROI for the target region.
    seed2_roi : sitk.Image
        Binary ROI for the excluded region.
    n_seeds : int
        Maximum seeds per region.

    Returns
    -------
    sitk.Image
        Binary segmentation of seed1 region (UInt8).
    """
    image = _ensure_same_grid(image, seed1_roi, sitk.sitkLinear)
    seed2_roi = _ensure_same_grid(seed2_roi, seed1_roi, sitk.sitkNearestNeighbor)
    img_f = sitk.Cast(image, sitk.sitkFloat32)

    def _sample_seeds(roi_img, max_n):
        arr = sitk.GetArrayFromImage(roi_img) > 0
        struct = np.ones((3, 3, 3), dtype=bool)
        eroded = binary_erosion(arr, structure=struct, iterations=1)
        if not np.any(eroded):
            eroded = arr
        indices = np.argwhere(eroded)
        if len(indices) > max_n:
            idx = np.linspace(0, len(indices) - 1, max_n, dtype=int)
            indices = indices[idx]
        return [tuple(int(v) for v in reversed(s)) for s in indices]

    seeds1 = _sample_seeds(seed1_roi, n_seeds)
    seeds2 = _sample_seeds(seed2_roi, n_seeds)

    if not seeds1 or not seeds2:
        return sitk.Cast(seed1_roi > 0, sitk.sitkUInt8)

    ic = sitk.IsolatedConnectedImageFilter()
    ic.SetSeed1(seeds1[0])
    ic.SetSeed2(seeds2[0])

    result = ic.Execute(img_f)
    return sitk.Cast(result, sitk.sitkUInt8)


# ============================================================================
# MORPHOLOGICAL WATERSHED FROM MARKERS
# ============================================================================

def morphological_watershed(
    image: sitk.Image,
    level: float = 0.1,
    fully_connected: bool = False,
) -> sitk.Image:
    """
    Apply morphological watershed segmentation (no markers).

    Produces a label image of watershed basins from the gradient of
    the input image.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    level : float
        Flooding level — higher values produce fewer basins (more merged).
    fully_connected : bool
        Use 26-connectivity (True) vs 6-connectivity (False).

    Returns
    -------
    sitk.Image
        Label image of watershed basins (UInt16).
    """
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    grad = sitk.GradientMagnitudeRecursiveGaussian(img_f, sigma=1.0)

    ws = sitk.MorphologicalWatershedImageFilter()
    ws.SetLevel(level)
    ws.SetFullyConnected(fully_connected)
    ws.SetMarkWatershedLine(False)

    return sitk.Cast(ws.Execute(grad), sitk.sitkUInt16)


def morphological_watershed_from_markers(
    image: sitk.Image,
    markers: sitk.Image,
    fully_connected: bool = False,
) -> sitk.Image:
    """
    Watershed segmentation driven by user-provided markers.

    Each connected marker region becomes a separate basin.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image (gradient is computed internally).
    markers : sitk.Image
        Integer marker image (each label seeds a basin; 0 = unlabelled).
    fully_connected : bool
        Use 26-connectivity (True) vs 6-connectivity (False).

    Returns
    -------
    sitk.Image
        Label image of watershed basins (UInt16).
    """
    image = _ensure_same_grid(image, markers, sitk.sitkLinear)
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    grad = sitk.GradientMagnitudeRecursiveGaussian(img_f, sigma=1.0)

    ws = sitk.MorphologicalWatershedFromMarkersImageFilter()
    ws.SetFullyConnected(fully_connected)
    ws.SetMarkWatershedLine(False)

    return sitk.Cast(
        ws.Execute(grad, sitk.Cast(markers, sitk.sitkUInt32)),
        sitk.sitkUInt16,
    )


# ============================================================================
# DISTANCE-CONSTRAINED CLIPPING
# ============================================================================

def constrain_by_distance(
    roi: sitk.Image,
    reference_roi: sitk.Image,
    max_distance_mm: float = 5.0,
    exclude_interior: bool = False,
) -> sitk.Image:
    """
    Clip a ROI to stay within a maximum distance from a reference ROI.

    Useful for constraining cartilage near bone, or limiting leakage.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI to constrain.
    reference_roi : sitk.Image
        Binary reference ROI.
    max_distance_mm : float
        Maximum allowed distance from the reference surface.
    exclude_interior : bool
        If True, also remove voxels inside the reference ROI
        (e.g., cartilage should not overlap bone).

    Returns
    -------
    sitk.Image
        Constrained binary ROI (UInt8).
    """
    reference_roi = _ensure_same_grid(reference_roi, roi, sitk.sitkNearestNeighbor)
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    ref_arr = sitk.GetArrayFromImage(reference_roi) > 0
    spacing = _spacing_zyx(roi)

    dist = distance_transform_edt(~ref_arr, sampling=spacing)
    valid = dist <= max_distance_mm

    if exclude_interior:
        valid = valid & (~ref_arr)

    result = roi_arr & valid
    return _arr_to_sitk_binary(result, roi)


# ============================================================================
# MASK SUBTRACTION / INTERSECTION
# ============================================================================

def subtract_mask(
    roi: sitk.Image,
    mask_to_remove: sitk.Image,
) -> sitk.Image:
    """
    Remove voxels from roi that overlap with mask_to_remove.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.
    mask_to_remove : sitk.Image
        Binary mask of voxels to remove.

    Returns
    -------
    sitk.Image
        ROI with overlapping voxels removed (UInt8).
    """
    mask_to_remove = _ensure_same_grid(mask_to_remove, roi, sitk.sitkNearestNeighbor)
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    mask_arr = sitk.GetArrayFromImage(mask_to_remove) > 0
    result = roi_arr & (~mask_arr)
    return _arr_to_sitk_binary(result, roi)


def intersect_masks(
    roi: sitk.Image,
    other: sitk.Image,
) -> sitk.Image:
    """
    Keep only voxels present in both ROIs.

    Parameters
    ----------
    roi : sitk.Image
        First binary ROI.
    other : sitk.Image
        Second binary ROI.

    Returns
    -------
    sitk.Image
        Intersection of the two ROIs (UInt8).
    """
    other = _ensure_same_grid(other, roi, sitk.sitkNearestNeighbor)
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    other_arr = sitk.GetArrayFromImage(other) > 0
    result = roi_arr & other_arr
    return _arr_to_sitk_binary(result, roi)


def union_masks(
    roi: sitk.Image,
    other: sitk.Image,
) -> sitk.Image:
    """
    Combine two binary ROIs (logical OR).

    Parameters
    ----------
    roi : sitk.Image
        First binary ROI.
    other : sitk.Image
        Second binary ROI.

    Returns
    -------
    sitk.Image
        Union of the two ROIs (UInt8).
    """
    other = _ensure_same_grid(other, roi, sitk.sitkNearestNeighbor)
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    other_arr = sitk.GetArrayFromImage(other) > 0
    result = roi_arr | other_arr
    return _arr_to_sitk_binary(result, roi)


# ============================================================================
# PREPROCESSING: N4 BIAS FIELD CORRECTION
# ============================================================================

def n4_bias_field_correction(
    image: sitk.Image,
    mask: sitk.Image = None,
    shrink_factor: int = 4,
    n_iterations: list = None,
    convergence_threshold: float = 0.001,
    spline_order: int = 3,
) -> sitk.Image:
    """
    Apply N4 bias field correction to an intensity image.

    Corrects for low-frequency intensity non-uniformity (e.g., from RF
    coil inhomogeneity in MRI).

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    mask : sitk.Image, optional
        Binary mask defining the region to use for bias estimation.
        If None, Otsu thresholding is used.
    shrink_factor : int
        Downsample factor for speed (default 4).
    n_iterations : list, optional
        Iterations per fitting level. Length controls number of fitting
        levels (default [50, 50, 50, 50] = 4 levels).
    convergence_threshold : float
        Convergence threshold.
    spline_order : int
        Order of the B-spline used for bias field estimation (default 3).

    Returns
    -------
    sitk.Image
        Bias-corrected image (Float32).
    """
    if n_iterations is None:
        n_iterations = [50, 50, 50, 50]

    img_f = sitk.Cast(image, sitk.sitkFloat32)

    if mask is None:
        mask = sitk.OtsuThreshold(img_f, 0, 1, 200)
    else:
        mask = _ensure_same_grid(mask, image, sitk.sitkNearestNeighbor)
        mask = sitk.Cast(mask > 0, sitk.sitkUInt8)

    # Shrink for speed
    shrunk_img = sitk.Shrink(img_f, [shrink_factor] * img_f.GetDimension())
    shrunk_mask = sitk.Shrink(mask, [shrink_factor] * mask.GetDimension())

    corrector = sitk.N4BiasFieldCorrectionImageFilter()
    corrector.SetMaximumNumberOfIterations(n_iterations)
    corrector.SetConvergenceThreshold(convergence_threshold)
    corrector.SetSplineOrder(spline_order)

    corrected_shrunk = corrector.Execute(shrunk_img, shrunk_mask)

    # Get log bias field and resample to full resolution
    log_bias = corrector.GetLogBiasFieldAsImage(img_f)
    corrected = img_f / sitk.Exp(log_bias)

    return corrected


# ============================================================================
# PREPROCESSING: ANISOTROPIC DIFFUSION SMOOTHING
# ============================================================================

def anisotropic_diffusion(
    image: sitk.Image,
    iterations: int = 5,
    time_step: float = 0.0625,
    conductance: float = 3.0,
) -> sitk.Image:
    """
    Apply curvature anisotropic diffusion smoothing.

    Smooths the image while preserving edges, useful as a preprocessing
    step before segmentation.

    Parameters
    ----------
    image : sitk.Image
        Scalar intensity image.
    iterations : int
        Number of diffusion iterations.
    time_step : float
        Time step per iteration (stability requires small values).
    conductance : float
        Conductance parameter — higher values smooth more aggressively
        across edges.

    Returns
    -------
    sitk.Image
        Smoothed image (Float32).
    """
    img_f = sitk.Cast(image, sitk.sitkFloat32)
    return sitk.CurvatureAnisotropicDiffusion(
        img_f,
        timeStep=time_step,
        conductanceParameter=conductance,
        numberOfIterations=iterations,
    )


# ============================================================================
# BINARY MORPHOLOGICAL OPERATIONS
# ============================================================================

def binary_erode(
    roi: sitk.Image,
    radius_mm: float = 1.0,
) -> sitk.Image:
    """
    Erode a binary ROI by a radius specified in mm.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.
    radius_mm : float
        Erosion radius in mm.

    Returns
    -------
    sitk.Image
        Eroded binary ROI (UInt8).
    """
    spacing = _spacing_zyx(roi)
    radius_voxels = [max(1, int(round(radius_mm / s))) for s in reversed(spacing)]
    binary = sitk.Cast(roi > 0, sitk.sitkUInt8)
    return sitk.BinaryErode(binary, radius_voxels)


def binary_dilate(
    roi: sitk.Image,
    radius_mm: float = 1.0,
) -> sitk.Image:
    """
    Dilate a binary ROI by a radius specified in mm.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.
    radius_mm : float
        Dilation radius in mm.

    Returns
    -------
    sitk.Image
        Dilated binary ROI (UInt8).
    """
    spacing = _spacing_zyx(roi)
    radius_voxels = [max(1, int(round(radius_mm / s))) for s in reversed(spacing)]
    binary = sitk.Cast(roi > 0, sitk.sitkUInt8)
    return sitk.BinaryDilate(binary, radius_voxels)


def binary_open(
    roi: sitk.Image,
    radius_mm: float = 1.0,
) -> sitk.Image:
    """
    Morphological opening (erosion followed by dilation) in mm.

    Removes small protrusions and disconnected fragments.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.
    radius_mm : float
        Structuring element radius in mm.

    Returns
    -------
    sitk.Image
        Opened binary ROI (UInt8).
    """
    spacing = _spacing_zyx(roi)
    radius_voxels = [max(1, int(round(radius_mm / s))) for s in reversed(spacing)]
    binary = sitk.Cast(roi > 0, sitk.sitkUInt8)
    return sitk.BinaryMorphologicalOpening(binary, radius_voxels)


def binary_close(
    roi: sitk.Image,
    radius_mm: float = 1.0,
) -> sitk.Image:
    """
    Morphological closing (dilation followed by erosion) in mm.

    Fills small holes and gaps in boundaries.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.
    radius_mm : float
        Structuring element radius in mm.

    Returns
    -------
    sitk.Image
        Closed binary ROI (UInt8).
    """
    spacing = _spacing_zyx(roi)
    radius_voxels = [max(1, int(round(radius_mm / s))) for s in reversed(spacing)]
    binary = sitk.Cast(roi > 0, sitk.sitkUInt8)
    return sitk.BinaryMorphologicalClosing(binary, radius_voxels)
