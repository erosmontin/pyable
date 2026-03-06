"""
Segmentation Refinement Utilities for Pyable

Generic functions for refining binary ROI masks using various strategies:
  - Watershed segmentation
  - Confidence-connected region growing
  - Geodesic active contour (level-set) boundary smoothing
  - Probability-gated expansion / shrinkage
  - STAPLE multi-segmentation fusion
  - Morphological utilities (fill holes, filter components, distance map, shell)

All functions operate on SimpleITK images and return SimpleITK images.
The corresponding methods on Roiable / LabelMapable are thin wrappers.

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
