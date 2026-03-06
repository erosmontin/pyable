"""
Metrics Module for Pyable

Provides comparison metrics, surface distance computation, morphometric
analysis, and quality scores for ROIs and label maps.

All functions operate on sitk.Image or numpy arrays and are used by the
Roiable / LabelMapable class methods.

Functions
---------
compare_segmentations : Full overlap + surface distance comparison
compute_overlap_metrics : Dice, Jaccard, volume similarity, FP/FN errors
compute_surface_distances : Symmetric surface distance stats (mean, RMS, HD95)
export_surface_with_error : Export meshes with per-vertex error for ParaView
compute_edge_alignment_score : How well ROI boundary aligns with image edges
compute_compactness_score : Volume / surface area ratio
compute_connectivity_score : Number of connected components
roi_principal_extents : PCA-based oriented bounding box extents
roi_max_feret : Maximum Feret diameter via convex hull
roi_max_thickness_inscribed : Max inscribed-sphere thickness (2 * max EDT)
roi_mean_thickness : Mean local thickness via EDT
build_label_priors : Build soft probability priors from a label map
"""

import numpy as np
import SimpleITK as sitk
from typing import Union, Optional, Dict, Tuple, List

from scipy.ndimage import (
    binary_dilation,
    binary_erosion,
    distance_transform_edt,
    label as nd_label,
)
from scipy.ndimage import gaussian_gradient_magnitude


# ============================================================================
# OVERLAP METRICS
# ============================================================================


def _ensure_same_geometry(ref: sitk.Image, mov: sitk.Image) -> sitk.Image:
    """Resample *mov* onto *ref* grid if geometries differ."""
    if (ref.GetSize() != mov.GetSize()
            or ref.GetSpacing() != mov.GetSpacing()
            or ref.GetOrigin() != mov.GetOrigin()
            or ref.GetDirection() != mov.GetDirection()):
        mov = sitk.Resample(
            mov, ref, sitk.Transform(),
            sitk.sitkNearestNeighbor, 0,
            useNearestNeighborExtrapolator=True,
        )
    return mov


def compute_overlap_metrics(
    reference: sitk.Image,
    test: sitk.Image,
) -> Dict[str, float]:
    """
    Overlap metrics between two binary / label images.

    Parameters
    ----------
    reference, test : sitk.Image
        Binary masks (>0 = foreground) on the same or different grids.

    Returns
    -------
    dict with Dice, Jaccard, VolumeSimilarity, FalseDiscoveryRate,
    FalseNegativeError, FalsePositiveError.
    """
    test = _ensure_same_geometry(reference, test)
    ov = sitk.LabelOverlapMeasuresImageFilter()
    ov.Execute(reference, test)
    return {
        'Dice': ov.GetDiceCoefficient(),
        'Jaccard': ov.GetJaccardCoefficient(),
        'VolumeSimilarity': ov.GetVolumeSimilarity(),
        'FalseDiscoveryRate': ov.GetFalseDiscoveryRate(),
        'FalseNegativeError': ov.GetFalseNegativeError(),
        'FalsePositiveError': ov.GetFalsePositiveError(),
    }


def compute_surface_distances(
    reference: sitk.Image,
    test: sitk.Image,
) -> Dict[str, Optional[float]]:
    """
    Symmetric surface distance metrics.

    Parameters
    ----------
    reference, test : sitk.Image
        Binary masks.

    Returns
    -------
    dict with MeanSurfaceDist_mm, RMSDist_mm, HD95_mm, SurfaceDice_1mm.
    """
    test = _ensure_same_geometry(reference, test)

    gt_surf = sitk.LabelContour(reference)
    pred_surf = sitk.LabelContour(test)

    dm_gt = sitk.Abs(sitk.SignedMaurerDistanceMap(
        gt_surf, squaredDistance=False, useImageSpacing=True))
    dm_pred = sitk.Abs(sitk.SignedMaurerDistanceMap(
        pred_surf, squaredDistance=False, useImageSpacing=True))

    a2b = sitk.GetArrayFromImage(sitk.Mask(dm_gt, pred_surf)).ravel()
    b2a = sitk.GetArrayFromImage(sitk.Mask(dm_pred, gt_surf)).ravel()

    a2b = a2b[a2b > 0]
    b2a = b2a[b2a > 0]
    concat = np.concatenate([a2b, b2a]) if (a2b.size + b2a.size) > 0 else np.array([])

    if concat.size == 0:
        return {
            'MeanSurfaceDist_mm': None,
            'RMSDist_mm': None,
            'HD95_mm': None,
            'SurfaceDice_1mm': None,
        }

    return {
        'MeanSurfaceDist_mm': float(np.mean(concat)),
        'RMSDist_mm': float(np.sqrt(np.mean(concat ** 2))),
        'HD95_mm': float(np.percentile(concat, 95)),
        'SurfaceDice_1mm': float(np.mean(concat < 1.0)),
    }


def compare_segmentations(
    reference: sitk.Image,
    test: sitk.Image,
) -> Dict[str, Optional[float]]:
    """
    Full comparison: overlap + surface distance metrics in one call.

    Parameters
    ----------
    reference, test : sitk.Image
        Binary masks.

    Returns
    -------
    dict merging overlap and surface distance metrics.
    """
    overlap = compute_overlap_metrics(reference, test)
    surface = compute_surface_distances(reference, test)
    return {**overlap, **surface}


# ============================================================================
# SURFACE EXPORT (for ParaView / 3D visualisation)
# ============================================================================


def export_surface_with_error(
    reference_arr: np.ndarray,
    test_arr: np.ndarray,
    spacing: Tuple[float, float, float] = (1.0, 1.0, 1.0),
    out_dir: str = '/tmp/surface_error',
    highlight_percentile: int = 95,
) -> None:
    """
    Export reference, test, and error-coloured meshes as PLY + VTP files.

    Requires ``pyvista`` and ``scikit-image``.

    Parameters
    ----------
    reference_arr, test_arr : np.ndarray
        Binary volumes in ZYX order (or sitk.Image, auto-converted).
    spacing : tuple
        Voxel size (x, y, z) in mm.
    out_dir : str
        Output directory.
    highlight_percentile : int
        Percentile threshold for the high-error sub-mesh.
    """
    import os
    from skimage import measure
    from scipy.spatial import cKDTree

    try:
        import pyvista as pv
    except ImportError:
        raise ImportError("pyvista is required for surface export: pip install pyvista")

    os.makedirs(out_dir, exist_ok=True)

    # Accept sitk.Image transparently
    if isinstance(reference_arr, sitk.Image):
        reference_arr = sitk.GetArrayFromImage(reference_arr)
    if isinstance(test_arr, sitk.Image):
        test_arr = sitk.GetArrayFromImage(test_arr)

    ref = np.asarray(reference_arr)
    tst = np.asarray(test_arr)

    marching_spacing = spacing[::-1]  # ZYX for marching_cubes

    verts_ref, faces_ref, _, _ = measure.marching_cubes(ref, level=0.5, spacing=marching_spacing)
    verts_tst, faces_tst, _, _ = measure.marching_cubes(tst, level=0.5, spacing=marching_spacing)

    tree_ref = cKDTree(verts_ref)
    dists_tst_to_ref, _ = tree_ref.query(verts_tst)

    def _encode_faces(faces):
        return np.hstack([np.full((faces.shape[0], 1), 3), faces]).astype(np.int64)

    mesh_ref = pv.PolyData(verts_ref, _encode_faces(faces_ref))
    mesh_tst = pv.PolyData(verts_tst, _encode_faces(faces_tst))

    mesh_err = mesh_tst.copy()
    mesh_err.point_data['surface_error_mm'] = dists_tst_to_ref.astype(np.float32)
    mesh_err.active_scalars_name = 'surface_error_mm'

    face_errs = dists_tst_to_ref[faces_tst].mean(axis=1).astype(np.float32)
    mesh_err.cell_data['face_error_mm'] = face_errs

    for name, m in [('Ref_surface', mesh_ref), ('Test_surface', mesh_tst), ('Error_surface', mesh_err)]:
        m.save(os.path.join(out_dir, f'{name}.ply'))
        m.save(os.path.join(out_dir, f'{name}.vtp'))

    thresh = float(np.percentile(face_errs, highlight_percentile))
    high_idx = np.where(face_errs >= thresh)[0]
    if high_idx.size > 0:
        sub_faces = faces_tst[high_idx]
        high_mesh = pv.PolyData(verts_tst, _encode_faces(sub_faces))
        high_mesh.point_data['surface_error_mm'] = mesh_err.point_data['surface_error_mm']
        high_mesh.cell_data['face_error_mm'] = face_errs[high_idx]
        high_mesh.active_scalars_name = 'surface_error_mm'
        high_mesh.save(os.path.join(out_dir, f'Error_highlight_p{highlight_percentile}.ply'))
        high_mesh.save(os.path.join(out_dir, f'Error_highlight_p{highlight_percentile}.vtp'))


# ============================================================================
# QUALITY SCORES (no ground truth needed)
# ============================================================================


def compute_edge_alignment_score(
    roi: sitk.Image,
    image: sitk.Image,
    band_width: int = 2,
    sigma_mm: float = 1.0,
) -> float:
    """
    Mean normalised gradient magnitude along the ROI surface band.

    Higher score → boundary aligns better with image edges.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.
    image : sitk.Image
        Grayscale image (e.g. MRI).
    band_width : int
        Erosion/dilation iterations for the surface band.
    sigma_mm : float
        Gaussian sigma in mm for gradient computation.

    Returns
    -------
    float
        Edge alignment score in [0, 1].
    """
    roi_arr = sitk.GetArrayFromImage(roi) > 0

    # Resample image to ROI grid if needed
    if image.GetSize() != roi.GetSize() or image.GetSpacing() != roi.GetSpacing():
        image = sitk.Resample(image, roi, sitk.Transform(), sitk.sitkLinear, 0.0)

    img_arr = sitk.GetArrayFromImage(image).astype(np.float32)
    spacing = np.array(roi.GetSpacing()[::-1])  # ZYX

    sigma_vox = sigma_mm / spacing
    grad_mag = gaussian_gradient_magnitude(img_arr, sigma=sigma_vox)

    grad_max = np.percentile(grad_mag[grad_mag > 0], 99) if np.any(grad_mag > 0) else 1.0
    grad_norm = np.clip(grad_mag / (grad_max + 1e-6), 0, 1)

    struct = np.ones((3, 3, 3), dtype=bool)
    dilated = binary_dilation(roi_arr, structure=struct, iterations=band_width)
    eroded = binary_erosion(roi_arr, structure=struct, iterations=band_width)
    boundary = dilated ^ eroded

    if not np.any(boundary):
        return 0.0

    return float(np.mean(grad_norm[boundary]))


def compute_compactness_score(roi: sitk.Image) -> float:
    """
    Compactness score: volume / surface-area.

    Higher → more compact / smoother shape.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.

    Returns
    -------
    float
    """
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    spacing = np.array(roi.GetSpacing()[::-1])
    voxel_vol = np.prod(spacing)

    if not np.any(roi_arr):
        return 0.0

    volume = np.sum(roi_arr) * voxel_vol

    struct = np.ones((3, 3, 3), dtype=bool)
    eroded = binary_erosion(roi_arr, structure=struct, iterations=1)
    boundary_voxels = np.sum(roi_arr & ~eroded)
    avg_spacing = np.mean(spacing)
    surface_area = boundary_voxels * (avg_spacing ** 2)

    if surface_area < 1e-6:
        return 1.0

    return float(volume / surface_area)


def compute_connectivity_score(roi: sitk.Image) -> int:
    """
    Number of connected components in a binary ROI.

    Parameters
    ----------
    roi : sitk.Image
        Binary ROI.

    Returns
    -------
    int
        Component count (ideally 1).
    """
    roi_arr = sitk.GetArrayFromImage(roi) > 0
    if not np.any(roi_arr):
        return 0
    _, n = nd_label(roi_arr)
    return int(n)


# ============================================================================
# MORPHOMETRIC ANALYSIS
# ============================================================================


def roi_principal_extents(
    mask: np.ndarray,
    spacing: Tuple[float, ...] = (1.0, 1.0, 1.0),
) -> Dict[str, float]:
    """
    PCA-based oriented bounding box extents.

    Parameters
    ----------
    mask : np.ndarray
        Binary mask in ZYX order.
    spacing : tuple
        Voxel spacing (sz, sy, sx) in mm.

    Returns
    -------
    dict with length, width, thickness (descending), eigenvectors, centre.
    """
    pts_zyx = np.column_stack(np.nonzero(mask))
    if pts_zyx.size == 0:
        return {'length': 0.0, 'width': 0.0, 'thickness': 0.0,
                'eigenvectors': None, 'centre': None}

    spacing = np.asarray(spacing)
    pts_phys = pts_zyx * spacing
    centre = pts_phys.mean(axis=0)
    X = pts_phys - centre
    cov = np.cov(X, rowvar=False)
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]
    evecs = evecs[:, order]
    proj = X @ evecs
    extents = proj.max(axis=0) - proj.min(axis=0)

    return {
        'length': float(extents[0]),
        'width': float(extents[1]),
        'thickness': float(extents[2]),
        'eigenvectors': evecs,
        'centre': centre,
    }


def roi_max_feret(
    mask: np.ndarray,
    spacing: Tuple[float, ...] = (1.0, 1.0, 1.0),
) -> float:
    """
    Maximum Feret diameter: max pairwise distance among convex-hull vertices.

    Parameters
    ----------
    mask : np.ndarray
        Binary mask in ZYX order.
    spacing : tuple
        Voxel spacing (sz, sy, sx).

    Returns
    -------
    float  Max Feret diameter in mm.
    """
    from scipy.spatial import ConvexHull, distance as sp_distance

    pts_zyx = np.column_stack(np.nonzero(mask))
    if pts_zyx.size == 0:
        return 0.0

    pts_phys = pts_zyx * np.asarray(spacing)
    if pts_phys.shape[0] < 4:
        return float(sp_distance.pdist(pts_phys).max()) if pts_phys.shape[0] > 1 else 0.0

    hull = ConvexHull(pts_phys)
    hull_pts = pts_phys[hull.vertices]
    return float(sp_distance.pdist(hull_pts).max())


def roi_max_thickness_inscribed(
    mask: np.ndarray,
    spacing: Tuple[float, ...] = (1.0, 1.0, 1.0),
) -> float:
    """
    Max inscribed-sphere thickness: ``2 * max(EDT)`` inside the mask.

    Parameters
    ----------
    mask : np.ndarray
        Binary mask in ZYX order.
    spacing : tuple
        Voxel spacing (sz, sy, sx).

    Returns
    -------
    float  Thickness in mm.
    """
    if mask.sum() == 0:
        return 0.0
    edt = distance_transform_edt(mask, sampling=spacing)
    return float(2.0 * edt.max())


def roi_mean_thickness(
    mask: np.ndarray,
    spacing: Tuple[float, ...] = (1.0, 1.0, 1.0),
) -> Dict[str, float]:
    """
    Local thickness statistics via EDT.

    Each interior voxel's thickness = ``2 * distance_to_boundary``.

    Parameters
    ----------
    mask : np.ndarray
        Binary mask in ZYX order.
    spacing : tuple
        Voxel spacing (sz, sy, sx).

    Returns
    -------
    dict with mean, median, std, max, p25, p75  (all in mm).
    """
    nan_result = dict(mean=np.nan, median=np.nan, std=np.nan,
                      max=np.nan, p25=np.nan, p75=np.nan)
    if mask.sum() == 0:
        return nan_result

    edt = distance_transform_edt(mask.astype(bool), sampling=spacing)
    thickness_map = 2.0 * edt[mask.astype(bool)]

    return {
        'mean': float(np.mean(thickness_map)),
        'median': float(np.median(thickness_map)),
        'std': float(np.std(thickness_map)),
        'max': float(np.max(thickness_map)),
        'p25': float(np.percentile(thickness_map, 25)),
        'p75': float(np.percentile(thickness_map, 75)),
    }


# ============================================================================
# LABEL-MAP PRIORS
# ============================================================================


def build_label_priors(
    label_image: sitk.Image,
    classes: Optional[List[int]] = None,
    tau: float = 0.8,
    blur_sigma_mm: float = 0.6,
    eps: float = 1e-6,
) -> Tuple[sitk.Image, List[int]]:
    """
    Build soft probability priors from a multi-label segmentation.

    For each foreground class, a sigmoid of the signed distance map is computed,
    optionally blurred. Returns a vector image whose channels are
    [background, class_1, class_2, …].

    Parameters
    ----------
    label_image : sitk.Image
        Integer label image.
    classes : list[int], optional
        Label values to include. If None, auto-detected from image.
    tau : float
        Sigmoid steepness parameter (smaller → sharper boundaries).
    blur_sigma_mm : float
        Gaussian blur sigma in mm for anti-aliasing (0 = no blur).
    eps : float
        Numeric floor to avoid log(0).

    Returns
    -------
    prior_vector : sitk.Image
        Vector image with one channel per class (including background channel 0).
    class_list : list[int]
        Ordered list ``[0, c1, c2, …]`` matching channels.
    """
    lab_np = sitk.GetArrayFromImage(label_image)

    if classes is None:
        classes = sorted(int(c) for c in np.unique(lab_np) if c != 0)

    fg = [c for c in classes if c != 0]
    var = float(blur_sigma_mm) ** 2 if blur_sigma_mm > 0 else 0.0

    priors = []
    for c in fg:
        bin_c = sitk.Equal(label_image, int(c))
        sd = sitk.SignedMaurerDistanceMap(
            bin_c, insideIsPositive=False,
            squaredDistance=False, useImageSpacing=True,
        )
        if var > 0:
            sd = sitk.DiscreteGaussian(sd, variance=var)
        d = sitk.GetArrayFromImage(sd)
        p = 1.0 / (1.0 + np.exp(d / float(tau)))
        priors.append(p.astype(np.float32))

    if priors:
        P_fg = np.stack(priors, axis=-1)
    else:
        P_fg = np.zeros(lab_np.shape + (0,), np.float32)

    P_bg = np.clip(1.0 - P_fg.sum(axis=-1), 0.0, 1.0)[..., None]
    P = np.concatenate([P_bg, P_fg], axis=-1)
    P = np.clip(P, eps, 1.0)
    P /= (P.sum(axis=-1, keepdims=True) + 1e-8)

    comps = [sitk.GetImageFromArray(P[..., i]) for i in range(P.shape[-1])]
    for im in comps:
        im.CopyInformation(label_image)
    prior_vec = sitk.Compose(comps)

    return prior_vec, [0] + fg


def combine_binary_masks_to_label(
    mask_images: List[sitk.Image],
    priority_order: Optional[List[int]] = None,
) -> sitk.Image:
    """
    Combine multiple binary masks into a single multi-label image.

    Parameters
    ----------
    mask_images : list[sitk.Image]
        Binary masks (same geometry). >0 = foreground.
    priority_order : list[int], optional
        Overwrite order (indices into mask_images). Last wins.
        If None, uses list order.

    Returns
    -------
    sitk.Image
        Integer label image (labels 1..K).
    """
    ref = mask_images[0]
    lab = np.zeros(sitk.GetArrayFromImage(ref).shape, dtype=np.int16)
    order = priority_order if priority_order is not None else list(range(len(mask_images)))

    for idx in order:
        arr = sitk.GetArrayFromImage(mask_images[idx]) > 0
        lab[arr] = idx + 1

    lab_img = sitk.GetImageFromArray(lab)
    lab_img.CopyInformation(ref)
    return lab_img


# ============================================================================
# DISTANCE MAP & SURFACE UTILITIES
# ============================================================================


def surface_mask(binary_mask: sitk.Image) -> sitk.Image:
    """
    Extract 1-voxel-thick surface contour from a binary mask.

    Parameters
    ----------
    binary_mask : sitk.Image
        Binary mask (>0 = foreground).

    Returns
    -------
    sitk.Image
        Binary image where 1 = surface voxels.
    """
    return sitk.LabelContour(sitk.Cast(binary_mask > 0, sitk.sitkUInt8))


def signed_distance_map(binary_mask: sitk.Image) -> sitk.Image:
    """
    Signed Maurer distance map from a binary mask.

    Negative inside, positive outside, in mm (physical spacing).

    Parameters
    ----------
    binary_mask : sitk.Image
        Binary mask (>0 = foreground).

    Returns
    -------
    sitk.Image
        Float32 signed distance map.
    """
    return sitk.Cast(
        sitk.SignedMaurerDistanceMap(
            sitk.Cast(binary_mask > 0, sitk.sitkUInt8),
            insideIsPositive=False,
            squaredDistance=False,
            useImageSpacing=True,
        ),
        sitk.sitkFloat32,
    )


def surface_distance_map_image(
    source_mask: sitk.Image,
    target_mask: sitk.Image,
) -> sitk.Image:
    """
    Image where SOURCE surface voxels contain distance (mm) to TARGET mask.

    All non-surface voxels are 0.  Useful for per-voxel error visualisation.

    Parameters
    ----------
    source_mask, target_mask : sitk.Image
        Binary masks (same geometry).

    Returns
    -------
    sitk.Image
        Float32 distance map (non-zero only on source surface).
    """
    target_mask = _ensure_same_geometry(source_mask, target_mask)
    src_surf = surface_mask(source_mask)
    dist_to_target = sitk.Abs(sitk.SignedMaurerDistanceMap(
        sitk.Cast(target_mask > 0, sitk.sitkUInt8),
        insideIsPositive=False,
        squaredDistance=False,
        useImageSpacing=True,
    ))
    out = sitk.Mask(dist_to_target, src_surf, outsideValue=0.0)
    return sitk.Cast(out, sitk.sitkFloat32)


# ============================================================================
# CHANGE MAPS (longitudinal pre/post analysis)
# ============================================================================


def compute_change_maps(
    pre_mask: sitk.Image,
    post_mask: sitk.Image,
) -> Dict[str, sitk.Image]:
    """
    Boolean change maps between pre and post binary masks.

    Parameters
    ----------
    pre_mask, post_mask : sitk.Image
        Binary masks on the same grid (post should already be registered).

    Returns
    -------
    dict with keys ``'removed'``, ``'added'``, ``'changed'`` (XOR),
    each a UInt8 binary image.
    """
    post_mask = _ensure_same_geometry(pre_mask, post_mask)
    pre = sitk.Cast(pre_mask > 0, sitk.sitkUInt8)
    post = sitk.Cast(post_mask > 0, sitk.sitkUInt8)

    removed = pre & sitk.Not(post)      # present pre, missing post
    added = post & sitk.Not(pre)         # missing pre, present post
    changed = pre ^ post                 # either removed or added

    return {
        'removed': sitk.Cast(removed, sitk.sitkUInt8),
        'added': sitk.Cast(added, sitk.sitkUInt8),
        'changed': sitk.Cast(changed, sitk.sitkUInt8),
    }


# ============================================================================
# COMPONENT SPLITTING
# ============================================================================


def split_connected_components(
    binary_mask: sitk.Image,
    sort_by: str = 'size',
) -> List[sitk.Image]:
    """
    Split a binary mask into its connected components.

    Parameters
    ----------
    binary_mask : sitk.Image
        Binary mask (>0 = foreground).
    sort_by : str
        ``'size'`` (descending voxel count) or ``'none'``.

    Returns
    -------
    list[sitk.Image]
        List of binary masks, one per component.
    """
    cc = sitk.ConnectedComponent(sitk.Cast(binary_mask > 0, sitk.sitkUInt8))
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(cc)
    labels = list(stats.GetLabels())

    if sort_by == 'size':
        labels = sorted(labels, key=lambda l: stats.GetNumberOfPixels(l), reverse=True)

    masks = []
    for lab in labels:
        m = sitk.Cast(sitk.Equal(cc, lab), sitk.sitkUInt8)
        masks.append(m)
    return masks


def split_left_right(
    binary_mask: sitk.Image,
    axis: int = 0,
) -> Tuple[sitk.Image, sitk.Image]:
    """
    Split a binary mask into left and right parts using connected-component
    centroids.

    The two largest components are identified.  "Left" is the one with the
    smaller centroid coordinate along *axis* (default: X = left-right in LPS).

    Parameters
    ----------
    binary_mask : sitk.Image
        Binary mask expected to contain ≥2 components.
    axis : int
        Physical axis index (0=X, 1=Y, 2=Z).  Default 0 (X).

    Returns
    -------
    (left, right) : tuple[sitk.Image, sitk.Image]
        Two binary masks.

    Raises
    ------
    RuntimeError
        If fewer than 2 components are found.
    """
    cc = sitk.ConnectedComponent(sitk.Cast(binary_mask > 0, sitk.sitkUInt8))
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(cc)
    labels = sorted(stats.GetLabels(),
                    key=lambda l: stats.GetNumberOfPixels(l), reverse=True)

    if len(labels) < 2:
        raise RuntimeError(
            f"Expected ≥2 components for left/right split, found {len(labels)}"
        )

    lab_a, lab_b = labels[0], labels[1]
    ca = np.array(stats.GetCentroid(lab_a))
    cb = np.array(stats.GetCentroid(lab_b))

    mask_a = sitk.Cast(sitk.Equal(cc, lab_a), sitk.sitkUInt8)
    mask_b = sitk.Cast(sitk.Equal(cc, lab_b), sitk.sitkUInt8)

    if ca[axis] <= cb[axis]:
        return mask_a, mask_b   # a = left, b = right
    else:
        return mask_b, mask_a


# ============================================================================
# RIGID REGISTRATION (mask-to-mask)
# ============================================================================


def rigid_register_masks(
    fixed_mask: sitk.Image,
    moving_mask: sitk.Image,
    method: str = 'rigid',
    iterations: int = 200,
) -> Tuple[sitk.Image, sitk.Transform]:
    """
    Register a moving binary mask onto a fixed mask using signed distance maps.

    Parameters
    ----------
    fixed_mask, moving_mask : sitk.Image
        Binary masks.
    method : str
        ``'rigid'`` (Euler3D) or ``'affine'``.
    iterations : int
        Max optimiser iterations.

    Returns
    -------
    (resampled, transform) : tuple
        ``resampled`` is the moving mask in fixed space (nearest-neighbor).
        ``transform`` is the registration transform.
    """
    fixed_sd = signed_distance_map(fixed_mask)
    moving_sd = signed_distance_map(moving_mask)

    if method == 'rigid':
        tfm = sitk.Euler3DTransform()
    elif method == 'affine':
        tfm = sitk.AffineTransform(fixed_sd.GetDimension())
    else:
        raise ValueError(f"Unknown method '{method}', choose 'rigid' or 'affine'")

    initial = sitk.CenteredTransformInitializer(
        fixed_sd, moving_sd, tfm,
        sitk.CenteredTransformInitializerFilter.GEOMETRY,
    )

    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMeanSquares()
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsRegularStepGradientDescent(
        learningRate=1.0,
        minStep=1e-3,
        numberOfIterations=iterations,
        gradientMagnitudeTolerance=1e-6,
    )
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetInitialTransform(initial, inPlace=False)
    final_transform = reg.Execute(fixed_sd, moving_sd)

    resampled = sitk.Resample(
        sitk.Cast(moving_mask > 0, sitk.sitkUInt8),
        fixed_mask,
        final_transform,
        sitk.sitkNearestNeighbor,
        0,
        sitk.sitkUInt8,
    )
    return resampled, final_transform


# ============================================================================
# VTP SURFACE EXPORT (VTK-based)
# ============================================================================


def export_surface_vtp(
    binary_mask: sitk.Image,
    out_path: str,
    distance_image: Optional[sitk.Image] = None,
    scalar_name: str = 'dist_mm',
) -> None:
    """
    Export a marching-cubes surface mesh as VTP, optionally coloured by a
    scalar distance image.

    Requires VTK (``pip install vtk``).

    Parameters
    ----------
    binary_mask : sitk.Image
        Binary mask (>0 = foreground).
    out_path : str
        Output ``.vtp`` file path.
    distance_image : sitk.Image, optional
        Float image sampled at each vertex to produce per-vertex scalars.
        Typically a distance-to-other-mask map.
    scalar_name : str
        Name of the scalar array in the VTP file.
    """
    try:
        import vtk
        from vtk.util import numpy_support
    except ImportError:
        raise ImportError("VTK is required for VTP export: pip install vtk")

    mask_u8 = sitk.Cast(binary_mask > 0, sitk.sitkUInt8)
    arr = sitk.GetArrayFromImage(mask_u8)  # z, y, x
    depth, height, width = arr.shape

    # Build vtkImageData
    vtk_img = vtk.vtkImageData()
    vtk_img.SetDimensions(width, height, depth)
    vtk_img.SetSpacing(mask_u8.GetSpacing())
    vtk_img.SetOrigin(mask_u8.GetOrigin())

    flat = arr.ravel(order='C')
    vtk_arr = numpy_support.numpy_to_vtk(
        num_array=flat, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR,
    )
    vtk_arr.SetName('mask')
    vtk_img.GetPointData().SetScalars(vtk_arr)

    # Marching cubes
    try:
        mc = vtk.vtkFlyingEdges3D()
    except AttributeError:
        mc = vtk.vtkMarchingCubes()
    mc.SetInputData(vtk_img)
    mc.SetValue(0, 0.5)
    mc.Update()

    # Clean + normals
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(mc.GetOutput())
    clean.Update()

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(clean.GetOutput())
    normals.AutoOrientNormalsOn()
    normals.ConsistencyOn()
    normals.SplittingOff()
    normals.Update()

    poly = normals.GetOutput()

    # Attach distance scalars if provided
    if distance_image is not None:
        dist_arr = sitk.GetArrayViewFromImage(distance_image)  # z, y, x
        spacing = distance_image.GetSpacing()
        origin = distance_image.GetOrigin()
        size = distance_image.GetSize()  # x, y, z

        npts = poly.GetNumberOfPoints()
        d = np.zeros(npts, dtype=np.float32)
        for i in range(npts):
            x, y, z = poly.GetPoint(i)
            ix = int(round((x - origin[0]) / spacing[0]))
            iy = int(round((y - origin[1]) / spacing[1]))
            iz = int(round((z - origin[2]) / spacing[2]))
            if 0 <= ix < size[0] and 0 <= iy < size[1] and 0 <= iz < size[2]:
                d[i] = float(dist_arr[iz, iy, ix])

        vtk_d = numpy_support.numpy_to_vtk(d, deep=True, array_type=vtk.VTK_FLOAT)
        vtk_d.SetName(scalar_name)
        poly.GetPointData().AddArray(vtk_d)
        poly.GetPointData().SetActiveScalars(scalar_name)

    # Write VTP
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(out_path)
    writer.SetInputData(poly)
    writer.Write()
