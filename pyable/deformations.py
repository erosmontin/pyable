"""
Deformation Module for Pyable

Provides utilities for applying registration transforms and displacement fields
to Imaginable objects. Supports multiple transform types (affine, B-spline,
displacement fields) from medical image registration workflows.

Key Features:
    - Load transforms from various file formats (.tfm, .h5, .txt)
    - Convert transforms to displacement fields
    - Apply deformations with geometry alignment
    - Invert displacement fields
    - Batch deformation operations
    - ROI-aware deformation with label preservation

Example:
    >>> img = SITKImaginable('image.nii.gz')
    >>> # Apply affine transform
    >>> img.applyTransform('transform.tfm', interpolator='linear')
    >>> # Or apply displacement field
    >>> img.applyDisplacementField('deformation.mha', reference='fixed.nii.gz')
    >>> img.write('deformed.nii.gz')
"""

import SimpleITK as sitk
import numpy as np
from typing import Union, Tuple, Optional, List
from pathlib import Path


# ============================================================================
# DISPLACEMENT FIELD UTILITIES
# ============================================================================


def initialize_deformation_field(
    reference_image: Union[sitk.Image, str],
    vector_components: int = 3,
) -> Tuple[sitk.Image, tuple, tuple, tuple, tuple]:
    """
    Create an empty (zero) deformation field from a reference image.

    Parameters
    ----------
    reference_image : sitk.Image or str
        Reference image defining the geometry, or path to image file
    vector_components : int, default=3
        Number of vector components (2 for 2D, 3 for 3D)

    Returns
    -------
    displacement_field : sitk.Image
        Zero-initialized vector image
    size : tuple
        Image size
    origin : tuple
        Image origin
    spacing : tuple
        Image spacing
    direction : tuple
        Image direction matrix

    Example
    -------
    >>> df, size, origin, spacing, direction = initialize_deformation_field('image.nii.gz')
    >>> print(f"Displacement field created with size {size}")
    """
    if isinstance(reference_image, str):
        reference_image = sitk.ReadImage(reference_image)

    size = reference_image.GetSize()
    origin = reference_image.GetOrigin()
    spacing = reference_image.GetSpacing()
    direction = reference_image.GetDirection()

    # Create vector image with specified components
    pixel_type = sitk.sitkVectorFloat64 if vector_components == 3 else sitk.sitkVectorFloat32
    displacement_field = sitk.Image(size, pixel_type, vector_components)

    # Copy geometry
    displacement_field.SetOrigin(origin)
    displacement_field.SetSpacing(spacing)
    displacement_field.SetDirection(direction)

    # Store metadata
    displacement_field.SetMetaData("TransformType", "DisplacementField")
    displacement_field.SetMetaData("TransformModel", "VectorFieldTransform")

    return displacement_field, size, origin, spacing, direction


def transform_to_displacement_field(
    transform: Union[sitk.Transform, str],
    output_size: tuple,
    output_origin: tuple,
    output_spacing: tuple,
    output_direction: tuple,
) -> sitk.Image:
    """
    Convert any SimpleITK transform to a displacement field.

    Useful for converting rigid, affine, B-spline, or other transforms
    into a dense displacement field for warping operations.

    Parameters
    ----------
    transform : sitk.Transform or str
        Transform object or path to transform file (.tfm, .h5)
    output_size : tuple
        Output displacement field size
    output_origin : tuple
        Output origin
    output_spacing : tuple
        Output spacing
    output_direction : tuple
        Output direction matrix

    Returns
    -------
    displacement_field : sitk.Image
        Dense displacement field matching output geometry

    Example
    -------
    >>> tfm = sitk.ReadTransform('affine.tfm')
    >>> df = transform_to_displacement_field(
    ...     tfm, (256, 256, 128), (0, 0, 0), (1, 1, 1), (1,0,0, 0,1,0, 0,0,1)
    ... )
    """
    if isinstance(transform, str):
        transform = sitk.ReadTransform(transform)

    displacement_field = sitk.TransformToDisplacementField(
        transform,
        size=output_size,
        outputOrigin=output_origin,
        outputSpacing=output_spacing,
        outputDirection=output_direction,
    )

    return displacement_field


def apply_deformation_field(
    image: Union[sitk.Image, str],
    displacement_field: Union[sitk.Image, str],
    target_image: Optional[Union[sitk.Image, str]] = None,
    interpolator: str = "linear",
    default_pixel_value: float = 0,
) -> sitk.Image:
    """
    Apply a displacement field to warp an image.

    Parameters
    ----------
    image : sitk.Image or str
        Image to warp
    displacement_field : sitk.Image or str
        Displacement field (vector image)
    target_image : sitk.Image or str, optional
        Target geometry reference. If None, uses input image geometry
    interpolator : str, default='linear'
        Interpolation method: 'linear', 'nearest', 'gaussian', 'bspline'
    default_pixel_value : float, default=0
        Value for pixels outside image domain

    Returns
    -------
    warped_image : sitk.Image
        Image warped by displacement field

    Example
    -------
    >>> img = sitk.ReadImage('moving.nii.gz')
    >>> df = sitk.ReadImage('displacement.mha')
    >>> warped = apply_deformation_field(img, df, interpolator='linear')
    """
    if isinstance(image, str):
        image = sitk.ReadImage(image)
    if isinstance(displacement_field, str):
        displacement_field = sitk.ReadImage(
            displacement_field, sitk.sitkVectorFloat64
        )

    # Set reference image
    reference_image = image
    if target_image is not None:
        if isinstance(target_image, str):
            target_image = sitk.ReadImage(target_image)
        reference_image = target_image

    # Configure warper
    warper = sitk.WarpImageFilter()
    warper.SetOutputParameteresFromImage(reference_image)

    # Set interpolator
    interpolators = {
        "linear": sitk.sitkLinear,
        "nearest": sitk.sitkNearestNeighbor,
        "gaussian": sitk.sitkGaussian,
        "bspline": sitk.sitkBSpline,
    }
    warper.SetInterpolator(interpolators.get(interpolator, sitk.sitkLinear))
    warper.SetEdgePaddingValue(default_pixel_value)

    warped_image = warper.Execute(image, displacement_field)
    return warped_image


def apply_transform(
    image: Union[sitk.Image, str],
    transform: Union[sitk.Transform, str],
    target_image: Optional[Union[sitk.Image, str]] = None,
    interpolator: str = "linear",
    default_pixel_value: float = 0,
) -> sitk.Image:
    """
    Apply a transform (affine, rigid, B-spline, etc.) to resample an image.

    Parameters
    ----------
    image : sitk.Image or str
        Image to transform
    transform : sitk.Transform or str
        Transform object or path to transform file (.tfm, .h5)
    target_image : sitk.Image or str, optional
        Target geometry reference. If None, uses input image geometry
    interpolator : str, default='linear'
        Interpolation method: 'linear', 'nearest', 'gaussian', 'bspline'
    default_pixel_value : float, default=0
        Value for pixels outside image domain

    Returns
    -------
    resampled_image : sitk.Image
        Image resampled via the transform

    Example
    -------
    >>> img = sitk.ReadImage('moving.nii.gz')
    >>> warped = apply_transform(img, 'transform.tfm', interpolator='linear')
    """
    if isinstance(image, str):
        image = sitk.ReadImage(image)
    if isinstance(transform, str):
        transform = sitk.ReadTransform(transform)

    # Set reference image
    reference_image = image
    if target_image is not None:
        if isinstance(target_image, str):
            target_image = sitk.ReadImage(target_image)
        reference_image = target_image

    # Configure resampler
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(reference_image)
    resampler.SetInterpolator(
        {
            "linear": sitk.sitkLinear,
            "nearest": sitk.sitkNearestNeighbor,
            "gaussian": sitk.sitkGaussian,
            "bspline": sitk.sitkBSpline,
        }.get(interpolator, sitk.sitkLinear)
    )
    resampler.SetDefaultPixelValue(default_pixel_value)
    resampler.SetTransform(transform)

    resampled_image = resampler.Execute(image)
    return resampled_image


# ============================================================================
# DISPLACEMENT FIELD INVERSION & GEOMETRY UTILITIES
# ============================================================================


def invert_displacement_field(
    displacement_field: Union[sitk.Image, str],
    target_image: Optional[Union[sitk.Image, str]] = None,
    max_iterations: int = 100,
    mean_error_tolerance: float = 1e-3,
) -> sitk.Image:
    """
    Invert a displacement field for reverse warping.

    Useful for forward-backward consistency checks and inverse warping.

    Parameters
    ----------
    displacement_field : sitk.Image or str
        Displacement field to invert
    target_image : sitk.Image or str, optional
        Target geometry. If None, uses displacement field geometry
    max_iterations : int, default=100
        Maximum iterations for inversion algorithm
    mean_error_tolerance : float, default=1e-3
        Tolerance for convergence

    Returns
    -------
    inverted_field : sitk.Image
        Inverted displacement field

    Example
    -------
    >>> df = sitk.ReadImage('forward_deform.mha')
    >>> df_inv = invert_displacement_field(df, max_iterations=100)
    """
    if isinstance(displacement_field, str):
        displacement_field = sitk.ReadImage(
            displacement_field, sitk.sitkVectorFloat64
        )

    # Use target geometry if provided
    if target_image is not None:
        if isinstance(target_image, str):
            target_image = sitk.ReadImage(target_image)
        size = target_image.GetSize()
        origin = target_image.GetOrigin()
        spacing = target_image.GetSpacing()
    else:
        size = displacement_field.GetSize()
        origin = displacement_field.GetOrigin()
        spacing = displacement_field.GetSpacing()

    # Invert using SimpleITK's built-in function
    inverted_field = sitk.InverseDisplacementField(
        displacement_field,
        size=size,
        outputOrigin=origin,
        outputSpacing=spacing,
        maxNumberOfIterations=max_iterations,
        meanErrorToleranceForInversion=mean_error_tolerance,
    )

    return inverted_field


def apply_inverted_deformation_field(
    image: Union[sitk.Image, str],
    displacement_field: Union[sitk.Image, str],
    target_image: Optional[Union[sitk.Image, str]] = None,
    interpolator: str = "linear",
    default_pixel_value: float = 0,
) -> sitk.Image:
    """
    Apply an inverted displacement field for reverse warping.

    Parameters
    ----------
    image : sitk.Image or str
        Image to warp
    displacement_field : sitk.Image or str
        Displacement field to invert and apply
    target_image : sitk.Image or str, optional
        Target geometry
    interpolator : str, default='linear'
        Interpolation method
    default_pixel_value : float, default=0
        Fill value

    Returns
    -------
    warped_image : sitk.Image
        Image warped by inverted displacement field
    """
    if isinstance(image, str):
        image = sitk.ReadImage(image)

    # Invert the field
    inverted_field = invert_displacement_field(displacement_field, target_image)

    # Apply inverted field
    return apply_deformation_field(
        image, inverted_field, target_image, interpolator, default_pixel_value
    )


def align_geometry(
    moving_image: Union[sitk.Image, str],
    reference_image: Union[sitk.Image, str],
) -> sitk.Image:
    """
    Align geometry (origin, spacing, direction) of moving image to reference.

    Useful for fixing displacement fields or images with incorrect metadata.

    Parameters
    ----------
    moving_image : sitk.Image or str
        Image with incorrect geometry
    reference_image : sitk.Image or str
        Image with correct geometry

    Returns
    -------
    aligned_image : sitk.Image
        Image with aligned geometry

    Example
    -------
    >>> df = sitk.ReadImage('deform.mha')  # May have lost geometry
    >>> fixed = sitk.ReadImage('fixed.nii.gz')
    >>> df_aligned = align_geometry(df, fixed)
    """
    if isinstance(moving_image, str):
        moving_image = sitk.ReadImage(moving_image, sitk.sitkVectorFloat64)
    if isinstance(reference_image, str):
        reference_image = sitk.ReadImage(reference_image)

    # Copy geometry
    moving_image.SetOrigin(reference_image.GetOrigin())
    moving_image.SetSpacing(reference_image.GetSpacing())
    moving_image.SetDirection(reference_image.GetDirection())

    return moving_image


def get_bspline_grid_info(transform: sitk.BSplineTransform) -> dict:
    """
    Extract grid information from a B-spline transform.

    Parameters
    ----------
    transform : sitk.BSplineTransform
        B-spline transform

    Returns
    -------
    grid_info : dict
        Dictionary with keys: 'origin', 'spacing', 'direction', 'mesh_size'
    """
    # GetTransformDomainPhysicalDimensions returns the total physical extent,
    # not the per-voxel spacing. Compute per-voxel spacing from extent / mesh_size.
    mesh_size = transform.GetTransformDomainMeshSize()
    phys_dims = transform.GetTransformDomainPhysicalDimensions()
    spacing = tuple(d / m for d, m in zip(phys_dims, mesh_size))
    return {
        "origin": transform.GetTransformDomainOrigin(),
        "physical_dimensions": phys_dims,
        "spacing": spacing,
        "direction": transform.GetTransformDomainDirection(),
        "mesh_size": mesh_size,
    }


def refine_bspline_grid(
    original_transform: sitk.BSplineTransform,
    new_mesh_size: Tuple[int, ...],
    reference_image: Optional[sitk.Image] = None,
) -> sitk.BSplineTransform:
    """
    Refine B-spline transform by resampling coefficients to new mesh size.

    Useful for mesh refinement in multi-resolution registration.

    Parameters
    ----------
    original_transform : sitk.BSplineTransform
        Original B-spline transform
    new_mesh_size : tuple
        New mesh size (e.g., (5, 5, 5))
    reference_image : sitk.Image, optional
        Reference image defining the domain. If None, one is constructed
        from the transform domain parameters.

    Returns
    -------
    refined_transform : sitk.BSplineTransform
        B-spline transform with refined mesh
    """
    # Get original coefficients
    original_coefficients = original_transform.GetCoefficientImages()

    # Build a reference image from the transform domain if not provided
    if reference_image is None:
        ndim = original_transform.GetDimension()
        phys_dims = original_transform.GetTransformDomainPhysicalDimensions()
        origin = original_transform.GetTransformDomainOrigin()
        direction = original_transform.GetTransformDomainDirection()
        mesh_size = original_transform.GetTransformDomainMeshSize()
        # Estimate a reasonable image size from the physical dimensions
        spacing = tuple(d / max(m, 1) for d, m in zip(phys_dims, mesh_size))
        size = [int(round(d / s)) for d, s in zip(phys_dims, spacing)]
        reference_image = sitk.Image(size, sitk.sitkFloat32)
        reference_image.SetOrigin(origin)
        reference_image.SetSpacing(spacing)
        reference_image.SetDirection(direction)

    # Create new B-spline transform with desired mesh size
    new_transform = sitk.BSplineTransformInitializer(
        reference_image,
        new_mesh_size,
    )

    # Extract grid information
    grid_info = get_bspline_grid_info(original_transform)

    # Resample coefficients to new grid
    resampler = sitk.ResampleImageFilter()
    new_coeff_ref = new_transform.GetCoefficientImages()[0]
    resampler.SetReferenceImage(new_coeff_ref)

    new_coefficients = [resampler.Execute(coef) for coef in original_coefficients]

    # Set new coefficients
    new_transform.SetCoefficientImages(new_coefficients)

    return new_transform


# ============================================================================
# COMPOSITE TRANSFORM & BATCH OPERATIONS
# ============================================================================


def create_composite_transform(
    transforms: List[Union[sitk.Transform, str]],
    inverse_flags: Optional[List[bool]] = None,
) -> sitk.CompositeTransform:
    """
    Create a composite transform from multiple transforms.

    Useful for applying multiple registration steps sequentially.

    Parameters
    ----------
    transforms : list
        List of transform objects or file paths
    inverse_flags : list, optional
        List of bools indicating whether to invert each transform

    Returns
    -------
    composite : sitk.CompositeTransform
        Composite transform

    Example
    -------
    >>> tfm1 = sitk.ReadTransform('rigid.tfm')
    >>> tfm2 = sitk.ReadTransform('deform.tfm')
    >>> composite = create_composite_transform([tfm1, tfm2])
    >>> warped = apply_transform(img, composite)
    """
    composite = sitk.CompositeTransform(3)

    for i, tfm in enumerate(transforms):
        if isinstance(tfm, str):
            tfm = sitk.ReadTransform(tfm)

        inverse = inverse_flags[i] if inverse_flags else False
        composite.AddTransform(tfm)
        if inverse:
            composite.FlattenTransformStack()
            # Note: SimpleITK composite transforms evaluate in reverse order

    return composite


def apply_multi_step_transform(
    image: Union[sitk.Image, str],
    transforms: List[Union[sitk.Transform, str]],
    target_image: Optional[Union[sitk.Image, str]] = None,
    interpolator: str = "linear",
) -> sitk.Image:
    """
    Apply multiple transforms sequentially to an image.

    Parameters
    ----------
    image : sitk.Image or str
        Image to transform
    transforms : list
        List of transforms (objects or file paths)
    target_image : sitk.Image or str, optional
        Target geometry
    interpolator : str, default='linear'
        Interpolation method

    Returns
    -------
    result : sitk.Image
        Image after all transforms applied

    Example
    -------
    >>> img = sitk.ReadImage('moving.nii.gz')
    >>> result = apply_multi_step_transform(
    ...     img,
    ...     ['rigid.tfm', 'deform.tfm'],
    ...     target_image='fixed.nii.gz'
    ... )
    """
    composite = create_composite_transform(transforms)
    return apply_transform(
        image, composite, target_image=target_image, interpolator=interpolator
    )


# ============================================================================
# LABEL-AWARE DEFORMATION (for ROI/Segmentation)
# ============================================================================


def apply_deformation_field_to_labels(
    label_image: Union[sitk.Image, str],
    displacement_field: Union[sitk.Image, str],
    target_image: Optional[Union[sitk.Image, str]] = None,
    interpolator: str = "nearest",
    default_label: int = 0,
) -> sitk.Image:
    """
    Apply displacement field to a label/segmentation map.

    Preserves label values using nearest-neighbor interpolation.

    Parameters
    ----------
    label_image : sitk.Image or str
        Label/segmentation image
    displacement_field : sitk.Image or str
        Displacement field
    target_image : sitk.Image or str, optional
        Target geometry
    interpolator : str, default='nearest'
        Use 'nearest' for labels
    default_label : int, default=0
        Background label value

    Returns
    -------
    warped_labels : sitk.Image
        Warped label image with preserved label values

    Example
    -------
    >>> labels = sitk.ReadImage('segmentation.nii.gz')
    >>> df = sitk.ReadImage('deformation.mha')
    >>> warped_labels = apply_deformation_field_to_labels(labels, df)
    """
    return apply_deformation_field(
        label_image,
        displacement_field,
        target_image=target_image,
        interpolator="nearest",
        default_pixel_value=default_label,
    )


def apply_transform_to_labels(
    label_image: Union[sitk.Image, str],
    transform: Union[sitk.Transform, str],
    target_image: Optional[Union[sitk.Image, str]] = None,
) -> sitk.Image:
    """
    Apply transform to a label/segmentation map.

    Preserves label values using nearest-neighbor interpolation.

    Parameters
    ----------
    label_image : sitk.Image or str
        Label/segmentation image
    transform : sitk.Transform or str
        Transform
    target_image : sitk.Image or str, optional
        Target geometry

    Returns
    -------
    warped_labels : sitk.Image
        Warped label image
    """
    return apply_transform(
        label_image, transform, target_image=target_image, interpolator="nearest"
    )
