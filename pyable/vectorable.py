"""
Vector Field and Time Series Image Classes for Pyable

Extends pyable's image processing capabilities to handle:
  - Vector Fields (displacement fields, velocity fields, etc.)
  - Time Series Images (4D temporal sequences)

Both classes maintain the same API as Imaginable (transformations, filters, 
resampling) while handling specialized data types.

Classes:
    - Vectorable: Handle 3D/2D vector fields (VectorFloat64, VectorFloat32)
    - TimeSeriesable: Handle 4D time-series images (multiple 3D frames)

Example:
    >>> # Vector field operations
    >>> vf = Vectorable('displacement.mha')
    >>> vf.scaleVector([2, 2, 2])
    >>> vf.applyGaussianSmoothing(sigma=1.0)
    >>> vf.write('smoothed_displacement.mha')
    >>> 
    >>> # Time series operations
    >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
    >>> print(f"Frames: {ts.getNumberOfFrames()}")
    >>> frame = ts.getFrame(0)
    >>> ts.applyFilterToAllFrames('gaussian', sigma=1.0)
    >>> ts.write('filtered_cardiac.nii.gz')
"""

import SimpleITK as sitk
import numpy as np
from pathlib import Path
from typing import Union, List, Tuple, Optional
import copy

try:
    from .imaginable import Imaginable
except ImportError:
    from imaginable import Imaginable


# ============================================================================
# VECTORABLE - Vector Field Class
# ============================================================================

class Vectorable(Imaginable):
    """
    Handle vector fields (displacement fields, velocity fields, vector sequences).
    
    Maintains all Imaginable capabilities (transformations, filters, I/O) while
    specializing for vector data (VectorFloat64, VectorFloat32).
    
    Key Methods:
        - Geometry operations: scaleVector(), rotateVector(), translateVector()
        - Vector-specific: getMagnitude(), getComponent(), setComponent()
        - Filtering: applyGaussianSmoothing(), normalize()
        - Statistics: getVectorStatistics(), getMeanVector()
    
    Example:
        >>> vf = Vectorable('displacement_field.mha')
        >>> # Get magnitude at each voxel
        >>> magnitude = vf.getMagnitude()
        >>> # Scale all vectors by 2x
        >>> vf.scaleVector([2, 2, 2])
        >>> # Apply smoothing
        >>> vf.applyGaussianSmoothing(sigma=1.5)
        >>> vf.write('processed_field.mha')
    """
    
    def __init__(self, filename: Optional[str] = None, image: Optional[sitk.Image] = None, 
                 verbose: bool = False):
        """
        Initialize Vectorable object.
        
        Parameters
        ----------
        filename : str, optional
            Path to vector field file (.mha, .nii.gz, .h5)
        image : sitk.Image, optional
            SimpleITK vector image
        verbose : bool, default=False
            Enable verbose output
        """
        super().__init__(filename, image, verbose)
        
        # Verify it's a vector image if provided
        current_image = self.getImage()
        if current_image is not None:
            if current_image.GetNumberOfComponentsPerPixel() == 1:
                raise ValueError("Vectorable expects vector images (components > 1), use Imaginable for scalar images")
    
    def getNumberOfComponents(self) -> int:
        """Get number of vector components (typically 2 or 3)."""
        return self.getImage().GetNumberOfComponentsPerPixel()
    
    def getMagnitude(self) -> 'Imaginable':
        """
        Get magnitude of each vector.
        
        Returns
        -------
        magnitude_image : Imaginable
            Scalar image with magnitude at each voxel
        
        Example
        -------
        >>> vf = Vectorable('displacement.mha')
        >>> mag = vf.getMagnitude()
        >>> print(f"Max displacement: {mag.getMaximumValue()}")
        """
        # Use SimpleITK's magnitude computation
        magnitude_filter = sitk.VectorMagnitudeImageFilter()
        magnitude_image = magnitude_filter.Execute(self.getImage())
        
        result = Imaginable(image=magnitude_image)
        return result
    
    def getComponent(self, component: int) -> 'Imaginable':
        """
        Extract a single vector component as scalar image.
        
        Parameters
        ----------
        component : int
            Component index (0, 1, 2, ...)
        
        Returns
        -------
        component_image : Imaginable
            Scalar image with selected component
        
        Example
        -------
        >>> vf = Vectorable('3d_displacement.mha')  # 3D vectors
        >>> x_component = vf.getComponent(0)
        >>> y_component = vf.getComponent(1)
        >>> z_component = vf.getComponent(2)
        """
        if component >= self.getNumberOfComponents():
            raise ValueError(f"Component {component} out of range (max: {self.getNumberOfComponents() - 1})")
        
        # Use VectorIndexSelectionCastImageFilter
        selector = sitk.VectorIndexSelectionCastImageFilter()
        selector.SetIndex(component)
        component_image = selector.Execute(self.getImage())
        
        result = Imaginable(image=component_image)
        return result
    
    def setComponent(self, component: int, scalar_image: Union[sitk.Image, 'Imaginable']) -> 'Vectorable':
        """
        Set a single component of the vector field.
        
        Parameters
        ----------
        component : int
            Component index to set
        scalar_image : sitk.Image or Imaginable
            Scalar image for this component
        
        Returns
        -------
        self : Vectorable
            Self for method chaining
        
        Example
        -------
        >>> vf = Vectorable('displacement.mha')
        >>> x_comp = Imaginable('x_component.nii.gz')
        >>> vf.setComponent(0, x_comp)
        """
        if isinstance(scalar_image, Imaginable):
            scalar_image = scalar_image.getImage()
        
        # Compose vectors from components
        components = [self.getComponent(i).getImage() 
                     for i in range(self.getNumberOfComponents())]
        components[component] = scalar_image
        
        vector_image = sitk.Compose(*components)
        return self.setImage(vector_image, f"component {component} updated")
    
    def scaleVector(self, factors: Union[List[float], Tuple[float, ...], float] = 1.0) -> 'Vectorable':
        """
        Scale all vectors by given factors.
        
        Parameters
        ----------
        factors : float or list of floats
            Scale factor(s). If float, applies to all components.
            If list, must match number of components.
        
        Returns
        -------
        self : Vectorable
            Self for method chaining
        
        Example
        -------
        >>> vf = Vectorable('displacement.mha')
        >>> vf.scaleVector([2.0, 2.0, 2.0])  # 2x scaling
        >>> vf.scaleVector(0.5)  # 0.5x scaling (all components)
        """
        if isinstance(factors, (int, float)):
            factors = [factors] * self.getNumberOfComponents()
        
        if len(factors) != self.getNumberOfComponents():
            raise ValueError(f"Expected {self.getNumberOfComponents()} factors, got {len(factors)}")
        
        # Scale using numpy for correctness
        arr = sitk.GetArrayFromImage(self.getImage())
        for i, f in enumerate(factors):
            arr[..., i] = arr[..., i] * f
        result = sitk.GetImageFromArray(arr, isVector=True)
        result.CopyInformation(self.getImage())
        
        return self.setImage(result, f"vectors scaled by {factors}")
    
    def normalize(self, target_magnitude: float = 1.0) -> 'Vectorable':
        """
        Normalize all vectors to unit length or target magnitude.
        
        Parameters
        ----------
        target_magnitude : float, default=1.0
            Target magnitude for all vectors
        
        Returns
        -------
        self : Vectorable
            Self for method chaining
        
        Example
        -------
        >>> vf = Vectorable('velocity_field.mha')
        >>> vf.normalize(target_magnitude=1.0)  # Unit vectors
        """
        magnitude_image = sitk.VectorMagnitudeImageFilter().Execute(self.getImage())
        
        # Use numpy for the whole operation to avoid geometry issues
        vector_array = sitk.GetArrayFromImage(self.getImage())
        magnitude_array = sitk.GetArrayFromImage(magnitude_image)
        
        # Avoid division by zero
        magnitude_array[magnitude_array == 0] = 1
        
        # Normalize
        with np.errstate(divide='ignore', invalid='ignore'):
            normalized = vector_array / magnitude_array[..., np.newaxis] * target_magnitude
            normalized = np.nan_to_num(normalized)
        
        result_image = sitk.GetImageFromArray(normalized, isVector=True)
        result_image.CopyInformation(self.getImage())
        
        return self.setImage(result_image, f"vectors normalized to {target_magnitude}")
    
    def applyGaussianSmoothing(self, sigma: Union[float, List[float]] = 1.0, use_spacing: bool = True) -> 'Vectorable':
        """
        Apply Gaussian smoothing to vector field.
        
        Parameters
        ----------
        sigma : float or list of floats, default=1.0
            Standard deviation for Gaussian kernel
        use_spacing : bool, default=True
            Use image spacing for kernel size
        
        Returns
        -------
        self : Vectorable
            Self for method chaining
        
        Example
        -------
        >>> vf = Vectorable('noisy_displacement.mha')
        >>> vf.applyGaussianSmoothing(sigma=[1.5, 1.5, 1.5])
        """
        if isinstance(sigma, (int, float)):
            sigma = [sigma] * self.getImageDimension()
        
        gaussian_filter = sitk.SmoothingRecursiveGaussianImageFilter()
        gaussian_filter.SetSigma(sigma)
        gaussian_filter.SetNormalizeAcrossScale(False)
        
        smoothed = gaussian_filter.Execute(self.getImage())
        return self.setImage(smoothed, f"gaussian smoothing (sigma={sigma})")
    
    def getVectorStatistics(self) -> dict:
        """
        Compute statistics on vector field magnitudes.
        
        Returns
        -------
        stats : dict
            Dictionary with 'mean', 'std', 'min', 'max' magnitude
        
        Example
        -------
        >>> vf = Vectorable('displacement.mha')
        >>> stats = vf.getVectorStatistics()
        >>> print(f"Mean displacement: {stats['mean']:.2f}")
        """
        magnitude = sitk.VectorMagnitudeImageFilter().Execute(self.getImage())
        magnitude_array = sitk.GetArrayFromImage(magnitude)
        
        return {
            'mean': float(np.mean(magnitude_array)),
            'std': float(np.std(magnitude_array)),
            'min': float(np.min(magnitude_array)),
            'max': float(np.max(magnitude_array)),
        }

    # Backward-compatible alias used by legacy callers.
    getStatistics = getVectorStatistics
    
    def getMeanVector(self) -> np.ndarray:
        """
        Get mean vector across entire field.
        
        Returns
        -------
        mean_vector : np.ndarray
            Mean of all vectors
        
        Example
        -------
        >>> vf = Vectorable('velocity_field.mha')
        >>> mean_vel = vf.getMeanVector()
        >>> print(f"Mean velocity: {mean_vel}")
        """
        vector_array = sitk.GetArrayFromImage(self.getImage())
        
        # Reshape to (n_voxels, n_components)
        original_shape = vector_array.shape
        if len(original_shape) == 4:  # 3D with components
            vector_array_flat = vector_array.reshape(-1, original_shape[-1])
        else:
            vector_array_flat = vector_array.reshape(-1, original_shape[-1])
        
        return np.mean(vector_array_flat, axis=0)
    
    def describe(self):
        """Print a concise summary of the vector field.

        Returns:
            dict: key vector field properties
        """
        info = super().describe() if hasattr(super(), 'describe') else {}
        try:
            info['num_components'] = self.getNumberOfComponents()
            stats = self.getVectorStatistics()
            if stats:
                for k, v in stats.items():
                    info[f'magnitude_{k}'] = v
        except Exception as e:
            info['vector_error'] = str(e)
        self._print_describe(info)
        return info

    def getDuplicate(self) -> 'Vectorable':
        """Create a copy of this vector field."""
        return Vectorable(image=sitk.Image(self.getImage()))
    
    def plotOverlay(self, overlay=None, alpha=0.5, title=None, component=None,
                   slice_idx=None, **kwargs):
        """
        Display vector field with optional overlay using interactive viewer.
        
        Parameters
        ----------
        overlay : sitk.Image or Imaginable, optional
            Overlay image (will be resampled to match this vector field)
        alpha : float, default=0.5
            Overlay opacity (0-1)
        title : str, optional
            Figure title
        component : int, optional
            Component to display (0=X, 1=Y, 2=Z). If None, shows selector GUI.
        slice_idx : int, optional
            Slice index for 3D images (middle slice if None)
        **kwargs
            Additional arguments passed to viewer
        
        Returns
        -------
        viewer : VectorPlotter
            Viewer instance
        
        Examples
        --------
        >>> vf = Vectorable('displacement_field.mha')
        >>> 
        >>> # Interactive component selector
        >>> vf.plotOverlay()
        >>> 
        >>> # Show specific component
        >>> vf.plotOverlay(component=0)  # X component
        >>> 
        >>> # With overlay
        >>> roi = Imaginable('roi_mask.nii.gz')
        >>> vf.plotOverlay(roi, component=1, alpha=0.6)  # Y component with overlay
        """
        try:
            from .plotable import plotOverlay
        except ImportError:
            from plotable import plotOverlay
        
        if title is None:
            title = "Vector Field Viewer"
        
        return plotOverlay(self, overlay=overlay, alpha=alpha, 
                          title=title, component=component, 
                          slice_idx=slice_idx, **kwargs)

    def viewInteractive(self, overlays=None, orientation=2, slice_idx=None, 
                       component=None, title=None, figsize=(14, 10), cmap='gray'):
        """
        Open an interactive GUI viewer with vector-specific controls.
        
        Features:
        - Orientation selection (axial, sagittal, coronal)
        - Slice navigation with auto-center
        - Vector component selection (X, Y, Z, Magnitude)
        - Multiple overlay layers with individual opacity control
        - Real-time component updates
        
        Parameters
        ----------
        overlays : Imaginable, array, or list, optional
            Single or multiple overlays to display with the vector field
        orientation : int, default=2
            Initial viewing orientation (0=axial, 1=sagittal, 2=coronal)
        slice_idx : int, optional
            Initial slice index. If None, uses center slice.
        component : int, optional
            Initial component to display (0=X, 1=Y, 2=Z)
        title : str, optional
            Window title. Auto-generated if None.
        figsize : tuple, default=(14, 10)
            Figure size in inches (width, height)
        cmap : str, default='gray'
            Colormap for display
            
        Returns
        -------
        InteractiveViewer
            Viewer instance
            
        Example
        -------
        >>> vf = Vectorable('displacement_field.mha')
        >>> img = Imaginable('reference_image.nii.gz')
        >>> vf.viewInteractive(overlays=img, component=0, orientation=2)
        """
        try:
            from .interactive_viewer import InteractiveViewer
        except ImportError:
            from interactive_viewer import InteractiveViewer
        
        if title is None:
            title = "Vector Field Viewer - Interactive"
        
        viewer = InteractiveViewer(self.getImage(), title=title, 
                                   figsize=figsize, cmap=cmap)
        viewer.current_orientation = orientation
        viewer.current_slice = slice_idx or viewer._get_center_slice(orientation)
        if component is not None:
            viewer.current_component = component
        
        # Add overlays
        if overlays is not None:
            if isinstance(overlays, (list, tuple)):
                viewer.add_overlays(overlays)
            else:
                viewer.add_overlay(overlays)
        
        viewer.show()
        return viewer


class TimeSeriesable(Imaginable):
    """
    Handle 4D time-series images (e.g., cardiac sequences, dynamic MRI).
    
    Treats the 4th dimension as time and provides:
      - Per-frame access and operations
      - Temporal filtering and interpolation
      - Frame averaging and statistics
      - Per-frame transformations
    
    All inherited Imaginable methods work on the full 4D image.
    Use frame-specific methods for individual timepoints.
    
    Key Methods:
        - Frame access: getFrame(), setFrame(), getNumberOfFrames()
        - Frame operations: applyFilterToAllFrames(), transformAllFrames()
        - Temporal: getTemporalMean(), getTemporalVariance()
        - Extraction: extractPhase(), getFrameRange()
    
    Example:
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> print(f"Frames: {ts.getNumberOfFrames()}")
        >>> 
        >>> # Get specific frame
        >>> diastole = ts.getFrame(0)
        >>> systole = ts.getFrame(5)
        >>> 
        >>> # Apply filter to all frames
        >>> ts.applyFilterToAllFrames('gaussian', sigma=1.0)
        >>> 
        >>> # Get temporal statistics
        >>> mean_vol = ts.getTemporalMean()
    """
    
    def __init__(self, filename: Optional[str] = None, image: Optional[sitk.Image] = None,
                 verbose: bool = False):
        """
        Initialize TimeSeriesable object.
        
        Parameters
        ----------
        filename : str, optional
            Path to 4D image file
        image : sitk.Image, optional
            SimpleITK 4D image
        verbose : bool, default=False
            Enable verbose output
        """
        super().__init__(filename, image, verbose)
        
        # Verify it's 4D
        current_image = self.getImage()
        if current_image is not None and current_image.GetDimension() != 4:
            raise ValueError(f"TimeSeriesable expects 4D images, got {current_image.GetDimension()}D")

    def _build_4d_image(self, array_4d: np.ndarray, start_frame: int = 0) -> sitk.Image:
        """Create a 4D scalar SimpleITK image while preserving time-series metadata."""
        image_4d = sitk.GetImageFromArray(array_4d, isVector=False)
        source = self.getImage()
        origin = list(source.GetOrigin())
        spacing = list(source.GetSpacing())
        if len(origin) > 3 and len(spacing) > 3:
            origin[3] = origin[3] + spacing[3] * start_frame
        image_4d.SetOrigin(tuple(origin))
        image_4d.SetSpacing(tuple(spacing))
        image_4d.SetDirection(source.GetDirection())
        return image_4d

    def _build_spatial_image(self, array_3d: np.ndarray) -> sitk.Image:
        """Create a 3D scalar image from a temporal reduction."""
        image_3d = sitk.GetImageFromArray(array_3d, isVector=False)
        source = self.getImage()
        direction_4d = np.asarray(source.GetDirection(), dtype=float).reshape((4, 4))
        image_3d.SetOrigin(tuple(source.GetOrigin()[:3]))
        image_3d.SetSpacing(tuple(source.GetSpacing()[:3]))
        image_3d.SetDirection(tuple(direction_4d[:3, :3].reshape(-1)))
        return image_3d
    
    def getNumberOfFrames(self) -> int:
        """
        Get number of time frames.
        
        Returns
        -------
        n_frames : int
            Number of frames in the time dimension
        """
        size = self.getImageSize()
        return size[3]  # 4th dimension
    
    def getFrame(self, frame_index: int) -> 'Imaginable':
        """
        Extract a single frame as 3D image.
        
        Parameters
        ----------
        frame_index : int
            Frame number (0-indexed)
        
        Returns
        -------
        frame : Imaginable
            3D image for the given frame
        
        Example
        -------
        >>> ts = TimeSeriesable('4d_cardiac.nii.gz')
        >>> frame0 = ts.getFrame(0)
        >>> frame5 = ts.getFrame(5)
        """
        if frame_index < 0 or frame_index >= self.getNumberOfFrames():
            raise IndexError(f"Frame {frame_index} out of range ({self.getNumberOfFrames()} frames)")
        
        # Use ExtractImageFilter to get single frame
        extractor = sitk.ExtractImageFilter()
        size = list(self.getImageSize())
        size[3] = 0  # Extract single slice in time
        
        index = [0, 0, 0, frame_index]
        extractor.SetSize(size)
        extractor.SetIndex(index)
        
        frame_image = extractor.Execute(self.getImage())
        # ExtractImageFilter already returns a proper 3D image with data and geometry
        
        result = Imaginable(image=frame_image)
        return result
    
    def setFrame(self, frame_index: int, frame_image: Union[sitk.Image, 'Imaginable']) -> 'TimeSeriesable':
        """
        Replace a single frame.
        
        Parameters
        ----------
        frame_index : int
            Frame number to replace
        frame_image : sitk.Image or Imaginable
            3D replacement image
        
        Returns
        -------
        self : TimeSeriesable
            Self for method chaining
        
        Example
        -------
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> processed_frame = Imaginable('frame0_processed.nii.gz')
        >>> ts.setFrame(0, processed_frame)
        """
        if isinstance(frame_image, Imaginable):
            frame_image = frame_image.getImage()

        if frame_index < 0 or frame_index >= self.getNumberOfFrames():
            raise IndexError(f"Frame {frame_index} out of range ({self.getNumberOfFrames()} frames)")
        
        # Get current 4D image as numpy
        current_4d = sitk.GetArrayFromImage(self.getImage())
        new_frame = sitk.GetArrayFromImage(frame_image)

        if new_frame.shape != current_4d[frame_index].shape:
            raise ValueError(
                f"Frame shape mismatch: expected {current_4d[frame_index].shape}, got {new_frame.shape}"
            )
        
        # Replace frame
        current_4d[frame_index] = new_frame
        
        # Convert back
        result_image = self._build_4d_image(current_4d)
        
        return self.setImage(result_image, f"frame {frame_index} replaced")
    
    def getFrameRange(self, start: int, end: int) -> 'TimeSeriesable':
        """
        Extract a contiguous range of frames.
        
        Parameters
        ----------
        start : int
            Starting frame (inclusive)
        end : int
            Ending frame (exclusive)
        
        Returns
        -------
        frame_range : TimeSeriesable
            New TimeSeriesable with subset of frames
        
        Example
        -------
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> systole_range = ts.getFrameRange(3, 7)  # Frames 3-6
        """
        if start < 0 or end > self.getNumberOfFrames() or start >= end:
            raise IndexError(
                f"Invalid frame range [{start}, {end}) for {self.getNumberOfFrames()} frames"
            )

        # Extract subset using numpy indexing
        array_4d = sitk.GetArrayFromImage(self.getImage())
        subset = array_4d[start:end].copy()
        
        result_image = self._build_4d_image(subset, start_frame=start)
        
        result = TimeSeriesable(image=result_image)
        return result
    
    def getTemporalMean(self) -> 'Imaginable':
        """
        Compute temporal mean (average across all frames).
        
        Returns
        -------
        mean_image : Imaginable
            3D image with mean across time
        
        Example
        -------
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> mean_vol = ts.getTemporalMean()
        >>> mean_vol.write('mean_volume.nii.gz')
        """
        array_4d = sitk.GetArrayFromImage(self.getImage())
        mean_array = np.mean(array_4d, axis=0)
        
        result_image = self._build_spatial_image(mean_array)
        
        result = Imaginable(image=result_image)
        return result
    
    def getTemporalVariance(self) -> 'Imaginable':
        """
        Compute temporal variance across frames.
        
        Returns
        -------
        variance_image : Imaginable
            3D image with variance across time
        """
        array_4d = sitk.GetArrayFromImage(self.getImage())
        var_array = np.var(array_4d, axis=0)
        
        result_image = self._build_spatial_image(var_array)
        
        result = Imaginable(image=result_image)
        return result
    
    def getTemporalStandardDeviation(self) -> 'Imaginable':
        """Compute temporal standard deviation across frames."""
        variance = self.getTemporalVariance()
        variance_img = variance.getImage()
        std_image = sitk.Sqrt(variance_img)
        return Imaginable(image=std_image)
    
    def applyFilterToAllFrames(self, filter_name: str, **kwargs) -> 'TimeSeriesable':
        """
        Apply a filter independently to each frame.
        
        Parameters
        ----------
        filter_name : str
            Name of filter: 'gaussian', 'median', 'bilateral'
        **kwargs
            Filter-specific parameters
        
        Returns
        -------
        self : TimeSeriesable
            Self for method chaining
        
        Example
        -------
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> ts.applyFilterToAllFrames('gaussian', sigma=1.0)
        >>> ts.applyFilterToAllFrames('median', radius=2)
        """
        n_frames = self.getNumberOfFrames()
        
        # Process each frame
        for i in range(n_frames):
            frame = self.getFrame(i)
            
            # Apply filter based on name
            if filter_name.lower() == 'gaussian':
                sigma = kwargs.get('sigma', 1.0)
                gaussian = sitk.SmoothingRecursiveGaussianImageFilter()
                gaussian.SetSigma(sigma)
                filtered = gaussian.Execute(frame.getImage())
            elif filter_name.lower() == 'median':
                radius = kwargs.get('radius', 1)
                median = sitk.MedianImageFilter()
                median.SetRadius(radius)
                filtered = median.Execute(frame.getImage())
            elif filter_name.lower() == 'bilateral':
                domain_sigma = kwargs.get('domain_sigma', 1.0)
                range_sigma = kwargs.get('range_sigma', 1.0)
                bilateral = sitk.BilateralImageFilter()
                bilateral.SetDomainSigma(domain_sigma)
                bilateral.SetRangeSigma(range_sigma)
                filtered = bilateral.Execute(frame.getImage())
            else:
                raise ValueError(f"Unknown filter: {filter_name}")
            
            frame.setImage(filtered)
            self.setFrame(i, frame)
        
        return self.setImage(self.getImage(), f"filter '{filter_name}' applied to all frames")
    
    def transformAllFrames(self, transform, interpolator: str = 'linear') -> 'TimeSeriesable':
        """
        Apply same transformation to all frames.
        
        Parameters
        ----------
        transform : str or sitk.Transform
            Transformation to apply
        interpolator : str, default='linear'
            Interpolation method
        
        Returns
        -------
        self : TimeSeriesable
            Self for method chaining
        
        Example
        -------
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> ts.transformAllFrames('transform.tfm')
        """
        from . import deformations
        
        n_frames = self.getNumberOfFrames()
        
        for i in range(n_frames):
            frame = self.getFrame(i)
            warped = deformations.apply_transform(
                frame.getImage(),
                transform,
                interpolator=interpolator
            )
            frame.setImage(warped)
            self.setFrame(i, frame)
        
        return self.setImage(self.getImage(), "transform applied to all frames")
    
    def extractPhase(self, phase_number: int) -> 'Imaginable':
        """
        Extract a specific cardiac phase (frame).
        
        Alias for getFrame() with semantic meaning for cardiac imaging.
        
        Parameters
        ----------
        phase_number : int
            Phase index (0-indexed)
        
        Returns
        -------
        phase_image : Imaginable
            3D image for the phase
        """
        return self.getFrame(phase_number)
    
    def describe(self):
        """Print a concise summary of the time series.

        Returns:
            dict: key time series properties
        """
        info = super().describe() if hasattr(super(), 'describe') else {}
        try:
            info['num_frames'] = self.getNumberOfFrames()
        except Exception as e:
            info['timeseries_error'] = str(e)
        self._print_describe(info)
        return info

    def getDuplicate(self) -> 'TimeSeriesable':
        """Create a copy of this time series."""
        return TimeSeriesable(image=sitk.Image(self.getImage()))
    
    def plotOverlay(self, overlay=None, alpha=0.5, title=None, frame=None, 
                   slice_idx=None, **kwargs):
        """
        Display time-series with optional overlay using interactive viewer.
        
        Parameters
        ----------
        overlay : sitk.Image or Imaginable, optional
            Overlay image (will be resampled to match current frame)
        alpha : float, default=0.5
            Overlay opacity (0-1)
        title : str, optional
            Figure title
        frame : int, optional
            Initial frame to display (0-based)
        slice_idx : int, optional
            Slice index for 3D frames (middle slice if None)
        **kwargs
            Additional arguments passed to viewer
        
        Returns
        -------
        viewer : TimeSeriesPlotter
            Viewer instance
        
        Examples
        --------
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> ts.plotOverlay()  # Interactive frame selector
        >>> ts.plotOverlay(frame=5)  # Start at frame 5
        >>> 
        >>> # With overlay
        >>> roi = Imaginable('roi_mask.nii.gz')
        >>> ts.plotOverlay(roi, alpha=0.6)
        """
        try:
            from .plotable import plotOverlay
        except ImportError:
            from plotable import plotOverlay
        
        if title is None:
            title = "Time Series Viewer"
        
        return plotOverlay(self, overlay=overlay, alpha=alpha, 
                          title=title, frame=frame, slice_idx=slice_idx, **kwargs)

    def viewInteractive(self, overlays=None, orientation=2, slice_idx=None, 
                       frame=None, title=None, figsize=(14, 10), cmap='gray'):
        """
        Open an interactive GUI viewer with time-series-specific controls.
        
        Features:
        - Orientation selection (axial, sagittal, coronal)
        - Slice navigation with auto-center
        - Time frame selection and navigation
        - Multiple overlay layers with individual opacity control
        - Real-time frame updates
        
        Parameters
        ----------
        overlays : Imaginable, array, or list, optional
            Single or multiple overlays to display with time series
        orientation : int, default=2
            Initial viewing orientation (0=axial, 1=sagittal, 2=coronal)
        slice_idx : int, optional
            Initial slice index. If None, uses center slice.
        frame : int, optional
            Initial frame to display. If None, shows first frame (0).
        title : str, optional
            Window title. Auto-generated if None.
        figsize : tuple, default=(14, 10)
            Figure size in inches (width, height)
        cmap : str, default='gray'
            Colormap for display
            
        Returns
        -------
        InteractiveViewer
            Viewer instance
            
        Example
        -------
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> roi = Imaginable('roi_mask.nii.gz')
        >>> ts.viewInteractive(overlays=roi, frame=5, orientation=0)
        """
        try:
            from .interactive_viewer import InteractiveViewer
        except ImportError:
            from interactive_viewer import InteractiveViewer
        
        if title is None:
            title = "Time Series Viewer - Interactive"
        
        viewer = InteractiveViewer(self.getImage(), title=title, 
                                   figsize=figsize, cmap=cmap)
        viewer.current_orientation = orientation
        viewer.current_slice = slice_idx or viewer._get_center_slice(orientation)
        if frame is not None:
            viewer.current_frame = frame
        
        # Add overlays
        if overlays is not None:
            if isinstance(overlays, (list, tuple)):
                viewer.add_overlays(overlays)
            else:
                viewer.add_overlay(overlays)
        
        viewer.show()
        return viewer


    pass
