"""
Interactive visualization and plotting module for pyable image classes.

Provides unified interface for plotting scalar images, vector fields, and time-series
with support for overlays, component selection, frame navigation, and interactive features.

Classes:
    - PlotViewer: Base viewer with overlay support and image saving
    - ScalarPlotter: For scalar (Imaginable) images
    - VectorPlotter: For vector fields (Vectorable) with component selection
    - TimeSeriesPlotter: For 4D time-series (TimeSeriesable) with frame navigation

Functions:
    - plotOverlay: Convenience function for plotting with overlays
    
Example:
    >>> from pyable_eros_montin import Imaginable, Vectorable
    >>> img = Imaginable('image.nii.gz')
    >>> img.plotOverlay('overlay.nii.gz')
    >>> 
    >>> vf = Vectorable('displacement.mha')
    >>> vf.plotOverlay(component=0)  # Show X-component with GUI
"""

import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from matplotlib.widgets import Slider, RadioButtons
from typing import Union, Optional, List, Tuple, Any
from pathlib import Path
import warnings


# ============================================================================
# PLOTVIEWER - Base Viewer Class
# ============================================================================

class PlotViewer:
    """
    Base viewer for image visualization with overlay support.
    
    Supports:
    - Single or multiple images
    - Overlay blending with alpha control
    - Slice selection for 3D images
    - Image statistics display
    - Save functionality
    
    Attributes
    ----------
    image : sitk.Image
        Primary image to display
    overlay : sitk.Image, optional
        Overlay image (automatically resampled to match primary)
    alpha : float
        Opacity of overlay (0-1)
    current_slice : int
        Current slice for 3D display
    cmap : str
        Colormap for primary image
    cmap_overlay : str
        Colormap for overlay
    """
    
    def __init__(self, image: sitk.Image, overlay: Optional[sitk.Image] = None,
                 title: str = "Image Viewer", cmap: str = 'gray', 
                 cmap_overlay: str = 'hot', figsize: Tuple[int, int] = (10, 8)):
        """
        Initialize PlotViewer.
        
        Parameters
        ----------
        image : sitk.Image
            Primary image (must be 2D or 3D scalar)
        overlay : sitk.Image, optional
            Overlay image (will be resampled to match primary)
        title : str, default="Image Viewer"
            Window title
        cmap : str, default='gray'
            Colormap for primary image
        cmap_overlay : str, default='hot'
            Colormap for overlay
        figsize : tuple, default=(10, 8)
            Figure size in inches (width, height)
        """
        self.image = image
        self.overlay = overlay
        self.title = title
        self.cmap = cmap
        self.cmap_overlay = cmap_overlay
        self.figsize = figsize
        self.alpha = 0.5
        self.current_slice = None
        self.fig = None
        self.ax = None
        self.slider = None
        
        # Validate and prepare images
        self._validate_images()
        self._resample_overlay()
    
    def _validate_images(self):
        """Validate that images are compatible."""
        if self.image.GetNumberOfComponentsPerPixel() != 1:
            raise ValueError("Primary image must be scalar (1 component)")
        
        if self.overlay is not None:
            if self.overlay.GetNumberOfComponentsPerPixel() != 1:
                raise ValueError("Overlay image must be scalar (1 component)")
            
            # Check dimension compatibility
            if self.image.GetDimension() != self.overlay.GetDimension():
                if not (self.image.GetDimension() == 3 and self.overlay.GetDimension() == 3):
                    raise ValueError(f"Image dimension mismatch: {self.image.GetDimension()} vs {self.overlay.GetDimension()}")
    
    def _resample_overlay(self):
        """Resample overlay to match primary image geometry."""
        if self.overlay is None:
            return
        
        # Check if already same size/spacing/origin
        if (self.overlay.GetSize() != self.image.GetSize() or
            self.overlay.GetSpacing() != self.image.GetSpacing() or
            self.overlay.GetOrigin() != self.image.GetOrigin()):
            
            resampler = sitk.ResampleImageFilter()
            resampler.SetReferenceImage(self.image)
            resampler.SetInterpolator(sitk.sitkLinear)
            self.overlay = resampler.Execute(self.overlay)
    
    def _get_2d_slice(self, image: sitk.Image, slice_idx: Optional[int] = None) -> np.ndarray:
        """
        Extract 2D slice from 2D or 3D image.
        
        Parameters
        ----------
        image : sitk.Image
            Source image (2D or 3D)
        slice_idx : int, optional
            Slice index for 3D images (uses middle slice if None)
        
        Returns
        -------
        slice_array : np.ndarray
            2D numpy array
        """
        if image.GetDimension() == 2:
            return sitk.GetArrayFromImage(image)
        
        elif image.GetDimension() == 3:
            if slice_idx is None:
                slice_idx = image.GetSize()[2] // 2
            
            self.current_slice = slice_idx
            
            # Extract slice in correct axis (SimpleITK uses ZYX ordering)
            extractor = sitk.ExtractImageFilter()
            size = list(image.GetSize())
            size[2] = 0  # Reduce Z dimension
            
            index = [0, 0, slice_idx]
            extractor.SetSize(size)
            extractor.SetIndex(index)
            
            slice_2d = extractor.Execute(image)
            return sitk.GetArrayFromImage(slice_2d)
        
        else:
            raise ValueError(f"Unsupported dimension: {image.GetDimension()}")
    
    def _normalize_array(self, arr: np.ndarray) -> np.ndarray:
        """Normalize array to [0, 1] range."""
        arr_min = np.min(arr)
        arr_max = np.max(arr)
        
        if arr_max == arr_min:
            return np.zeros_like(arr)
        
        return (arr - arr_min) / (arr_max - arr_min)
    
    def _on_slider_changed(self, val):
        """Callback for slice slider."""
        self.current_slice = int(val)
        self._update_display()
    
    def _update_display(self):
        """Update the displayed image and overlay."""
        self.ax.clear()
        
        # Get image slices
        img_array = self._get_2d_slice(self.image, self.current_slice)
        img_array = self._normalize_array(img_array)
        
        # Display primary image
        self.ax.imshow(img_array, cmap=self.cmap, origin='lower')
        
        # Overlay if available
        if self.overlay is not None:
            overlay_array = self._get_2d_slice(self.overlay, self.current_slice)
            overlay_array = self._normalize_array(overlay_array)
            
            self.ax.imshow(overlay_array, cmap=self.cmap_overlay, 
                          alpha=self.alpha, origin='lower')
        
        self.ax.set_title(f"{self.title}")
        if self.current_slice is not None:
            self.ax.set_title(f"{self.title} (Slice: {self.current_slice})")
        
        self.fig.canvas.draw_idle()
    
    def show(self, slice_idx: Optional[int] = None, alpha: float = 0.5):
        """
        Display the image viewer.
        
        Parameters
        ----------
        slice_idx : int, optional
            Initial slice index for 3D images
        alpha : float, default=0.5
            Initial overlay opacity (0-1)
        
        Examples
        --------
        >>> viewer = PlotViewer(image, overlay=overlay_image)
        >>> viewer.show(slice_idx=50, alpha=0.7)
        """
        self.alpha = alpha
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        
        # Set initial slice
        if self.image.GetDimension() == 3:
            if slice_idx is None:
                slice_idx = self.image.GetSize()[2] // 2
            self.current_slice = slice_idx
        
        self._update_display()
        
        # Add slider for 3D images
        if self.image.GetDimension() == 3:
            ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03])
            self.slider = Slider(ax_slider, 'Slice', 0, 
                                self.image.GetSize()[2] - 1,
                                valinit=self.current_slice, valstep=1)
            self.slider.on_changed(self._on_slider_changed)
        
        plt.tight_layout(rect=[0, 0.08, 1, 1])
        plt.show()
    
    def saveFigure(self, output_path: str, dpi: int = 150):
        """
        Save current figure to file.
        
        Parameters
        ----------
        output_path : str
            Path for saved image
        dpi : int, default=150
            Resolution in DPI
        """
        if self.fig is None:
            raise RuntimeError("No figure to save. Call show() first.")
        
        self.fig.savefig(output_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved figure to {output_path}")


# ============================================================================
# SCALARPLOTTER - Scalar Image Plotter
# ============================================================================

class ScalarPlotter(PlotViewer):
    """
    Specialized viewer for scalar (Imaginable) images.
    
    Features:
    - Single or multiple scalar images
    - Overlay with adjustable transparency
    - Automatic histogram equalization option
    - Statistics display
    
    Example
    -------
    >>> from pyable_eros_montin import Imaginable
    >>> img = Imaginable('image.nii.gz')
    >>> overlay = Imaginable('segmentation.nii.gz')
    >>> plotter = ScalarPlotter(img.getImage(), overlay.getImage())
    >>> plotter.show(alpha=0.6)
    """
    
    def __init__(self, image: Union[sitk.Image, 'Imaginable'],
                 overlay: Optional[Union[sitk.Image, 'Imaginable']] = None,
                 title: str = "Scalar Image Viewer", **kwargs):
        """
        Initialize ScalarPlotter.
        
        Parameters
        ----------
        image : sitk.Image or Imaginable
            Primary image
        overlay : sitk.Image or Imaginable, optional
            Overlay image
        title : str, default="Scalar Image Viewer"
            Window title
        **kwargs
            Additional arguments for PlotViewer
        """
        # Convert Imaginable to sitk.Image
        if hasattr(image, 'getImage'):
            image = image.getImage()
        if overlay is not None and hasattr(overlay, 'getImage'):
            overlay = overlay.getImage()
        
        super().__init__(image, overlay, title=title, **kwargs)


# ============================================================================
# VECTORPLOTTER - Vector Field Plotter
# ============================================================================

class VectorPlotter(PlotViewer):
    """
    Interactive viewer for vector fields (Vectorable).
    
    Features:
    - Component selection via GUI (X, Y, Z)
    - Vector magnitude display
    - Component-specific overlay support
    - Arrow overlay for vector field visualization
    
    Example
    -------
    >>> from pyable_eros_montin import Vectorable
    >>> vf = Vectorable('displacement_field.mha')
    >>> vf.plotOverlay()  # Interactive component selector
    """
    
    def __init__(self, vector_image: Union[sitk.Image, 'Vectorable'],
                 overlay: Optional[Union[sitk.Image, 'Imaginable']] = None,
                 title: str = "Vector Field Viewer", 
                 show_vectors: bool = False, **kwargs):
        """
        Initialize VectorPlotter.
        
        Parameters
        ----------
        vector_image : sitk.Image or Vectorable
            Vector field image
        overlay : sitk.Image or Imaginable, optional
            Overlay image
        title : str, default="Vector Field Viewer"
            Window title
        show_vectors : bool, default=False
            Display vector arrows on image
        **kwargs
            Additional arguments for PlotViewer
        """
        # Convert Vectorable to sitk.Image
        if hasattr(vector_image, 'getImage'):
            self.vectorable = vector_image
            vector_image = vector_image.getImage()
        else:
            self.vectorable = None
        
        if overlay is not None and hasattr(overlay, 'getImage'):
            overlay = overlay.getImage()
        
        # Validate vector image
        if vector_image.GetNumberOfComponentsPerPixel() < 2:
            raise ValueError("VectorPlotter requires vector image (components >= 2)")
        
        self.n_components = vector_image.GetNumberOfComponentsPerPixel()
        self.show_vectors = show_vectors
        self.current_component = 0
        self.vector_image = vector_image
        
        # Extract first component as scalar for initialization
        selector = sitk.VectorIndexSelectionCastImageFilter()
        selector.SetIndex(0)
        scalar_image = selector.Execute(vector_image)
        
        super().__init__(scalar_image, overlay, title=title, **kwargs)
        self.image = scalar_image
    
    def _update_component(self, component_idx: int):
        """Update displayed component."""
        self.current_component = component_idx
        
        selector = sitk.VectorIndexSelectionCastImageFilter()
        selector.SetIndex(component_idx)
        self.image = selector.Execute(self.vector_image)
        
        self._update_display()
    
    def _update_magnitude(self):
        """Display vector magnitude."""
        mag_filter = sitk.VectorMagnitudeImageFilter()
        self.image = mag_filter.Execute(self.vector_image)
        self._update_display()
    
    def show(self, component: int = 0, alpha: float = 0.5, 
             slice_idx: Optional[int] = None):
        """
        Display vector field viewer with component selector.
        
        Parameters
        ----------
        component : int, default=0
            Initial component to display (0=X, 1=Y, 2=Z)
        alpha : float, default=0.5
            Overlay opacity
        slice_idx : int, optional
            Initial slice for 3D images
        
        Example
        -------
        >>> vf = Vectorable('displacement.mha')
        >>> plotter = VectorPlotter(vf)
        >>> plotter.show(component=0)  # Start with X component
        """
        self.alpha = alpha
        self.fig, (self.ax, ax_radio) = plt.subplots(1, 2, 
                                                     figsize=(self.figsize[0] + 2, self.figsize[1]),
                                                     gridspec_kw={'width_ratios': [4, 1]})
        
        # Update to initial component
        self.current_component = component
        self._update_component(component)
        
        if self.image.GetDimension() == 3:
            if slice_idx is None:
                slice_idx = self.image.GetSize()[2] // 2
            self.current_slice = slice_idx
        
        self._update_display()
        
        # Component selector radio buttons
        component_names = ['X', 'Y', 'Z', 'Magnitude']
        component_names = component_names[:self.n_components] + ['Magnitude']
        
        def on_component_select(label):
            if label == 'Magnitude':
                self._update_magnitude()
            else:
                idx = ['X', 'Y', 'Z'].index(label)
                self._update_component(idx)
        
        radio = RadioButtons(ax_radio, component_names, active=component)
        radio.on_clicked(on_component_select)
        
        # Add slice slider for 3D
        if self.image.GetDimension() == 3:
            ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03])
            self.slider = Slider(ax_slider, 'Slice', 0, 
                                self.image.GetSize()[2] - 1,
                                valinit=self.current_slice, valstep=1)
            self.slider.on_changed(self._on_slider_changed)
        
        plt.tight_layout(rect=[0, 0.08, 1, 1])
        plt.show()


# ============================================================================
# TIMESERIESPLOTTER - Time Series Plotter
# ============================================================================

class TimeSeriesPlotter(PlotViewer):
    """
    Interactive viewer for time-series images (TimeSeriesable).
    
    Features:
    - Frame navigation via slider
    - Frame number display
    - Temporal statistics
    - Per-frame overlay support
    - Playback simulation (frame stepping)
    
    Example
    -------
    >>> from pyable_eros_montin import TimeSeriesable
    >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
    >>> ts.plotOverlay()  # Interactive frame selector
    """
    
    def __init__(self, timeseries_image: Union[sitk.Image, 'TimeSeriesable'],
                 overlay: Optional[Union[sitk.Image, 'Imaginable']] = None,
                 title: str = "Time Series Viewer", **kwargs):
        """
        Initialize TimeSeriesPlotter.
        
        Parameters
        ----------
        timeseries_image : sitk.Image or TimeSeriesable
            4D time-series image
        overlay : sitk.Image or Imaginable, optional
            Overlay image (must be 3D)
        title : str, default="Time Series Viewer"
            Window title
        **kwargs
            Additional arguments for PlotViewer
        """
        # Convert TimeSeriesable to sitk.Image
        if hasattr(timeseries_image, 'getImage'):
            self.timeseriesable = timeseries_image
            timeseries_image = timeseries_image.getImage()
        else:
            self.timeseriesable = None
        
        if overlay is not None and hasattr(overlay, 'getImage'):
            overlay = overlay.getImage()
        
        # Validate 4D image
        if timeseries_image.GetDimension() != 4:
            raise ValueError(f"TimeSeriesPlotter requires 4D image, got {timeseries_image.GetDimension()}D")
        
        self.n_frames = timeseries_image.GetSize()[3]
        self.current_frame = 0
        self.timeseries_image = timeseries_image
        
        # Extract first frame as 3D scalar for initialization
        extractor = sitk.ExtractImageFilter()
        size = list(timeseries_image.GetSize())
        size[3] = 0
        extractor.SetSize(size)
        extractor.SetIndex([0, 0, 0, 0])
        
        frame_image = extractor.Execute(timeseries_image)
        
        super().__init__(frame_image, overlay, title=title, **kwargs)
        self.image = frame_image
    
    def _update_frame(self, frame_idx: int):
        """Update displayed frame."""
        self.current_frame = frame_idx
        
        extractor = sitk.ExtractImageFilter()
        size = list(self.timeseries_image.GetSize())
        size[3] = 0
        extractor.SetSize(size)
        extractor.SetIndex([0, 0, 0, frame_idx])
        
        self.image = extractor.Execute(self.timeseries_image)
        self._update_display()
    
    def _on_frame_slider_changed(self, val):
        """Callback for frame slider."""
        frame_idx = int(val)
        self._update_frame(frame_idx)
    
    def _update_display(self):
        """Update display with frame number."""
        self.ax.clear()
        
        # Get 2D slice from current frame
        img_array = self._get_2d_slice(self.image, self.current_slice)
        img_array = self._normalize_array(img_array)
        
        # Display
        self.ax.imshow(img_array, cmap=self.cmap, origin='lower')
        
        # Overlay
        if self.overlay is not None:
            overlay_array = self._get_2d_slice(self.overlay, self.current_slice)
            overlay_array = self._normalize_array(overlay_array)
            self.ax.imshow(overlay_array, cmap=self.cmap_overlay, 
                          alpha=self.alpha, origin='lower')
        
        title = f"{self.title} (Frame {self.current_frame}/{self.n_frames - 1})"
        if self.current_slice is not None:
            title += f", Slice {self.current_slice}"
        self.ax.set_title(title)
        
        self.fig.canvas.draw_idle()
    
    def show(self, frame: int = 0, alpha: float = 0.5, 
             slice_idx: Optional[int] = None):
        """
        Display time-series viewer with frame selector.
        
        Parameters
        ----------
        frame : int, default=0
            Initial frame to display
        alpha : float, default=0.5
            Overlay opacity
        slice_idx : int, optional
            Initial slice for 3D frames
        
        Example
        -------
        >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
        >>> plotter = TimeSeriesPlotter(ts)
        >>> plotter.show(frame=0)  # Start with first frame
        """
        self.alpha = alpha
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        
        # Update to initial frame
        self._update_frame(frame)
        
        if self.image.GetDimension() == 3:
            if slice_idx is None:
                slice_idx = self.image.GetSize()[2] // 2
            self.current_slice = slice_idx
        
        self._update_display()
        
        # Frame slider
        ax_frame = plt.axes([0.2, 0.10, 0.6, 0.03])
        frame_slider = Slider(ax_frame, 'Frame', 0, self.n_frames - 1,
                             valinit=frame, valstep=1)
        frame_slider.on_changed(self._on_frame_slider_changed)
        
        # Slice slider (if 3D frames)
        if self.image.GetDimension() == 3:
            ax_slice = plt.axes([0.2, 0.05, 0.6, 0.03])
            self.slider = Slider(ax_slice, 'Slice', 0, 
                                self.image.GetSize()[2] - 1,
                                valinit=self.current_slice, valstep=1)
            self.slider.on_changed(self._on_slider_changed)
        
        plt.tight_layout(rect=[0, 0.15, 1, 1])
        plt.show()


# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================

def plotOverlay(image: Union[sitk.Image, 'Imaginable', 'Vectorable', 'TimeSeriesable'],
                overlay: Optional[Union[sitk.Image, 'Imaginable', 'Vectorable', 'TimeSeriesable']] = None,
                alpha: float = 0.5,
                title: Optional[str] = None,
                component: Optional[int] = None,
                frame: Optional[int] = None,
                slice_idx: Optional[int] = None,
                show_vectors: bool = False,
                **kwargs) -> PlotViewer:
    """
    Convenience function for plotting images with overlays.
    
    Automatically selects appropriate viewer (Scalar, Vector, or TimeSeries)
    based on image type. Can be called directly or as method on able classes.
    
    Parameters
    ----------
    image : sitk.Image or Imaginable or Vectorable or TimeSeriesable
        Primary image to display
    overlay : same types, optional
        Overlay image (automatically resampled to match primary)
    alpha : float, default=0.5
        Overlay opacity (0-1)
    title : str, optional
        Figure title (auto-generated if None)
    component : int, optional
        Component to display for vector fields (0=X, 1=Y, 2=Z)
    frame : int, optional
        Frame to display for time-series (0-based)
    slice_idx : int, optional
        Slice index for 3D images (middle slice if None)
    show_vectors : bool, default=False
        Show vector arrows on magnitude image
    **kwargs
        Additional arguments passed to viewer
    
    Returns
    -------
    viewer : PlotViewer
        Viewer instance
    
    Examples
    --------
    **Scalar images:**
    
    >>> from pyable_eros_montin import Imaginable, plotOverlay
    >>> img = Imaginable('image.nii.gz')
    >>> overlay = Imaginable('segmentation.nii.gz')
    >>> plotOverlay(img, overlay, alpha=0.6)
    
    **Via method on able:**
    
    >>> img.plotOverlay(overlay, alpha=0.6)
    
    **Vector fields with component selection:**
    
    >>> from pyable_eros_montin import Vectorable
    >>> vf = Vectorable('displacement.mha')
    >>> vf.plotOverlay()  # Interactive component selector
    >>> vf.plotOverlay(component=0)  # Show X component directly
    
    **Time series with frame selection:**
    
    >>> from pyable_eros_montin import TimeSeriesable
    >>> ts = TimeSeriesable('cardiac_4d.nii.gz')
    >>> ts.plotOverlay()  # Interactive frame selector
    >>> ts.plotOverlay(frame=5)  # Show 6th frame directly
    """
    # Determine image type and create appropriate viewer
    is_vector = False
    is_timeseries = False
    
    if hasattr(image, 'getNumberOfComponents'):
        # Imaginable or Vectorable
        n_comp = image.getNumberOfComponents() if hasattr(image, 'getNumberOfComponents') else 1
        if n_comp > 1:
            is_vector = True
    
    if hasattr(image, 'getNumberOfFrames'):
        is_timeseries = True
    
    # Set default title
    if title is None:
        if is_vector:
            title = "Vector Field Viewer"
        elif is_timeseries:
            title = "Time Series Viewer"
        else:
            title = "Image Viewer"
    
    # Create appropriate viewer
    if is_timeseries:
        viewer = TimeSeriesPlotter(image, overlay, title=title, **kwargs)
        viewer.show(frame=frame or 0, alpha=alpha, slice_idx=slice_idx)
    
    elif is_vector:
        viewer = VectorPlotter(image, overlay, title=title, show_vectors=show_vectors, **kwargs)
        viewer.show(component=component or 0, alpha=alpha, slice_idx=slice_idx)
    
    else:
        viewer = ScalarPlotter(image, overlay, title=title, **kwargs)
        viewer.show(slice_idx=slice_idx, alpha=alpha)
    
    return viewer
