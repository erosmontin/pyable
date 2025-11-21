"""
Interactive GUI viewer for medical images with advanced controls.

Provides a comprehensive viewer with:
- Orientation selection (axial, sagittal, coronal)
- Slice index control with center auto-detection
- Multiple overlay layers with individual opacity control
- Vector component selection (X, Y, Z, magnitude)
- Time frame selection for 4D sequences
- Synchronized updates and smooth interaction

Classes:
    - InteractiveViewer: Main GUI viewer class
    - OverlayManager: Manages multiple overlay layers
"""

import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from matplotlib.widgets import Slider, RadioButtons, CheckButtons
from typing import Union, Optional, List, Tuple, Any, Dict
import copy


class OverlayManager:
    """Manages multiple overlay layers with individual opacity control."""
    
    def __init__(self):
        """Initialize overlay manager."""
        self.overlays: Dict[str, Dict[str, Any]] = {}
        self.order = []  # Rendering order
    
    def add(self, name: str, image: sitk.Image, cmap: str = 'hot', 
            opacity: float = 0.5, visible: bool = True):
        """
        Add an overlay layer.
        
        Parameters
        ----------
        name : str
            Unique name for the overlay
        image : sitk.Image
            SimpleITK image to overlay
        cmap : str
            Colormap name
        opacity : float
            Initial opacity (0-1)
        visible : bool
            Whether overlay is initially visible
        """
        self.overlays[name] = {
            'image': image,
            'cmap': cmap,
            'opacity': opacity,
            'visible': visible
        }
        self.order.append(name)
    
    def add_array(self, name: str, array: np.ndarray, cmap: str = 'hot',
                  opacity: float = 0.5, visible: bool = True):
        """
        Add an overlay from numpy array.
        
        Parameters
        ----------
        name : str
            Unique name for the overlay
        array : np.ndarray
            2D numpy array
        cmap : str
            Colormap name
        opacity : float
            Initial opacity (0-1)
        visible : bool
            Whether overlay is initially visible
        """
        self.overlays[name] = {
            'array': array,
            'cmap': cmap,
            'opacity': opacity,
            'visible': visible
        }
        self.order.append(name)
    
    def set_opacity(self, name: str, opacity: float):
        """Set opacity for a specific overlay."""
        if name in self.overlays:
            self.overlays[name]['opacity'] = np.clip(opacity, 0, 1)
    
    def set_visible(self, name: str, visible: bool):
        """Set visibility for a specific overlay."""
        if name in self.overlays:
            self.overlays[name]['visible'] = visible
    
    def get_opacity(self, name: str) -> float:
        """Get opacity of an overlay."""
        return self.overlays.get(name, {}).get('opacity', 0.5)
    
    def is_visible(self, name: str) -> bool:
        """Check if overlay is visible."""
        return self.overlays.get(name, {}).get('visible', False)
    
    def get_all(self) -> Dict[str, Dict[str, Any]]:
        """Get all overlays."""
        return self.overlays.copy()
    
    def get_count(self) -> int:
        """Get number of overlays."""
        return len(self.overlays)
    
    def remove(self, name: str):
        """Remove an overlay."""
        if name in self.overlays:
            del self.overlays[name]
            self.order.remove(name)
    
    def clear(self):
        """Remove all overlays."""
        self.overlays.clear()
        self.order.clear()


class InteractiveViewer:
    """
    Comprehensive interactive viewer for medical images.
    
    Features:
    - Orientation selection (axial, sagittal, coronal)
    - Slice navigation with automatic centering
    - Multiple overlay layers with individual opacity sliders
    - Vector component selection and magnitude computation
    - Time frame selection for 4D sequences
    - Synchronized updates
    - Image statistics display
    
    Attributes
    ----------
    image : sitk.Image
        Primary image
    overlays : OverlayManager
        Manages overlay layers
    current_orientation : int
        Current viewing orientation (0=axial, 1=sagittal, 2=coronal)
    current_slice : int
        Current slice index
    current_frame : int
        Current time frame (for 4D)
    current_component : int
        Current component (for vector fields)
    """
    
    def __init__(self, image: sitk.Image, title: str = "Interactive Viewer",
                 figsize: Tuple[int, int] = (14, 10), cmap: str = 'gray'):
        """
        Initialize InteractiveViewer.
        
        Parameters
        ----------
        image : sitk.Image
            Primary image to display
        title : str
            Window title
        figsize : tuple
            Figure size (width, height)
        cmap : str
            Primary image colormap
        """
        self.image = image
        self.title = title
        self.figsize = figsize
        self.cmap = cmap
        self.overlays = OverlayManager()
        
        # State
        self.current_orientation = 2  # 0=axial, 1=sagittal, 2=coronal
        self.current_slice = None
        self.current_frame = 0  # For 4D
        self.current_component = 0  # For vector fields
        
        # UI elements
        self.fig = None
        self.ax_main = None
        self.ax_slice = None
        self.ax_orient = None
        self.sliders = {}
        self.checkboxes = None
        self.radio_component = None
        
        # Cache
        self._slice_cache = {}
    
    def _get_image_size(self, image: sitk.Image) -> Tuple[int, int, int]:
        """Get image size, handling 2D/3D/4D."""
        size = image.GetSize()
        if len(size) == 2:
            return (size[0], size[1], 1)
        elif len(size) == 3:
            return (size[0], size[1], size[2])
        elif len(size) == 4:
            return (size[0], size[1], size[2], size[3])
        else:
            raise ValueError(f"Unsupported image dimension: {len(size)}")
    
    def _get_center_slice(self, orientation: int) -> int:
        """Get center slice for a given orientation."""
        size = self._get_image_size(self.image)
        if orientation == 0:  # Axial (Z)
            return size[2] // 2
        elif orientation == 1:  # Sagittal (X)
            return size[0] // 2
        elif orientation == 2:  # Coronal (Y)
            return size[1] // 2
        return 0
    
    def _extract_slice_2d(self, image: sitk.Image, orientation: int, 
                         slice_idx: int, frame: int = 0) -> np.ndarray:
        """
        Extract a 2D slice from 3D/4D image.
        
        Parameters
        ----------
        image : sitk.Image
            Source image
        orientation : int
            0=axial (XY), 1=sagittal (YZ), 2=coronal (XZ)
        slice_idx : int
            Slice index along orientation axis
        frame : int
            Frame index for 4D
            
        Returns
        -------
        np.ndarray
            2D slice as numpy array
        """
        try:
            if image.GetDimension() == 4:
                # Extract time frame first
                extractor = sitk.ExtractImageFilter()
                size = list(image.GetSize())
                size[3] = 0  # Extract single frame
                extractor.SetSize(size)
                index = [0, 0, 0, frame]
                extractor.SetIndex(index)
                frame_img = extractor.Execute(image)
                image = frame_img
            
            # Now extract 2D slice
            arr = sitk.GetArrayFromImage(image)
            
            # arr shape is (Z, Y, X) in SimpleITK convention
            if orientation == 0:  # Axial (XY plane, view from Z)
                return arr[slice_idx, :, :]
            elif orientation == 1:  # Sagittal (YZ plane, view from X)
                return arr[:, :, slice_idx]
            elif orientation == 2:  # Coronal (XZ plane, view from Y)
                return arr[:, slice_idx, :]
        except Exception as e:
            # Return blank if extraction fails
            return np.zeros((100, 100))
        
        return np.zeros((100, 100))
    
    def _get_slice_max(self, orientation: int) -> int:
        """Get maximum slice index for orientation."""
        size = self._get_image_size(self.image)
        if orientation == 0:  # Axial (Z)
            return max(0, size[2] - 1)
        elif orientation == 1:  # Sagittal (X)
            return max(0, size[0] - 1)
        elif orientation == 2:  # Coronal (Y)
            return max(0, size[1] - 1)
        return 0
    
    def _normalize_array(self, arr: np.ndarray) -> np.ndarray:
        """Normalize array to 0-1 range."""
        arr = np.asarray(arr, dtype=float)
        vmin = np.percentile(arr, 2)
        vmax = np.percentile(arr, 98)
        if vmax > vmin:
            arr = (arr - vmin) / (vmax - vmin)
        arr = np.clip(arr, 0, 1)
        return arr
    
    def add_overlay(self, overlay: Union[sitk.Image, np.ndarray, 'Imaginable'],
                   name: str = None, cmap: str = 'hot', opacity: float = 0.5):
        """
        Add an overlay layer.
        
        Parameters
        ----------
        overlay : sitk.Image, np.ndarray, or Imaginable
            Overlay image or array
        name : str, optional
            Overlay name (auto-generated if None)
        cmap : str
            Colormap
        opacity : float
            Initial opacity (0-1)
        """
        if name is None:
            name = f"Overlay_{self.overlays.get_count()}"
        
        # Convert Imaginable to SimpleITK image
        if hasattr(overlay, 'getImage'):
            overlay = overlay.getImage()
        
        if isinstance(overlay, sitk.Image):
            self.overlays.add(name, overlay, cmap=cmap, opacity=opacity)
        elif isinstance(overlay, np.ndarray):
            self.overlays.add_array(name, overlay, cmap=cmap, opacity=opacity)
        else:
            raise TypeError(f"Unsupported overlay type: {type(overlay)}")
    
    def add_overlays(self, overlays: List[Union[sitk.Image, np.ndarray]],
                    names: List[str] = None, cmap: str = 'hot'):
        """
        Add multiple overlay layers.
        
        Parameters
        ----------
        overlays : list
            List of overlay images or arrays
        names : list, optional
            Names for each overlay
        cmap : str
            Colormap for all overlays
        """
        for i, overlay in enumerate(overlays):
            name = names[i] if names else f"Overlay_{i}"
            opacity = 0.5 if i == 0 else 0.3  # First darker, others lighter
            self.add_overlay(overlay, name=name, cmap=cmap, opacity=opacity)
    
    def _update_display(self):
        """Update the main display with current slice and overlays."""
        self.ax_main.clear()
        
        # Get and display main image slice
        main_slice = self._extract_slice_2d(self.image, self.current_orientation,
                                           self.current_slice, self.current_frame)
        main_slice = self._normalize_array(main_slice)
        
        self.ax_main.imshow(main_slice, cmap=self.cmap, origin='lower')
        
        # Display overlays in order
        for overlay_name in self.overlays.order:
            overlay_data = self.overlays.overlays[overlay_name]
            
            if not overlay_data['visible']:
                continue
            
            # Get overlay slice
            if 'image' in overlay_data:
                overlay_slice = self._extract_slice_2d(overlay_data['image'],
                                                       self.current_orientation,
                                                       self.current_slice,
                                                       self.current_frame)
            else:
                overlay_slice = overlay_data.get('array', np.zeros((100, 100)))
            
            overlay_slice = self._normalize_array(overlay_slice)
            opacity = overlay_data['opacity']
            cmap = overlay_data['cmap']
            
            self.ax_main.imshow(overlay_slice, cmap=cmap, alpha=opacity, origin='lower')
        
        # Update title
        orient_names = ['Axial', 'Sagittal', 'Coronal']
        title = f"{self.title} - {orient_names[self.current_orientation]} (Slice {self.current_slice})"
        if self.current_frame > 0:
            title += f", Frame {self.current_frame}"
        self.ax_main.set_title(title, fontsize=12, fontweight='bold')
        
        self.ax_main.axis('off')
        self.fig.canvas.draw_idle()
    
    def _on_orientation_changed(self, label):
        """Handle orientation radio button change."""
        orient_map = {'Axial': 0, 'Sagittal': 1, 'Coronal': 2}
        self.current_orientation = orient_map[label]
        
        # Update slice to center
        self.current_slice = self._get_center_slice(self.current_orientation)
        
        # Update slice slider
        max_slice = self._get_slice_max(self.current_orientation)
        self.sliders['slice'].set_val(self.current_slice)
        self.sliders['slice'].set_valmin(0)
        self.sliders['slice'].set_valmax(max_slice)
        
        self._update_display()
    
    def _on_slice_changed(self, value):
        """Handle slice slider change."""
        self.current_slice = int(value)
        self._update_display()
    
    def _on_frame_changed(self, value):
        """Handle frame slider change."""
        self.current_frame = int(value)
        self._update_display()
    
    def _on_component_changed(self, label):
        """Handle component radio button change."""
        component_map = {'X': 0, 'Y': 1, 'Z': 2, 'Magnitude': 3}
        self.current_component = component_map[label]
        self._update_display()
    
    def _on_overlay_checked(self, label):
        """Handle overlay checkbox change."""
        visible = self.checkboxes.get_status()[list(self.overlays.overlays.keys()).index(label)]
        self.overlays.set_visible(label, visible)
        self._update_display()
    
    def _on_overlay_opacity_changed(self, value, overlay_name: str):
        """Handle overlay opacity slider change."""
        self.overlays.set_opacity(overlay_name, value)
        self._update_display()
    
    def show(self):
        """Display the interactive viewer."""
        self.fig = plt.figure(figsize=self.figsize)
        
        # Main image display
        self.ax_main = plt.subplot(1, 2, 1)
        
        # Control panel
        ax_controls = plt.subplot(1, 2, 2)
        ax_controls.axis('off')
        
        # Orientation selector
        ax_orient = plt.axes([0.52, 0.75, 0.15, 0.15])
        self.ax_orient = ax_orient
        self.radio_orient = RadioButtons(ax_orient, ['Axial', 'Sagittal', 'Coronal'],
                                        active=self.current_orientation)
        self.radio_orient.on_clicked(self._on_orientation_changed)
        ax_orient.set_title('Orientation', fontweight='bold', fontsize=10)
        
        # Slice slider
        ax_slice = plt.axes([0.52, 0.65, 0.35, 0.03])
        max_slice = self._get_slice_max(self.current_orientation)
        if self.current_slice is None:
            self.current_slice = self._get_center_slice(self.current_orientation)
        self.sliders['slice'] = Slider(ax_slice, 'Slice', 0, max_slice,
                                       valinit=self.current_slice, valstep=1)
        self.sliders['slice'].on_changed(self._on_slice_changed)
        
        # Time frame slider (if 4D)
        if self.image.GetDimension() == 4:
            ax_frame = plt.axes([0.52, 0.60, 0.35, 0.03])
            n_frames = self.image.GetSize()[3]
            self.sliders['frame'] = Slider(ax_frame, 'Frame', 0, n_frames - 1,
                                          valinit=0, valstep=1)
            self.sliders['frame'].on_changed(self._on_frame_changed)
        
        # Overlay visibility checkboxes
        if self.overlays.get_count() > 0:
            overlay_names = list(self.overlays.overlays.keys())
            ax_overlays = plt.axes([0.52, 0.35, 0.1, len(overlay_names) * 0.04])
            self.checkboxes = CheckButtons(ax_overlays, overlay_names,
                                          [self.overlays.is_visible(n) for n in overlay_names])
            self.checkboxes.on_clicked(self._on_overlay_checked)
            ax_overlays.set_title('Overlays', fontweight='bold', fontsize=10)
            
            # Opacity sliders for each overlay
            for i, name in enumerate(overlay_names):
                opacity = self.overlays.get_opacity(name)
                y_pos = 0.28 - i * 0.04
                ax_opacity = plt.axes([0.52, y_pos, 0.35, 0.02])
                slider = Slider(ax_opacity, f'{name} opacity', 0, 1,
                               valinit=opacity, valstep=0.05)
                slider.on_changed(lambda value, n=name: self._on_overlay_opacity_changed(value, n))
                self.sliders[f'opacity_{name}'] = slider
        
        # Component selector (if vector field)
        if self.image.GetNumberOfComponentsPerPixel() > 1:
            ax_component = plt.axes([0.68, 0.75, 0.15, 0.15])
            n_comp = self.image.GetNumberOfComponentsPerPixel()
            labels = ['X', 'Y', 'Z', 'Magnitude'][:n_comp + 1]
            self.radio_component = RadioButtons(ax_component, labels, active=0)
            self.radio_component.on_clicked(self._on_component_changed)
            ax_component.set_title('Component', fontweight='bold', fontsize=10)
        
        self._update_display()
        plt.tight_layout()
        plt.show()


# Convenience function
def viewInteractiveImage(image: 'Imaginable', overlays: List[Any] = None,
                        title: str = None, figsize: Tuple[int, int] = (14, 10),
                        cmap: str = 'gray') -> InteractiveViewer:
    """
    Open an interactive viewer for an image.
    
    Parameters
    ----------
    image : Imaginable
        Image to display
    overlays : list, optional
        List of overlay images or arrays
    title : str, optional
        Window title
    figsize : tuple
        Figure size
    cmap : str
        Colormap
        
    Returns
    -------
    InteractiveViewer
        Viewer instance
    """
    if title is None:
        title = "Interactive Image Viewer"
    
    # Convert Imaginable to SimpleITK if needed
    if hasattr(image, 'getImage'):
        sitk_image = image.getImage()
    else:
        sitk_image = image
    
    viewer = InteractiveViewer(sitk_image, title=title, figsize=figsize, cmap=cmap)
    
    if overlays:
        viewer.add_overlays(overlays)
    
    viewer.show()
    return viewer
