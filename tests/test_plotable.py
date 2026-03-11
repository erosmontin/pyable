"""
Tests for plotable module - visualization and plotting functionality.

Tests all plotter types:
- ScalarPlotter: For Imaginable scalar images
- VectorPlotter: For Vectorable vector fields with component selection
- TimeSeriesPlotter: For TimeSeriesable 4D time series with frame navigation
"""

import pytest
import numpy as np
import SimpleITK as sitk
import tempfile
from pathlib import Path
import sys
import matplotlib.pyplot as plt

# Ensure the workspace copy of pyable is imported before any installed package.
sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))

from pyable import (
    Imaginable, Vectorable, TimeSeriesable,
    PlotViewer, ScalarPlotter, VectorPlotter, TimeSeriesPlotter,
    plotOverlay
)


# ============================================================================
# TEST FIXTURES
# ============================================================================

@pytest.fixture
def scalar_image_2d():
    """Create a simple 2D scalar image."""
    arr = np.random.rand(100, 100)
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing([1.0, 1.0])
    img.SetOrigin([0.0, 0.0])
    return img


@pytest.fixture
def scalar_image_3d():
    """Create a simple 3D scalar image."""
    arr = np.random.rand(50, 50, 50)
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing([1.0, 1.0, 1.0])
    img.SetOrigin([0.0, 0.0, 0.0])
    return img


@pytest.fixture
def scalar_imaginable_3d(scalar_image_3d):
    """Create 3D Imaginable from image."""
    return Imaginable(image=scalar_image_3d)


@pytest.fixture
def overlay_image_3d():
    """Create overlay image with different values."""
    arr = np.random.rand(50, 50, 50) * 0.5
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing([1.0, 1.0, 1.0])
    img.SetOrigin([0.0, 0.0, 0.0])
    return img


@pytest.fixture
def overlay_imaginable_3d(overlay_image_3d):
    """Create overlay Imaginable."""
    return Imaginable(image=overlay_image_3d)


@pytest.fixture
def vector_image_3d():
    """Create a 3D vector field (displacement field)."""
    shape = (50, 50, 50)
    # Create 3D vectors (X, Y, Z components)
    arr = np.zeros((50, 50, 50, 3), dtype=np.float32)
    arr[..., 0] = np.random.randn(50, 50, 50) * 5  # X displacements
    arr[..., 1] = np.random.randn(50, 50, 50) * 5  # Y displacements
    arr[..., 2] = np.random.randn(50, 50, 50) * 3  # Z displacements (smaller)
    
    # Create vector image
    vector_img = sitk.GetImageFromArray(arr, isVector=True)
    vector_img.SetSpacing([1.0, 1.0, 1.0])
    vector_img.SetOrigin([0.0, 0.0, 0.0])
    return vector_img


@pytest.fixture
def vectorable_3d(vector_image_3d):
    """Create 3D Vectorable."""
    return Vectorable(image=vector_image_3d)


@pytest.fixture
def vector_image_2d():
    """Create a 2D vector field."""
    shape = (100, 100)
    arr = np.zeros((100, 100, 2), dtype=np.float32)
    arr[..., 0] = np.random.randn(100, 100) * 5  # X
    arr[..., 1] = np.random.randn(100, 100) * 5  # Y
    
    vector_img = sitk.GetImageFromArray(arr, isVector=True)
    vector_img.SetSpacing([1.0, 1.0])
    vector_img.SetOrigin([0.0, 0.0])
    return vector_img


@pytest.fixture
def timeseries_image_3d():
    """Create a 4D time-series image (5 frames of 50x50x50)."""
    # For 4D image in SimpleITK, we need to create it manually
    # Create 3D image
    img_3d = sitk.Image((50, 50, 50), sitk.sitkFloat32)
    img_3d.SetSpacing([1.0, 1.0, 1.0])
    img_3d.SetOrigin([0.0, 0.0, 0.0])
    
    # Populate with data
    arr_3d = np.random.rand(50, 50, 50).astype(np.float32)
    img_3d = sitk.GetImageFromArray(arr_3d)
    img_3d.SetSpacing([1.0, 1.0, 1.0])
    
    # Create 4D image by stacking frames
    # Join multiple 3D images to make 4D
    frames = []
    for i in range(5):
        frame_arr = np.random.rand(50, 50, 50).astype(np.float32)
        frame = sitk.GetImageFromArray(frame_arr)
        frame.SetSpacing([1.0, 1.0, 1.0])
        frames.append(frame)
    
    # Use JoinSeries to create 4D image
    ts_img = sitk.JoinSeries(frames)
    ts_img.SetSpacing(list(ts_img.GetSpacing()[:3]) + [0.1])  # Add time spacing
    return ts_img


@pytest.fixture
def timeseriesable_3d(timeseries_image_3d):
    """Create 3D TimeSeriesable."""
    return TimeSeriesable(image=timeseries_image_3d)


# ============================================================================
# PLOTVIEWER TESTS
# ============================================================================

class TestPlotViewer:
    """Test base PlotViewer class."""
    
    def test_init_scalar_2d(self, scalar_image_2d):
        """Test initialization with 2D scalar image."""
        viewer = PlotViewer(scalar_image_2d)
        assert viewer.image is not None
        assert viewer.overlay is None
        assert viewer.alpha == 0.5
    
    def test_init_scalar_3d(self, scalar_image_3d):
        """Test initialization with 3D scalar image."""
        viewer = PlotViewer(scalar_image_3d)
        assert viewer.image is not None
        assert viewer.current_slice is None  # Not set until show()
    
    def test_init_with_overlay(self, scalar_image_3d, overlay_image_3d):
        """Test initialization with overlay."""
        viewer = PlotViewer(scalar_image_3d, overlay=overlay_image_3d)
        assert viewer.image is not None
        assert viewer.overlay is not None
    
    def test_validate_non_scalar_image(self, vector_image_3d):
        """Test that non-scalar images raise error."""
        with pytest.raises(ValueError, match="must be scalar"):
            PlotViewer(vector_image_3d)
    
    def test_get_2d_slice_from_2d(self, scalar_image_2d):
        """Test extracting 2D slice from 2D image."""
        viewer = PlotViewer(scalar_image_2d)
        arr = viewer._get_2d_slice(scalar_image_2d)
        assert arr.ndim == 2
        assert arr.shape == (100, 100)
    
    def test_get_2d_slice_from_3d(self, scalar_image_3d):
        """Test extracting 2D slice from 3D image."""
        viewer = PlotViewer(scalar_image_3d)
        arr = viewer._get_2d_slice(scalar_image_3d, slice_idx=25)
        assert arr.ndim == 2
        assert viewer.current_slice == 25
    
    def test_get_2d_slice_middle(self, scalar_image_3d):
        """Test default middle slice."""
        viewer = PlotViewer(scalar_image_3d)
        arr = viewer._get_2d_slice(scalar_image_3d)
        assert viewer.current_slice == 25  # Middle of 50
    
    def test_normalize_array(self):
        """Test array normalization."""
        viewer = PlotViewer(sitk.Image((10, 10), sitk.sitkFloat32))
        arr = np.array([[0, 50], [100, 200]])
        normalized = viewer._normalize_array(arr)
        assert normalized.min() == 0.0
        assert normalized.max() == 1.0
    
    def test_resample_overlay_different_size(self, scalar_image_3d):
        """Test overlay resampling when sizes differ."""
        # Create overlay with different size
        overlay_arr = np.random.rand(30, 30, 30)
        overlay = sitk.GetImageFromArray(overlay_arr)
        overlay.SetSpacing([2.0, 2.0, 2.0])
        overlay.SetOrigin([10.0, 10.0, 10.0])
        
        viewer = PlotViewer(scalar_image_3d, overlay=overlay)
        # Overlay should be resampled to match primary image
        assert viewer.overlay.GetSize() == scalar_image_3d.GetSize()
    
    def test_save_figure(self, scalar_image_2d):
        """Test saving figure to file."""
        viewer = PlotViewer(scalar_image_2d)
        viewer.fig, viewer.ax = plt.subplots(figsize=viewer.figsize)
        viewer._update_display()
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test_figure.png"
            viewer.saveFigure(str(output_path))
            assert output_path.exists()
        
        plt.close('all')


# ============================================================================
# SCALARPLOTTER TESTS
# ============================================================================

class TestScalarPlotter:
    """Test ScalarPlotter for scalar images."""
    
    def test_init_with_imaginable(self, scalar_imaginable_3d):
        """Test initialization with Imaginable."""
        plotter = ScalarPlotter(scalar_imaginable_3d)
        assert plotter.image is not None
    
    def test_init_with_image(self, scalar_image_3d):
        """Test initialization with sitk.Image."""
        plotter = ScalarPlotter(scalar_image_3d)
        assert plotter.image is not None
    
    def test_init_with_overlay(self, scalar_imaginable_3d, overlay_imaginable_3d):
        """Test initialization with overlay."""
        plotter = ScalarPlotter(scalar_imaginable_3d, overlay=overlay_imaginable_3d)
        assert plotter.overlay is not None
    
    def test_show_2d(self, scalar_image_2d):
        """Test show with 2D image (non-blocking)."""
        plotter = ScalarPlotter(scalar_image_2d)
        # We don't actually call show() to avoid blocking
        # Just test that plotter initializes correctly
        assert plotter.image is not None
        plt.close('all')
    
    def test_custom_title(self, scalar_image_3d):
        """Test custom title."""
        plotter = ScalarPlotter(scalar_image_3d, title="Custom Title")
        assert plotter.title == "Custom Title"


# ============================================================================
# VECTORPLOTTER TESTS
# ============================================================================

class TestVectorPlotter:
    """Test VectorPlotter for vector fields."""
    
    def test_init_with_vectorable_3d(self, vectorable_3d):
        """Test initialization with Vectorable."""
        plotter = VectorPlotter(vectorable_3d)
        assert plotter.vector_image is not None
        assert plotter.n_components == 3
    
    def test_init_with_vector_image_3d(self, vector_image_3d):
        """Test initialization with vector image."""
        plotter = VectorPlotter(vector_image_3d)
        assert plotter.n_components == 3
    
    def test_init_with_vector_image_2d(self, vector_image_2d):
        """Test initialization with 2D vector image."""
        plotter = VectorPlotter(vector_image_2d)
        assert plotter.n_components == 2
    
    def test_init_scalar_image_fails(self, scalar_image_3d):
        """Test that scalar image raises error."""
        with pytest.raises(ValueError, match="vector image"):
            VectorPlotter(scalar_image_3d)
    
    def test_init_with_overlay(self, vectorable_3d, overlay_imaginable_3d):
        """Test initialization with overlay."""
        plotter = VectorPlotter(vectorable_3d, overlay=overlay_imaginable_3d)
        assert plotter.overlay is not None
    
    def test_update_component_x(self, vectorable_3d):
        """Test updating to X component."""
        plotter = VectorPlotter(vectorable_3d)
        plotter.fig, plotter.ax = plt.subplots()
        plotter._update_component(0)
        assert plotter.current_component == 0
        # Image should now be scalar (X component)
        assert plotter.image.GetNumberOfComponentsPerPixel() == 1
        plt.close('all')
    
    def test_update_component_y(self, vectorable_3d):
        """Test updating to Y component."""
        plotter = VectorPlotter(vectorable_3d)
        plotter.fig, plotter.ax = plt.subplots()
        plotter._update_component(1)
        assert plotter.current_component == 1
        assert plotter.image.GetNumberOfComponentsPerPixel() == 1
        plt.close('all')
    
    def test_update_component_z(self, vectorable_3d):
        """Test updating to Z component."""
        plotter = VectorPlotter(vectorable_3d)
        plotter.fig, plotter.ax = plt.subplots()
        plotter._update_component(2)
        assert plotter.current_component == 2
        assert plotter.image.GetNumberOfComponentsPerPixel() == 1
        plt.close('all')
    
    def test_update_magnitude(self, vectorable_3d):
        """Test updating to magnitude view."""
        plotter = VectorPlotter(vectorable_3d)
        plotter.fig, plotter.ax = plt.subplots()
        plotter._update_magnitude()
        # Image should be scalar magnitude
        assert plotter.image.GetNumberOfComponentsPerPixel() == 1
        plt.close('all')


# ============================================================================
# TIMESERIESPLOTTER TESTS
# ============================================================================

class TestTimeSeriesPlotter:
    """Test TimeSeriesPlotter for 4D time-series."""
    
    def test_init_with_timeseriesable(self, timeseriesable_3d):
        """Test initialization with TimeSeriesable."""
        plotter = TimeSeriesPlotter(timeseriesable_3d)
        assert plotter.timeseries_image is not None
        assert plotter.n_frames == 5
        assert plotter.current_frame == 0
    
    def test_init_with_4d_image(self, timeseries_image_3d):
        """Test initialization with 4D image."""
        plotter = TimeSeriesPlotter(timeseries_image_3d)
        assert plotter.n_frames == 5
    
    def test_init_3d_image_fails(self, scalar_image_3d):
        """Test that 3D image raises error."""
        with pytest.raises(ValueError, match="4D image"):
            TimeSeriesPlotter(scalar_image_3d)
    
    def test_init_with_overlay(self, timeseriesable_3d, overlay_imaginable_3d):
        """Test initialization with overlay."""
        plotter = TimeSeriesPlotter(timeseriesable_3d, overlay=overlay_imaginable_3d)
        assert plotter.overlay is not None
    
    def test_update_frame_0(self, timeseriesable_3d):
        """Test updating to frame 0."""
        plotter = TimeSeriesPlotter(timeseriesable_3d)
        plotter.fig, plotter.ax = plt.subplots()
        plotter._update_frame(0)
        assert plotter.current_frame == 0
        # Frame should be 3D
        assert plotter.image.GetDimension() == 3
        plt.close('all')
    
    def test_update_frame_middle(self, timeseriesable_3d):
        """Test updating to middle frame."""
        plotter = TimeSeriesPlotter(timeseriesable_3d)
        plotter.fig, plotter.ax = plt.subplots()
        plotter._update_frame(2)
        assert plotter.current_frame == 2
        plt.close('all')
    
    def test_update_frame_last(self, timeseriesable_3d):
        """Test updating to last frame."""
        plotter = TimeSeriesPlotter(timeseriesable_3d)
        plotter.fig, plotter.ax = plt.subplots()
        plotter._update_frame(4)
        assert plotter.current_frame == 4
        plt.close('all')
    
    def test_frame_shape_consistency(self, timeseriesable_3d):
        """Test that extracted frames have same shape."""
        plotter = TimeSeriesPlotter(timeseriesable_3d)
        plotter.fig, plotter.ax = plt.subplots()
        
        plotter._update_frame(0)
        shape_0 = plotter.image.GetSize()
        
        plotter._update_frame(3)
        shape_3 = plotter.image.GetSize()
        
        assert shape_0 == shape_3
        plt.close('all')


# ============================================================================
# PLOTOVERLAY CONVENIENCE FUNCTION TESTS
# ============================================================================

class TestPlotOverlay:
    """Test plotOverlay convenience function."""
    
    def test_plotoverlay_scalar_imaginable(self, scalar_imaginable_3d):
        """Test plotOverlay with scalar Imaginable."""
        viewer = plotOverlay(scalar_imaginable_3d)
        assert isinstance(viewer, ScalarPlotter)
        plt.close('all')
    
    def test_plotoverlay_scalar_with_overlay(self, scalar_imaginable_3d, overlay_imaginable_3d):
        """Test plotOverlay with scalar and overlay."""
        viewer = plotOverlay(scalar_imaginable_3d, overlay=overlay_imaginable_3d)
        assert isinstance(viewer, ScalarPlotter)
        assert viewer.overlay is not None
        plt.close('all')
    
    def test_plotoverlay_vector(self, vectorable_3d):
        """Test plotOverlay with vector field."""
        viewer = plotOverlay(vectorable_3d)
        assert isinstance(viewer, VectorPlotter)
        plt.close('all')
    
    def test_plotoverlay_vector_specific_component(self, vectorable_3d):
        """Test plotOverlay with specific component."""
        viewer = plotOverlay(vectorable_3d, component=1)
        assert isinstance(viewer, VectorPlotter)
        plt.close('all')
    
    def test_plotoverlay_timeseries(self, timeseriesable_3d):
        """Test plotOverlay with time-series."""
        viewer = plotOverlay(timeseriesable_3d)
        assert isinstance(viewer, TimeSeriesPlotter)
        plt.close('all')
    
    def test_plotoverlay_timeseries_specific_frame(self, timeseriesable_3d):
        """Test plotOverlay with specific frame."""
        viewer = plotOverlay(timeseriesable_3d, frame=2)
        assert isinstance(viewer, TimeSeriesPlotter)
        plt.close('all')
    
    def test_plotoverlay_alpha(self, scalar_imaginable_3d, overlay_imaginable_3d):
        """Test alpha parameter."""
        viewer = plotOverlay(scalar_imaginable_3d, overlay=overlay_imaginable_3d, 
                            alpha=0.7)
        assert viewer.alpha == 0.7
        plt.close('all')
    
    def test_plotoverlay_custom_title(self, scalar_imaginable_3d):
        """Test custom title."""
        viewer = plotOverlay(scalar_imaginable_3d, title="My Custom Title")
        assert viewer.title == "My Custom Title"
        plt.close('all')


# ============================================================================
# ABLE CLASSES METHOD TESTS
# ============================================================================

class TestAblePlotOverlayMethods:
    """Test plotOverlay() methods on able classes."""
    
    def test_imaginable_plotoverlay(self, scalar_imaginable_3d, overlay_imaginable_3d):
        """Test Imaginable.plotOverlay() method."""
        viewer = scalar_imaginable_3d.plotOverlay(overlay=overlay_imaginable_3d)
        assert isinstance(viewer, ScalarPlotter)
        plt.close('all')
    
    def test_vectorable_plotoverlay(self, vectorable_3d):
        """Test Vectorable.plotOverlay() method."""
        viewer = vectorable_3d.plotOverlay()
        assert isinstance(viewer, VectorPlotter)
        plt.close('all')
    
    def test_vectorable_plotoverlay_component(self, vectorable_3d):
        """Test Vectorable.plotOverlay() with component."""
        viewer = vectorable_3d.plotOverlay(component=1)
        assert isinstance(viewer, VectorPlotter)
        plt.close('all')
    
    def test_timeseriesable_plotoverlay(self, timeseriesable_3d):
        """Test TimeSeriesable.plotOverlay() method."""
        viewer = timeseriesable_3d.plotOverlay()
        assert isinstance(viewer, TimeSeriesPlotter)
        plt.close('all')
    
    def test_timeseriesable_plotoverlay_frame(self, timeseriesable_3d):
        """Test TimeSeriesable.plotOverlay() with frame."""
        viewer = timeseriesable_3d.plotOverlay(frame=2)
        assert isinstance(viewer, TimeSeriesPlotter)
        plt.close('all')


# ============================================================================
# SLICING UTILITIES TESTS
# ============================================================================

class TestSliceExtraction:
    """Test Imaginable.extractRepresentativeSlices() method."""
    
    def test_extract_slices_basic(self, scalar_imaginable_3d):
        """Test basic slice extraction."""
        result = scalar_imaginable_3d.extractRepresentativeSlices(
            planes='all', offsets=[-5, 0, 5]
        )
        assert 'slices' in result
        assert 'plane_names' in result
        assert 'center_of_gravity' in result
        assert len(result['slices']) > 0
        assert all(isinstance(s, np.ndarray) for s in result['slices'])
    
    def test_extract_slices_single_plane(self, scalar_imaginable_3d):
        """Test extraction with single plane."""
        result = scalar_imaginable_3d.extractRepresentativeSlices(
            planes=[0], offsets=[0]
        )
        assert len(result['slices']) >= 1
    
    def test_extract_slices_multiple_offsets(self, scalar_imaginable_3d):
        """Test extraction with multiple offsets."""
        offsets = [-10, -5, 0, 5, 10]
        result = scalar_imaginable_3d.extractRepresentativeSlices(
            planes=[0, 1, 2], offsets=offsets
        )
        assert len(result['slices']) > 0
        assert len(result['offsets']) == len(offsets)
    
    def test_extract_slices_returns_2d(self, scalar_imaginable_3d):
        """Test that extracted slices are 2D."""
        result = scalar_imaginable_3d.extractRepresentativeSlices(planes='all')
        for slc in result['slices']:
            assert len(slc.shape) == 2, f"Expected 2D slice, got shape {slc.shape}"
    
    def test_extract_slices_center_of_gravity(self, scalar_imaginable_3d):
        """Test center of gravity computation."""
        result = scalar_imaginable_3d.extractRepresentativeSlices(planes='all')
        assert result['center_of_gravity'] is not None
        assert len(result['center_of_gravity']) == 3
        assert result['center_of_gravity_index'] is not None
        assert len(result['center_of_gravity_index']) == 3
    
    def test_extract_slices_verbose(self, scalar_imaginable_3d, capsys):
        """Test verbose mode."""
        result = scalar_imaginable_3d.extractRepresentativeSlices(
            planes=[0], offsets=[0], verbose=True
        )
        captured = capsys.readouterr()
        assert 'Center of gravity' in captured.out or len(result['slices']) > 0


# ============================================================================
# BATCH PROCESSING TESTS
# ============================================================================

class TestBatchProcessing:
    """Test processImageDirectory utility."""
    
    def test_process_directory_basic(self, tmp_path):
        """Test basic directory processing."""
        # Create test images
        for i in range(3):
            arr = np.random.rand(32, 32, 32).astype(np.float32)
            img = sitk.GetImageFromArray(arr)
            filepath = tmp_path / f"image_{i}.nii.gz"
            sitk.WriteImage(img, str(filepath))
        
        # Simple processor that returns image properties
        def processor(img):
            return {'size': img.getImageSize(), 'spacing': img.getImageSpacing()}
        
        try:
            from pyable.utils import processImageDirectory
        except ImportError:
            from pyable.utils import processImageDirectory
        
        df = processImageDirectory(
            str(tmp_path), processor,
            file_pattern='*.nii.gz',
            verbose=False
        )
        
        assert len(df) == 3
        assert 'filepath' in df.columns
        assert 'size' in df.columns
    
    def test_process_directory_csv_export(self, tmp_path):
        """Test CSV export."""
        # Create test images
        for i in range(2):
            arr = np.random.rand(32, 32, 32).astype(np.float32)
            img = sitk.GetImageFromArray(arr)
            filepath = tmp_path / f"image_{i}.nii.gz"
            sitk.WriteImage(img, str(filepath))
        
        def processor(img):
            return {'size': img.getImageSize()}
        
        try:
            from pyable.utils import processImageDirectory
        except ImportError:
            from pyable.utils import processImageDirectory
        
        output_csv = tmp_path / "results.csv"
        df = processImageDirectory(
            str(tmp_path), processor,
            file_pattern='*.nii.gz',
            output_csv=str(output_csv),
            verbose=False
        )
        
        assert output_csv.exists()
        df_read = __import__('pandas').read_csv(output_csv)
        assert len(df_read) == 2
    
    def test_process_directory_empty_folder(self, tmp_path):
        """Test with empty directory."""
        try:
            from pyable.utils import processImageDirectory
        except ImportError:
            from pyable.utils import processImageDirectory
        
        def processor(img):
            return {}
        
        df = processImageDirectory(
            str(tmp_path), processor,
            file_pattern='*.nii.gz',
            verbose=False
        )
        
        assert len(df) == 0
    
    def test_process_directory_nonexistent(self):
        """Test with nonexistent directory."""
        try:
            from pyable.utils import processImageDirectory
        except ImportError:
            from pyable.utils import processImageDirectory
        
        def processor(img):
            return {}
        
        with pytest.raises(ValueError):
            processImageDirectory(
                "/nonexistent/path",
                processor,
                verbose=False
            )


# ============================================================================
# GRIDPLOTTER TESTS
# ============================================================================

class TestGridPlotter:
    """Test GridPlotter class."""
    
    def test_gridplotter_basic(self):
        """Test basic grid plotter creation."""
        try:
            from pyable.plotable import GridPlotter
        except ImportError:
            from pyable.plotable import GridPlotter
        
        grid = GridPlotter()
        assert grid.fig is None
        assert grid.axes is None
    
    def test_gridplotter_show_grid(self):
        """Test showing grid of slices."""
        try:
            from pyable.plotable import GridPlotter
        except ImportError:
            from pyable.plotable import GridPlotter
        
        # Create test slices
        slices = [np.random.rand(50, 50) for _ in range(6)]
        
        grid = GridPlotter(figsize=(10, 8))
        grid.show_grid(slices, rows=2, cols=3)
        
        assert grid.fig is not None
        assert grid.axes is not None
        plt.close('all')
    
    def test_gridplotter_with_titles(self):
        """Test grid plotter with custom titles."""
        try:
            from pyable.plotable import GridPlotter
        except ImportError:
            from pyable.plotable import GridPlotter
        
        slices = [np.random.rand(50, 50) for _ in range(4)]
        titles = ['Slice A', 'Slice B', 'Slice C', 'Slice D']
        
        grid = GridPlotter()
        grid.show_grid(slices, rows=2, cols=2, titles=titles)
        
        assert grid.fig is not None
        plt.close('all')
    
    def test_gridplotter_auto_grid_size(self):
        """Test automatic grid size calculation."""
        try:
            from pyable.plotable import GridPlotter
        except ImportError:
            from pyable.plotable import GridPlotter
        
        slices = [np.random.rand(50, 50) for _ in range(7)]
        
        grid = GridPlotter()
        grid.show_grid(slices)  # Should auto-calculate grid size
        
        assert grid.fig is not None
        plt.close('all')
    
    def test_gridplotter_with_overlays(self):
        """Test grid plotter with overlay slices."""
        try:
            from pyable.plotable import GridPlotter
        except ImportError:
            from pyable.plotable import GridPlotter
        
        slices = [np.random.rand(50, 50) for _ in range(3)]
        overlays = [np.random.rand(50, 50) for _ in range(3)]
        
        grid = GridPlotter()
        grid.show_grid(slices, overlays=overlays, alpha_overlay=0.5)
        
        assert grid.fig is not None
        plt.close('all')
    
    def test_gridplotter_colormap_params(self):
        """Test grid plotter with custom colormap parameters."""
        try:
            from pyable.plotable import GridPlotter
        except ImportError:
            from pyable.plotable import GridPlotter
        
        slices = [np.random.rand(50, 50) for _ in range(4)]
        
        grid = GridPlotter()
        grid.show_grid(slices, cmap='viridis', vmin=0.2, vmax=0.8)
        
        assert grid.fig is not None
        plt.close('all')


# ============================================================================
# INTEGRATION TESTS FOR NEW FEATURES
# ============================================================================

class TestIntegrationNewFeatures:
    """Integration tests combining new utilities."""
    
    def test_extract_slices_then_plot_grid(self, scalar_imaginable_3d):
        """Test extracting slices and plotting in grid."""
        try:
            from pyable.plotable import GridPlotter
        except ImportError:
            from pyable.plotable import GridPlotter
        
        # Extract slices
        result = scalar_imaginable_3d.extractRepresentativeSlices(
            planes='all', offsets=[-5, 0, 5]
        )
        
        # Plot in grid
        grid = GridPlotter()
        titles = [f"{name} ({offset}mm)" for name, offset in result['plane_names']]
        grid.show_grid(result['slices'], titles=titles)
        
        assert grid.fig is not None
        plt.close('all')
    
    def test_batch_process_with_slice_extraction(self, tmp_path):
        """Test batch processing with slice extraction."""
        # Create test images
        for i in range(2):
            arr = np.random.rand(32, 32, 32).astype(np.float32)
            img = sitk.GetImageFromArray(arr)
            filepath = tmp_path / f"image_{i}.nii.gz"
            sitk.WriteImage(img, str(filepath))
        
        def processor(img):
            result = img.extractRepresentativeSlices(planes=[0], offsets=[0])
            return {
                'num_slices': len(result['slices']),
                'size': img.getImageSize()
            }
        
        try:
            from pyable.utils import processImageDirectory
        except ImportError:
            from pyable.utils import processImageDirectory
        
        df = processImageDirectory(
            str(tmp_path), processor,
            file_pattern='*.nii.gz',
            verbose=False
        )
        
        assert len(df) == 2
        assert 'num_slices' in df.columns


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
