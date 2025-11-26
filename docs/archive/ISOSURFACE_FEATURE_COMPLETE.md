# Archived: Isosurface Feature Complete Checklist

This archive preserves the checklist that was used during implementation and validation of the isosurface rendering work.

Checklist:

- [x] sitk2vtk / vtk2sitk conversion
- [x] Imaginable.renderIsosurface for scalar images
- [x] ROI label-preserving isosurface generation
- [x] Vector magnitude isosurface generation
- [x] Color/opacity support
- [x] show=False batch mode returning vtk actors
- [x] Tests: `test_isosurface.py` passing
- [x] Examples: `example_isosurface.py`, `example_show_roi_vtk.py`

Validation notes:
- Tests passed in the working environment. Examples run in non-interactive mode.
- High-polygon meshes observed when running on aparc+aseg volumes (normal for whole-brain segmentation). Consider adding decimation or multi-resolution exporting if file sizes become problematic.
