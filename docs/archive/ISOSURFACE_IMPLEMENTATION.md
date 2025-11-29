# Archived: Isosurface Implementation

This archived file mirrors the earlier top-level `ISOSURFACE_IMPLEMENTATION.md` and documents the VTK-based isosurface implementation used by `pyable.meshable` and `pyable.imaginable.renderIsosurface`.

Summary:

- Conversion helpers: `sitk2vtk()` and `vtk2sitk()` in `meshable.py`.
- Isosurface generator: `Imaginable.renderIsosurface()` uses VTK's marching cubes through a sitk->vtk pipeline.
- Features: scalar isosurfaces for continuous images, labeled ROI boundary extraction (label-preserving NN resampling before conversion), vector magnitude isosurface generation, color/opacity parameterization, batch rendering mode (`show=False`) returning VTK actors.
- Tests: `test_isosurface.py` covers scalar, ROI, vector magnitude and color modes.

Notes:
- See `docs/ISOSURFACE_RENDERING_GUIDE.md` for user-facing guidance and examples.
