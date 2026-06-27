# Pyable Functionality Example Battery

This directory contains self-contained examples that exercise the main pyable APIs using synthetic 3D data.

Run the full suite:

```powershell
conda run -n able python examples\functionality_battery\run_all_functionality_examples.py
```

Run one example:

```powershell
conda run -n able python examples\functionality_battery\05_roiable_morphology_and_shape.py
```

Outputs are written to:

```text
examples\_outputs\functionality_battery
```

The suite covers:

- Synthetic dataset generation
- `Imaginable` image loading, metadata, NumPy conversion, statistics, and I/O
- Geometry edits, resampling, padding, cropping, and transforms
- Image math, filtering, edge maps, and sharpening
- Thresholding, Otsu, multi-Otsu, connected thresholding, and watershed segmentation
- `Roiable` morphology, shells, distance maps, shape metrics, and comparisons
- `LabelMapable` label extraction, label editing, priors, and ROI conversion
- Deformation field and transform utilities
- `Vectorable` components, magnitude, scaling, normalization, smoothing, and statistics
- `TimeSeriesable` frames, frame ranges, temporal mean/variance/std, and frame filtering
- Overlay PNG/report exports
- VTK/ParaView exports and SimpleITK/VTK conversion
- Batch directory processing with `processImageDirectory`