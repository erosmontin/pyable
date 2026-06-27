# ParaView Example Battery

These examples export the same hip image, ROI, and label map in several VTK formats for ParaView.

Run one example:

```powershell
conda run -n able python examples\paraview_battery\03_roi_surface_smoothing_levels.py
```

Run the full battery:

```powershell
conda run -n able python examples\paraview_battery\run_all_paraview_examples.py
```

Outputs are written to:

```text
examples\_outputs\paraview_battery
```

The examples cover:

- Image volume export as `.vts`
- ROI surface export as `.vtp`
- ROI smoothing comparisons
- ROI voxel-mask volume export as `.vts`
- Label-map per-label surface export
- Label-map smoothing comparisons
- Label-map multiblock export as `.vtm`
- LPS, RAS, FSL, and raw index coordinate spaces