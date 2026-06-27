from _common import make_image, make_labelmap, make_roi, output_dir, print_created


image = make_image()
roi = make_roi()
labelmap = make_labelmap()
out = output_dir("11_plotting_and_overlay_exports")

single_overlay_path = out / "single_roi_overlay.png"
grid_overlay_path = out / "labelmap_overlay_grid.png"
report_path = out / "overlay_report.png"

image.overlayAbleImage(
    roi,
    axis=2,
    index=image.getImageCenterIndex()[2],
    title="Synthetic ROI overlay",
    save=str(single_overlay_path),
)

image.overlayAbleImage(
    labelmap,
    axis=2,
    index=image.getImageCenterIndex()[2],
    slice_offsets=[-6, 0, 6],
    title="Synthetic labelmap overlay grid",
    save=str(grid_overlay_path),
)

image.overlayReport(
    labelmap,
    save=str(report_path),
    show=False,
    title="Synthetic Overlay Report",
)

print_created(single_overlay_path, grid_overlay_path, report_path)