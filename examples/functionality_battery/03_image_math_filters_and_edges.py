from _common import image_copy, make_image, output_dir, print_created, write_json


image = make_image()
out = output_dir("03_image_math_filters_and_edges")

scaled = image_copy(image).add(25).multiply(0.5)
scaled_path = scaled.writeImageAs(str(out / "add_then_multiply.nii.gz"))

smoothed = image.smoothAnisotropic(iterations=3, time_step=0.04, conductance=2.0)
smoothed_path = smoothed.writeImageAs(str(out / "anisotropic_smoothed.nii.gz"))

edges = image.getEdgeMap(sigma=1.0)
edges_path = edges.writeImageAs(str(out / "edge_map.nii.gz"))

sharpened = image_copy(image).sharpen()
sharpened_path = sharpened.writeImageAs(str(out / "sharpened.nii.gz"))

summary_path = write_json(
    out / "filter_summary.json",
    {
        "original_mean": image.getMeanValue(),
        "scaled_mean": scaled.getMeanValue(),
        "smoothed_mean": smoothed.getMeanValue(),
        "edge_max": edges.getMaximumValue(),
        "sharpened_max": sharpened.getMaximumValue(),
    },
)

print_created(scaled_path, smoothed_path, edges_path, sharpened_path, summary_path)