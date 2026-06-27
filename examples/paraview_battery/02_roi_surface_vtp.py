from _common import align_to_image, load_image, load_roi, output_dir, print_paths


image = load_image()
roi = align_to_image(load_roi(), image)
out = output_dir("02_roi_surface_vtp")

image_output = image.writeParaView(
    out / "image.vts",
    space="lps",
    array_name="intensity",
)

surface_output = roi.writeParaView(
    out / "acetabular_cartilage_surface.vtp",
    mode="surface",
    space="lps",
)

print_paths(image_output, surface_output)