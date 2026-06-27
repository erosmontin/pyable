from _common import load_image, output_dir, print_paths


image = load_image()
out = output_dir("01_image_volume_vts")

full_resolution = image.writeParaView(
    out / "image_lps_stride_1.vts",
    space="lps",
    stride=1,
    array_name="intensity",
)

preview_resolution = image.writeParaView(
    out / "image_lps_stride_2.vts",
    space="lps",
    stride=2,
    array_name="intensity",
)

print_paths(full_resolution, preview_resolution)