from _common import align_to_image, load_image, load_roi, output_dir, print_paths


image = load_image()
roi = align_to_image(load_roi(), image)
out = output_dir("04_roi_volume_vts")

volume_output = roi.writeParaView(
    out / "acetabular_cartilage_mask_volume.vts",
    mode="volume",
    space="lps",
    array_name="mask",
)

print_paths(volume_output)