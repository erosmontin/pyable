from _common import align_to_image, load_image, load_labelmap, output_dir, print_paths


image = load_image()
labelmap = align_to_image(load_labelmap(), image)
out = output_dir("07_labelmap_multiblock_vtm")

multiblock_output = labelmap.writeParaView(
    out / "hip_labels_smoothed.vtm",
    space="lps",
    smooth_iterations=20,
    relaxation_factor=0.1,
)

print_paths(multiblock_output)