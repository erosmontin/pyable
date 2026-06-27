from _common import LABEL_NAMES, align_to_image, load_image, load_labelmap
from _common import output_dir, print_paths


image = load_image()
labelmap = align_to_image(load_labelmap(), image)
out = output_dir("05_labelmap_per_label_vtp")

image_output = image.writeParaView(
    out / "image.vts",
    space="lps",
    array_name="intensity",
)

label_outputs = labelmap.writeParaView(
    out / "labels",
    space="lps",
    prefix="hip",
    label_names=LABEL_NAMES,
)

print_paths(image_output, *label_outputs.values())