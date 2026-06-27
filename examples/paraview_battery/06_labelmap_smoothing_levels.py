from _common import LABEL_NAMES, align_to_image, load_image, load_labelmap
from _common import output_dir, print_paths


image = load_image()
labelmap = align_to_image(load_labelmap(), image)
out = output_dir("06_labelmap_smoothing_levels")

created = []
for iterations in (0, 20, 50, 100):
    label_outputs = labelmap.writeParaView(
        out / f"smooth_{iterations:03d}",
        space="lps",
        prefix="hip",
        label_names=LABEL_NAMES,
        smooth_iterations=iterations,
        relaxation_factor=0.1,
    )
    created.extend(label_outputs.values())

print_paths(*created)