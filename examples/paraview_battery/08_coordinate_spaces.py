from _common import align_to_image, load_image, load_roi, output_dir, print_paths


image = load_image()
roi = align_to_image(load_roi(), image)
out = output_dir("08_coordinate_spaces")

created = []
for space in ("lps", "ras", "fsl", "index"):
    created.append(
        roi.writeParaView(
            out / f"acetabular_cartilage_{space}.vtp",
            mode="surface",
            space=space,
            smooth_iterations=20,
            relaxation_factor=0.1,
        )
    )

print_paths(*created)