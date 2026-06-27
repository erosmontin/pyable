from _common import align_to_image, load_image, load_roi, output_dir, print_paths


image = load_image()
roi = align_to_image(load_roi(), image)
out = output_dir("03_roi_surface_smoothing_levels")

created = []
for iterations in (0, 10, 20, 50):
    created.append(
        roi.writeParaView(
            out / f"acetabular_cartilage_smooth_{iterations:03d}.vtp",
            mode="surface",
            space="lps",
            smooth_iterations=iterations,
            relaxation_factor=0.1,
        )
    )

print_paths(*created)