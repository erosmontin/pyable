from _common import make_roi, make_shifted_roi, output_dir, print_created, write_json


reference = make_roi()
test = make_shifted_roi(shift_xyz=(2.0, 0.0, 0.0))
out = output_dir("06_roiable_comparison_and_change_maps")

metrics = test.compareTo(reference)
surface = test.getSurfaceDistances(reference)
change_maps = reference.getChangeMaps(test)

paths = [
    reference.writeImageAs(str(out / "reference_roi.nii.gz")),
    test.writeImageAs(str(out / "shifted_test_roi.nii.gz")),
    change_maps["removed"].writeImageAs(str(out / "removed_voxels.nii.gz")),
    change_maps["added"].writeImageAs(str(out / "added_voxels.nii.gz")),
    change_maps["changed"].writeImageAs(str(out / "changed_voxels.nii.gz")),
]

summary_path = write_json(
    out / "roi_comparison_summary.json",
    {
        "overlap_and_surface_metrics": metrics,
        "surface_distances": surface,
        "compactness": reference.getCompactnessScore(),
        "connected_components": reference.getConnectedComponentCount(),
        "edge_alignment_on_self_image": reference.getEdgeAlignmentScore(reference),
    },
)
paths.append(summary_path)

print_created(*paths)