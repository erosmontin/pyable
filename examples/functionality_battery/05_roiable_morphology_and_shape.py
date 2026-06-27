from _common import make_roi, output_dir, print_created, roi_copy, write_json


roi = make_roi()
out = output_dir("05_roiable_morphology_and_shape")

dilated = roi_copy(roi).dilateMM(radius_mm=1.6)
eroded = roi_copy(roi).erodeMM(radius_mm=1.6)
opened = roi_copy(roi).openMM(radius_mm=1.0)
closed = roi_copy(roi).closeMM(radius_mm=1.0)
shell = roi.getShell(width_mm=2.4)
distance = roi.getDistanceMap()

paths = [
    dilated.writeImageAs(str(out / "roi_dilated_1p6mm.nii.gz")),
    eroded.writeImageAs(str(out / "roi_eroded_1p6mm.nii.gz")),
    opened.writeImageAs(str(out / "roi_opened_1p0mm.nii.gz")),
    closed.writeImageAs(str(out / "roi_closed_1p0mm.nii.gz")),
    shell.writeImageAs(str(out / "roi_shell_2p4mm.nii.gz")),
    distance.writeImageAs(str(out / "roi_distance_map.nii.gz")),
]

summary_path = write_json(
    out / "roi_shape_summary.json",
    {
        "centroid_xyz_mm": roi.getCentroidCoordinates(),
        "center_of_gravity_xyz_mm": roi.getCenterOfGravityCoordinates(),
        "bounding_box_zyx": roi.getBoundingBox(exclude=[0]),
        "original_voxels": roi.getNumberOfNonZeroVoxels(),
        "dilated_voxels": dilated.getNumberOfNonZeroVoxels(),
        "eroded_voxels": eroded.getNumberOfNonZeroVoxels(),
        "morphometrics": roi.getMorphometrics(),
        "principal_extents_mm": roi.getPrincipalExtents(),
        "max_feret_mm": roi.getMaxFeret(),
        "thickness_stats_mm": roi.getThicknessStats(),
    },
)
paths.append(summary_path)

print_created(*paths)