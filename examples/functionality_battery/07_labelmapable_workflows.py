from _common import labelmap_copy, make_labelmap, output_dir, print_created, write_json


labelmap = make_labelmap()
out = output_dir("07_labelmapable_workflows")

labels = labelmap.getLabels()
roi_1 = labelmap.extractLabel(1)
roi_2 = labelmap.extractLabel(2)
combined = labelmap.combineBinaryMasks([roi_1, roi_2])
priors, classes = labelmap.buildPriors(tau=0.8, blur_sigma_mm=0.4)
union_rois = labelmap.toRoiable()

edited = labelmap_copy(labelmap)
edited.setLabel(3, roi_1.getShell(width_mm=1.6))

paths = [
    roi_1.writeImageAs(str(out / "label_1_as_roi.nii.gz")),
    roi_2.writeImageAs(str(out / "label_2_as_roi.nii.gz")),
    combined.writeImageAs(str(out / "combined_labelmap_from_rois.nii.gz")),
    priors.writeImageAs(str(out / "soft_label_priors.nii.gz")),
    edited.writeImageAs(str(out / "edited_label_3_shell.nii.gz")),
]

summary_path = write_json(
    out / "labelmap_summary.json",
    {
        "labels": labels,
        "centroids_xyz_mm": labelmap.getCentroidCoordinatesPerLabel(),
        "center_of_gravity_xyz_mm": labelmap.getCenterOfGravityCoordinatesPerLabel(),
        "prior_classes": classes,
        "number_of_rois_from_labelmap": len(union_rois),
        "combined_labels": combined.getLabels(),
    },
)
paths.append(summary_path)

print_created(*paths)