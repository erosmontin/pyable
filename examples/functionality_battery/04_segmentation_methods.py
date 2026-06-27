from _common import make_image, make_seed_roi, output_dir, print_created, write_json


image = make_image()
seed = make_seed_roi()
out = output_dir("04_segmentation_methods")

manual = image.segmentThreshold(lower=95, upper=220)
otsu = image.segmentOtsu(n_bins=128)
multi = image.segmentMultiOtsu(n_thresholds=2, n_bins=128)
connected = image.segmentConnectedThreshold(seed, lower=95, upper=220, n_seeds=50)
watershed = image.segmentMorphologicalWatershed(level=0.08, fully_connected=False)

manual_path = manual.writeImageAs(str(out / "manual_threshold_roi.nii.gz"))
otsu_path = otsu.writeImageAs(str(out / "otsu_roi.nii.gz"))
multi_path = multi.writeImageAs(str(out / "multi_otsu_labelmap.nii.gz"))
connected_path = connected.writeImageAs(str(out / "connected_threshold_roi.nii.gz"))
watershed_path = watershed.writeImageAs(str(out / "morphological_watershed_labelmap.nii.gz"))

summary_path = write_json(
    out / "segmentation_summary.json",
    {
        "manual_voxels": manual.getNumberOfNonZeroVoxels(),
        "otsu_voxels": otsu.getNumberOfNonZeroVoxels(),
        "multi_otsu_labels": multi.getLabels(),
        "connected_voxels": connected.getNumberOfNonZeroVoxels(),
        "watershed_unique_values": sorted(int(v) for v in watershed.getImageUniqueValues(exclude=[])),
    },
)

print_created(manual_path, otsu_path, multi_path, connected_path, watershed_path, summary_path)