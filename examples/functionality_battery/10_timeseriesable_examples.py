from _common import make_timeseries, output_dir, print_created, write_json


timeseries = make_timeseries()
out = output_dir("10_timeseriesable_examples")

frame_0 = timeseries.getFrame(0)
frame_3 = timeseries.getFrame(3)
subset = timeseries.getFrameRange(1, 4)
mean_image = timeseries.getTemporalMean()
variance_image = timeseries.getTemporalVariance()
std_image = timeseries.getTemporalStandardDeviation()
filtered = timeseries.getDuplicate().applyFilterToAllFrames("gaussian", sigma=0.6)

paths = [
    frame_0.writeImageAs(str(out / "frame_00.nii.gz")),
    frame_3.writeImageAs(str(out / "frame_03.nii.gz")),
    subset.writeImageAs(str(out / "frames_01_to_03.nii.gz")),
    mean_image.writeImageAs(str(out / "temporal_mean.nii.gz")),
    variance_image.writeImageAs(str(out / "temporal_variance.nii.gz")),
    std_image.writeImageAs(str(out / "temporal_std.nii.gz")),
    filtered.writeImageAs(str(out / "timeseries_gaussian_filtered.nii.gz")),
]

summary_path = write_json(
    out / "timeseries_summary.json",
    {
        "number_of_frames": timeseries.getNumberOfFrames(),
        "full_size_xyzt": timeseries.getImageSize(),
        "subset_number_of_frames": subset.getNumberOfFrames(),
        "frame_0_mean": frame_0.getMeanValue(),
        "frame_3_mean": frame_3.getMeanValue(),
        "temporal_mean_image_mean": mean_image.getMeanValue(),
    },
)
paths.append(summary_path)

print_created(*paths)