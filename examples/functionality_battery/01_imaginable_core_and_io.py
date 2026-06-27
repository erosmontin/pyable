from _common import make_image, output_dir, print_created, write_json


image = make_image()
out = output_dir("01_imaginable_core_and_io")

copy_path = image.writeImageAs(str(out / "synthetic_image_copy.nii.gz"))

center_index = image.getImageCenterIndex()
center_point = image.getImageCenterCoordinate()
round_trip_index = image.getIndexFromCoordinates(center_point)

summary = {
    "size_xyz": image.getImageSize(),
    "spacing_xyz_mm": image.getImageSpacing(),
    "origin_xyz_mm": image.getImageOrigin(),
    "direction": image.getImageDirection(),
    "pixel_type": image.getImagePixelTypeAsString(),
    "numpy_zyx_shape": image.getImageAsNumpy().shape,
    "numpy_xyz_shape": image.getImageAsNumpyXYZ().shape,
    "center_index_xyz": center_index,
    "center_point_xyz_mm": center_point,
    "round_trip_index_xyz": round_trip_index,
    "min": image.getMinimumValue(),
    "mean": image.getMeanValue(),
    "max": image.getMaximumValue(),
    "std": image.getStdValue(),
    "voxel_volume_mm3": image.getVoxelVolume(),
    "physical_volume_mm3": image.getVolume(),
}
summary_path = write_json(out / "image_summary.json", summary)

print_created(copy_path, summary_path)