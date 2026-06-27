from _common import image_copy, make_image, output_dir, print_created, write_json


image = make_image()
out = output_dir("02_geometry_resampling_and_transforms")

resized = image_copy(image).changeImageSize([48, 48, 30])
resized_path = resized.writeImageAs(str(out / "resized_48x48x30.nii.gz"))

resampled = image_copy(image).changeImageSpacing([1.2, 1.2, 1.6])
resampled_path = resampled.writeImageAs(str(out / "resampled_spacing_1p2_1p2_1p6.nii.gz"))

padded = image_copy(image).padImage([4, 4, 2], [4, 4, 2], padding_value=0)
padded_path = padded.writeImageAs(str(out / "padded.nii.gz"))

cropped = image_copy(image).cropImage([8, 8, 4], [56, 56, 36])
cropped_path = cropped.writeImageAs(str(out / "cropped.nii.gz"))

translated = image_copy(image).translateImage([2.0, -1.0, 0.0])
translated_path = translated.writeImageAs(str(out / "translated.nii.gz"))

summary_path = write_json(
    out / "geometry_summary.json",
    {
        "original_size": image.getImageSize(),
        "resized_size": resized.getImageSize(),
        "resampled_spacing": resampled.getImageSpacing(),
        "padded_size": padded.getImageSize(),
        "cropped_size": cropped.getImageSize(),
        "translated_origin": translated.getImageOrigin(),
    },
)

print_created(resized_path, resampled_path, padded_path, cropped_path, translated_path, summary_path)