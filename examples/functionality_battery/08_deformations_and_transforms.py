import SimpleITK as sitk

from _common import make_image, make_labelmap, make_roi, output_dir, print_created, write_json
from pyable import deformations


image = make_image()
roi = make_roi()
labelmap = make_labelmap()
out = output_dir("08_deformations_and_transforms")

field, size, origin, spacing, direction = deformations.initialize_deformation_field(image.getImage())
translation = sitk.TranslationTransform(3, (2.0, -1.0, 0.0))
field_from_transform = deformations.transform_to_displacement_field(
    translation,
    size,
    origin,
    spacing,
    direction,
)

warped_image = image.applyTransform(translation)
warped_roi = roi.applyTransformToROI(translation)
warped_labelmap = labelmap.applyTransformToLabelMap(translation)
field_warped_image = image.warpImage(field_from_transform)

paths = [
    warped_image.writeImageAs(str(out / "image_translated_by_transform.nii.gz")),
    warped_roi.writeImageAs(str(out / "roi_translated_nearest_neighbor.nii.gz")),
    warped_labelmap.writeImageAs(str(out / "labelmap_translated_nearest_neighbor.nii.gz")),
    field_warped_image.writeImageAs(str(out / "image_warped_by_displacement_field.nii.gz")),
]

summary_path = write_json(
    out / "deformation_summary.json",
    {
        "field_size": field.GetSize(),
        "field_components": field.GetNumberOfComponentsPerPixel(),
        "field_from_transform_size": field_from_transform.GetSize(),
        "warped_image_size": warped_image.getImageSize(),
        "warped_roi_voxels": warped_roi.getNumberOfNonZeroVoxels(),
        "warped_labelmap_labels": warped_labelmap.getLabels(),
    },
)
paths.append(summary_path)

print_created(*paths)