from pathlib import Path

from pyable import Imaginable, Roiable,LabelMapable


image_path = Path(
    r"C:\Users\montie01\MYDATA\Labrum_Cartilage_Segmentations"
    r"\Seg 1\IM-0004-0001.dcm (1).nii"
)

segmentation_path = Path(
    r"C:\Users\montie01\RESULTS"
    r"\cartilage_refinement_final_all_features"
    r"\Seg_1\IM-0004-0001.dcm__1"
    r"\quality_augmented_joint_full_R.nii.gz"
)

output_directory = Path(__file__).parent / "_outputs" / "paraviewfix_labelmappable"
output_directory.mkdir(parents=True, exist_ok=True)


# Fail immediately with a clear error if either path is wrong.
if not image_path.is_file():
    raise FileNotFoundError(
        f"Image file does not exist:\n{image_path}"
    )

if not segmentation_path.is_file():
    raise FileNotFoundError(
        f"Segmentation file does not exist:\n{segmentation_path}"
    )

print("Image:", image_path)
print("Segmentation:", segmentation_path)


image = Imaginable(filename=str(image_path))
segmentation = LabelMapable(filename=str(segmentation_path))


# Additional safety checks.
if not image.isImageSet():
    raise RuntimeError("The image was not loaded into Imaginable.")

if not segmentation.isImageSet():
    raise RuntimeError("The segmentation was not loaded into LabelMapable.")


# Put the segmentation onto exactly the same image grid.
if not segmentation.isImaginableInTheSameSpace(image):
    print("Resampling segmentation onto the image geometry...")
    segmentation.resampleOnTargetImage(image)


image_output = image.writeParaView(
    str(output_directory / "image.vts"),
    space="lps",
    array_name="intensity",
)


label_files = segmentation.writeParaView(
    output_directory / "labels",
    space="lps",
    prefix="hip",
    label_names={
        1: "acetabular_cartilage",
        2: "femur",
        3: "femoral_cartilage",
        4: "acetabulum",
    },
    smooth_iterations=200,
    relaxation_factor=0.1,
)
print("Created:")
print(image_output)
print(label_files)