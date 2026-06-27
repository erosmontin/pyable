from pathlib import Path
import sys


EXAMPLES_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[2]

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

IMAGE_PATH = Path(
    r"C:\Users\montie01\MYDATA\Labrum_Cartilage_Segmentations"
    r"\Seg 1\IM-0004-0001.dcm (1).nii"
)

ROI_PATH = Path(
    r"C:\Users\montie01\RESULTS"
    r"\cartilage_refinement_final_all_features"
    r"\Seg_1\IM-0004-0001.dcm__1"
    r"\quality_augmented_acetabular_cartilage_R.nii.gz"
)

LABELMAP_PATH = Path(
    r"C:\Users\montie01\RESULTS"
    r"\cartilage_refinement_final_all_features"
    r"\Seg_1\IM-0004-0001.dcm__1"
    r"\quality_augmented_joint_full_R.nii.gz"
)

OUTPUT_ROOT = EXAMPLES_DIR / "_outputs" / "paraview_battery"

LABEL_NAMES = {
    1: "acetabular_cartilage",
    2: "femur",
    3: "femoral_cartilage",
    4: "acetabulum",
}


def require_file(path, label):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"{label} file does not exist:\n{path}")
    return path


def ensure_inputs():
    require_file(IMAGE_PATH, "Image")
    require_file(ROI_PATH, "ROI")
    require_file(LABELMAP_PATH, "Label map")


def output_dir(name):
    directory = OUTPUT_ROOT / name
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def load_image():
    ensure_inputs()
    from pyable import Imaginable

    image = Imaginable(filename=str(IMAGE_PATH))
    if not image.isImageSet():
        raise RuntimeError("The image was not loaded into Imaginable.")
    return image


def load_roi():
    ensure_inputs()
    from pyable import Roiable

    roi = Roiable(filename=str(ROI_PATH))
    if not roi.isImageSet():
        raise RuntimeError("The ROI was not loaded into Roiable.")
    return roi


def load_labelmap():
    ensure_inputs()
    from pyable import LabelMapable

    labelmap = LabelMapable(filename=str(LABELMAP_PATH))
    if not labelmap.isImageSet():
        raise RuntimeError("The label map was not loaded into LabelMapable.")
    return labelmap


def align_to_image(segmentation, image):
    if not segmentation.isImaginableInTheSameSpace(image):
        print("Resampling segmentation onto the image geometry...")
        segmentation.resampleOnTargetImage(image)
    return segmentation


def print_paths(*paths):
    print("Created:")
    for path in paths:
        print(path)