from pathlib import Path
import json
import sys

import numpy as np
import SimpleITK as sitk


EXAMPLES_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = EXAMPLES_DIR / "_outputs" / "functionality_battery"
DATA_DIR = OUTPUT_ROOT / "synthetic_data"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SPACING = (0.8, 0.8, 1.2)
ORIGIN = (-20.0, -25.0, -18.0)
DIRECTION = (1.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 1.0)


def output_dir(name):
    directory = OUTPUT_ROOT / name
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def set_3d_info(image):
    image.SetSpacing(SPACING)
    image.SetOrigin(ORIGIN)
    image.SetDirection(DIRECTION)
    return image


def synthetic_arrays():
    z, y, x = np.indices((40, 64, 64), dtype=np.float32)

    sphere = ((x - 30.0) ** 2 / 13.0 ** 2
              + (y - 32.0) ** 2 / 15.0 ** 2
              + (z - 20.0) ** 2 / 9.0 ** 2) <= 1.0

    second = ((x - 43.0) ** 2 / 7.0 ** 2
              + (y - 25.0) ** 2 / 10.0 ** 2
              + (z - 24.0) ** 2 / 6.0 ** 2) <= 1.0

    third = ((x - 20.0) ** 2 / 6.0 ** 2
             + (y - 43.0) ** 2 / 8.0 ** 2
             + (z - 15.0) ** 2 / 5.0 ** 2) <= 1.0

    rng = np.random.default_rng(7)
    background = 20.0 + 0.6 * x + 0.3 * y + 0.8 * z
    image = background.astype(np.float32)
    image += 110.0 * sphere.astype(np.float32)
    image += 65.0 * second.astype(np.float32)
    image += 45.0 * third.astype(np.float32)
    image += rng.normal(0.0, 4.0, image.shape).astype(np.float32)

    roi = sphere.astype(np.uint8)

    labelmap = np.zeros_like(roi, dtype=np.uint8)
    labelmap[sphere] = 1
    labelmap[second] = 2
    labelmap[third] = 3

    return image, roi, labelmap


def make_image():
    from pyable import Imaginable

    image_array, _, _ = synthetic_arrays()
    image = sitk.GetImageFromArray(image_array)
    set_3d_info(image)
    return Imaginable(image=image)


def make_roi():
    from pyable import Roiable

    _, roi_array, _ = synthetic_arrays()
    image = sitk.GetImageFromArray(roi_array)
    set_3d_info(image)
    return Roiable(image=image)


def make_labelmap():
    from pyable import LabelMapable

    _, _, label_array = synthetic_arrays()
    image = sitk.GetImageFromArray(label_array)
    set_3d_info(image)
    return LabelMapable(image=image)


def make_shifted_roi(shift_xyz=(2.0, 0.0, 0.0)):
    from pyable import Roiable

    roi = make_roi()
    transform = sitk.TranslationTransform(3, tuple(float(v) for v in shift_xyz))
    shifted = sitk.Resample(
        roi.getImage(),
        roi.getImage(),
        transform,
        sitk.sitkNearestNeighbor,
        0,
        sitk.sitkUInt8,
    )
    return Roiable(image=shifted)


def make_seed_roi():
    from pyable import Roiable

    image_array, _, _ = synthetic_arrays()
    seed = np.zeros_like(image_array, dtype=np.uint8)
    seed[18:22, 29:35, 27:33] = 1
    image = sitk.GetImageFromArray(seed)
    set_3d_info(image)
    return Roiable(image=image)


def make_vector_field():
    from pyable import Vectorable

    z, y, x = np.indices((40, 64, 64), dtype=np.float32)
    vector = np.zeros((40, 64, 64, 3), dtype=np.float32)
    vector[..., 0] = 0.8 * np.sin(x / 12.0)
    vector[..., 1] = 0.6 * np.cos(y / 11.0)
    vector[..., 2] = 0.4 * np.sin(z / 8.0)

    image = sitk.GetImageFromArray(vector, isVector=True)
    set_3d_info(image)
    return Vectorable(image=image)


def make_timeseries():
    from pyable import TimeSeriesable

    base, _, _ = synthetic_arrays()
    frames = []
    for t in range(5):
        frames.append(base + 8.0 * np.sin(t / 4.0 * np.pi) + t * 1.5)
    array_4d = np.stack(frames, axis=0).astype(np.float32)

    image = sitk.GetImageFromArray(array_4d, isVector=False)
    image.SetSpacing((*SPACING, 1.0))
    image.SetOrigin((*ORIGIN, 0.0))
    image.SetDirection(tuple(np.eye(4).reshape(-1)))
    return TimeSeriesable(image=image)


def save_synthetic_dataset():
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    image = make_image()
    roi = make_roi()
    labelmap = make_labelmap()
    vector = make_vector_field()
    timeseries = make_timeseries()

    paths = {
        "image": DATA_DIR / "synthetic_image.nii.gz",
        "roi": DATA_DIR / "synthetic_roi.nii.gz",
        "labelmap": DATA_DIR / "synthetic_labelmap.nii.gz",
        "vector": DATA_DIR / "synthetic_vector_field.mha",
        "timeseries": DATA_DIR / "synthetic_timeseries.nii.gz",
    }

    image.writeImageAs(str(paths["image"]))
    roi.writeImageAs(str(paths["roi"]))
    labelmap.writeImageAs(str(paths["labelmap"]))
    vector.writeImageAs(str(paths["vector"]))
    timeseries.writeImageAs(str(paths["timeseries"]))

    return paths


def roi_copy(roi):
    from pyable import Roiable
    return Roiable(image=sitk.Image(roi.getImage()))


def image_copy(image):
    from pyable import Imaginable
    return Imaginable(image=sitk.Image(image.getImage()))


def labelmap_copy(labelmap):
    from pyable import LabelMapable
    return LabelMapable(image=sitk.Image(labelmap.getImage()))


def json_ready(value):
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(data), indent=2), encoding="utf-8")
    return str(path)


def print_created(*paths):
    print("Created:")
    for path in paths:
        print(path)