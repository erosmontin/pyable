import base64
import sys
import types
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import importlib.util
import matplotlib.pyplot as plt
import numpy as np
import pytest

# Ensure the workspace copy of pyable is imported before any installed package.
REPO_ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(REPO_ROOT))

sys.modules.setdefault("imaginable", types.SimpleNamespace())
spec = importlib.util.spec_from_file_location(
    "pyable_utils_under_test",
    REPO_ROOT / "pyable" / "utils.py",
)
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)


def test_make_axis_index_pairs_treats_multi_axis_index_as_point_by_default():
    pairs, multi_slice = utils.makeAxisIndexPairs([0, 1, 2], [10, 20, 30])

    assert multi_slice is True
    assert pairs == [
        (0, 10),
        (1, 20),
        (2, 30),
    ]


def test_make_axis_index_pairs_can_use_axis_major_cartesian_product():
    pairs, multi_slice = utils.makeAxisIndexPairs(
        [0, 1, 2],
        [10, 20, 30],
        index_mode="cartesian",
    )

    assert multi_slice is True
    assert pairs == [
        (0, 10), (0, 20), (0, 30),
        (1, 10), (1, 20), (1, 30),
        (2, 10), (2, 20), (2, 30),
    ]


def test_make_axis_index_pairs_supports_multi_axis_point_with_slice_offsets():
    pairs, multi_slice = utils.makeAxisIndexPairs(
        [0, 2],
        [10, 20, 30],
        slice_offsets=[-1, 0, 1],
    )

    assert multi_slice is True
    assert pairs == [
        (0, 9), (0, 10), (0, 11),
        (2, 29), (2, 30), (2, 31),
    ]


def test_make_axis_index_pairs_supports_multi_axis_scalar_with_slice_offsets():
    pairs, multi_slice = utils.makeAxisIndexPairs([0, 2], 20, slice_offsets=[-1, 0, 1])

    assert multi_slice is True
    assert pairs == [
        (0, 19), (0, 20), (0, 21),
        (2, 19), (2, 20), (2, 21),
    ]


def test_make_axis_index_pairs_rejects_slice_offsets_with_cartesian_multi_index():
    with pytest.raises(ValueError, match="slice_offsets"):
        utils.makeAxisIndexPairs(
            0,
            [20, 21],
            slice_offsets=[-1, 0, 1],
            index_mode="cartesian",
        )


def test_overlay_numpy_image_and_labelmap_adds_colorbar_to_current_axes():
    plt.close("all")
    image = np.arange(100, dtype=float).reshape(10, 10)
    labelmap = np.zeros((10, 10), dtype=int)
    labelmap[2:5, 3:7] = 1

    utils.overlayNumpyImageAndNumpyLabelmap(
        image,
        labelmap,
        show=False,
        labelmap_name="Segment",
    )

    fig = plt.gcf()
    assert len(fig.axes) == 2
    assert fig.axes[1].get_ylabel() == "Segment"

    plt.close(fig)


def test_overlay_numpy_image_and_labelmap_uses_supplied_axes():
    plt.close("all")
    fig, ax = plt.subplots()
    image = np.arange(25, dtype=float).reshape(5, 5)
    labelmap = np.zeros((5, 5), dtype=int)
    labelmap[1:4, 1:4] = 2

    result = utils.overlayNumpyImageAndNumpyLabelmap(
        image,
        labelmap,
        ax=ax,
        show=False,
        colorbar=True,
        labelmap_name="ROI",
    )

    assert result["axis"] is ax
    assert result["figure"] is fig
    assert len(fig.axes) == 2
    assert fig.axes[1].get_ylabel() == "ROI"

    plt.close(fig)


def test_overlay_numpy_image_and_labelmap_grid_uses_square_layout_and_titles():
    plt.close("all")
    images = [np.full((5, 5), i, dtype=float) for i in range(5)]
    labelmaps = []
    for i in range(5):
        labelmap = np.zeros((5, 5), dtype=int)
        labelmap[1:4, 1:4] = i + 1
        labelmaps.append(labelmap)
    titles = [f"Slice {i}" for i in range(5)]

    result = utils.overlayNumpyImageAndNumpyLabelmapGrid(
        images,
        labelmaps,
        titles=titles,
        show=False,
    )

    axes = result["axes"]
    axes_flat = axes.ravel()
    assert axes.shape == (3, 3)
    assert [ax.get_title() for ax in axes_flat[:5]] == titles
    assert all(len(ax.images) == 2 for ax in axes_flat[:5])
    assert all(len(ax.images) == 0 for ax in axes_flat[5:])

    plt.close(result["figure"])


def test_overlay_numpy_image_and_labelmap_to_image_returns_rgba_array():
    image = np.zeros((5, 5), dtype=float)
    labelmap = np.zeros((5, 5), dtype=int)
    labelmap[2, 2] = 1

    rgba = utils.overlayNumpyImageAndNumpyLabelmapToImage(
        image,
        labelmap,
        alpha_value=1.0,
        image_vmin=0,
        image_vmax=1,
        labelmap_vmin=0,
        labelmap_vmax=1,
        origin="upper",
    )

    assert rgba.shape == (5, 5, 4)
    assert rgba.dtype == np.uint8
    assert np.all(rgba[..., 3] == 255)
    assert not np.array_equal(rgba[2, 2], rgba[0, 0])


def test_overlay_numpy_image_and_labelmap_to_image_can_return_base64_png():
    image = np.zeros((5, 5), dtype=float)
    labelmap = np.zeros((5, 5), dtype=int)
    labelmap[1:4, 1:4] = 1

    encoded = utils.overlayNumpyImageAndNumpyLabelmapToImage(
        image,
        labelmap,
        as_base64=True,
    )

    assert base64.b64decode(encoded).startswith(b"\x89PNG\r\n\x1a\n")

    data_uri = utils.overlayNumpyImageAndNumpyLabelmapToImage(
        image,
        labelmap,
        as_base64=True,
        data_uri=True,
    )

    assert data_uri.startswith("data:image/png;base64,")


def test_overlay_numpy_image_and_labelmap_to_image_can_add_title_and_save():
    image = np.zeros((5, 7), dtype=float)
    labelmap = np.zeros((5, 7), dtype=int)
    labelmap[2, 3] = 1
    output = Path(__file__).with_name("overlay_test_output.png")

    try:
        rgba = utils.overlayNumpyImageAndNumpyLabelmapToImage(
            image,
            labelmap,
            title="Center slice",
            save=output,
            origin="upper",
        )

        assert rgba.shape[0] > image.shape[0]
        assert rgba.shape[1] >= image.shape[1]
        assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    finally:
        output.unlink(missing_ok=True)


def test_overlay_numpy_image_and_labelmap_grid_to_image_makes_tight_montage():
    images = [np.full((5, 6), i, dtype=float) for i in range(5)]
    labelmaps = []
    for i in range(5):
        labelmap = np.zeros((5, 6), dtype=int)
        labelmap[1:4, 2:5] = i + 1
        labelmaps.append(labelmap)

    rgba = utils.overlayNumpyImageAndNumpyLabelmapGridToImage(
        images,
        labelmaps,
        ncols=3,
        tile_gap=0,
        origin="upper",
    )

    assert rgba.shape == (10, 18, 4)
    assert rgba.dtype == np.uint8

    encoded = utils.overlayNumpyImageAndNumpyLabelmapGridToImage(
        images,
        labelmaps,
        ncols=3,
        title="2.5D context",
        titles=[f"S{i}" for i in range(5)],
        as_base64=True,
    )

    assert base64.b64decode(encoded).startswith(b"\x89PNG\r\n\x1a\n")
