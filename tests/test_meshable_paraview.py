"""
Unit tests for pyable.meshable ParaView and VTK/SimpleITK utilities.

Run from the repository root:

    python -m pytest tests/test_meshable_paraview.py -v

The tests use unittest syntax, so they can also be run with:

    python -m unittest tests.test_meshable_paraview -v
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import vtk
from vtk.util import numpy_support

from pyable.meshable import (
    binary_mask_to_polydata,
    get_voxel_to_space_matrix,
    label_map_to_multiblock,
    sitk2vtk,
    sitk_to_structured_grid,
    vtk2sitk,
    write_vtk_dataset,
)


def make_scalar_image() -> tuple[sitk.Image, np.ndarray]:
    """Create a small scalar image with non-default geometry."""
    array = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)

    image = sitk.GetImageFromArray(array)
    image.SetSpacing((1.25, 1.5, 2.0))
    image.SetOrigin((10.0, 20.0, 30.0))
    image.SetDirection(
        (
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        )
    )

    return image, array


def make_vector_image() -> tuple[sitk.Image, np.ndarray]:
    """Create a small three-component vector image."""
    array = np.zeros((2, 3, 4, 3), dtype=np.float32)
    array[..., 0] = 1.0
    array[..., 1] = 2.0
    array[..., 2] = 3.0

    image = sitk.GetImageFromArray(array, isVector=True)
    image.SetSpacing((1.25, 1.5, 2.0))
    image.SetOrigin((10.0, 20.0, 30.0))
    image.SetDirection(
        (
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        )
    )

    return image, array



def get_vtk_dimensions(dataset) -> tuple[int, int, int]:
    """
    Return VTK dataset dimensions across VTK Python binding variants.

    Some VTK builds expose GetDimensions() with no arguments, while others
    require a mutable 3-element list.
    """
    try:
        return tuple(int(v) for v in dataset.GetDimensions())
    except TypeError:
        dimensions = [0, 0, 0]
        dataset.GetDimensions(dimensions)
        return tuple(int(v) for v in dimensions)


class TestMeshableParaView(unittest.TestCase):
    def test_scalar_sitk_vtk_round_trip(self):
        image, expected_array = make_scalar_image()

        vtk_image = sitk2vtk(image)
        recovered = vtk2sitk(vtk_image)

        recovered_array = sitk.GetArrayFromImage(recovered)

        np.testing.assert_array_equal(recovered_array, expected_array)
        np.testing.assert_allclose(recovered.GetSpacing(), image.GetSpacing())
        np.testing.assert_allclose(recovered.GetOrigin(), image.GetOrigin())
        np.testing.assert_allclose(recovered.GetDirection(), image.GetDirection())

    def test_vector_sitk_vtk_round_trip(self):
        image, expected_array = make_vector_image()

        vtk_image = sitk2vtk(image)
        recovered = vtk2sitk(vtk_image)

        recovered_array = sitk.GetArrayFromImage(recovered)

        self.assertEqual(recovered.GetNumberOfComponentsPerPixel(), 3)
        np.testing.assert_array_equal(recovered_array, expected_array)
        np.testing.assert_allclose(recovered.GetSpacing(), image.GetSpacing())
        np.testing.assert_allclose(recovered.GetOrigin(), image.GetOrigin())
        np.testing.assert_allclose(recovered.GetDirection(), image.GetDirection())

    def test_structured_grid_dimensions_and_scalar_order(self):
        image, expected_array = make_scalar_image()

        grid = sitk_to_structured_grid(
            image,
            space="lps",
            array_name="intensity",
        )

        nx, ny, nz = image.GetSize()

        self.assertEqual(get_vtk_dimensions(grid), (nx, ny, nz))
        self.assertEqual(grid.GetNumberOfPoints(), nx * ny * nz)

        vtk_scalars = grid.GetPointData().GetScalars()
        self.assertIsNotNone(vtk_scalars)
        self.assertEqual(vtk_scalars.GetName(), "intensity")

        recovered_values = numpy_support.vtk_to_numpy(vtk_scalars)
        recovered_array = recovered_values.reshape(
            expected_array.shape,
            order="C",
        )

        np.testing.assert_array_equal(recovered_array, expected_array)

    def test_structured_grid_lps_coordinates(self):
        image, _ = make_scalar_image()

        grid = sitk_to_structured_grid(image, space="lps")

        nx, ny, _ = image.GetSize()
        i, j, k = 3, 2, 1

        point_id = i + nx * (j + ny * k)

        expected = image.TransformContinuousIndexToPhysicalPoint(
            (float(i), float(j), float(k))
        )
        actual = grid.GetPoint(point_id)

        np.testing.assert_allclose(actual, expected, atol=1e-5)

    def test_ras_coordinate_conversion(self):
        image, _ = make_scalar_image()

        lps_grid = sitk_to_structured_grid(image, space="lps")
        ras_grid = sitk_to_structured_grid(image, space="ras")

        lps_point = np.asarray(lps_grid.GetPoint(0))
        ras_point = np.asarray(ras_grid.GetPoint(0))

        expected_ras = np.asarray(
            [-lps_point[0], -lps_point[1], lps_point[2]]
        )

        np.testing.assert_allclose(ras_point, expected_ras, atol=1e-5)

    def test_fsl_coordinate_conversion_for_identity_direction(self):
        image, _ = make_scalar_image()
        matrix = get_voxel_to_space_matrix(image, space="fsl")

        nx = image.GetSize()[0]
        sx, sy, sz = image.GetSpacing()

        first = matrix @ np.asarray([0.0, 0.0, 0.0, 1.0])
        last_x = matrix @ np.asarray([float(nx - 1), 0.0, 0.0, 1.0])

        # Identity direction has positive determinant, so FSL reverses X.
        np.testing.assert_allclose(
            first[:3],
            ((nx - 1) * sx, 0.0, 0.0),
            atol=1e-6,
        )
        np.testing.assert_allclose(
            last_x[:3],
            (0.0, 0.0, 0.0),
            atol=1e-6,
        )

        # Check that Y and Z retain their scaled-voxel directions.
        y_point = matrix @ np.asarray([0.0, 1.0, 0.0, 1.0])
        z_point = matrix @ np.asarray([0.0, 0.0, 1.0, 1.0])

        self.assertAlmostEqual(y_point[1], sy)
        self.assertAlmostEqual(z_point[2], sz)

    def test_vector_structured_grid(self):
        image, expected_array = make_vector_image()

        grid = sitk_to_structured_grid(
            image,
            space="lps",
            array_name="displacement",
            data_role="vectors",
        )

        vtk_vectors = grid.GetPointData().GetVectors()

        self.assertIsNotNone(vtk_vectors)
        self.assertEqual(vtk_vectors.GetName(), "displacement")
        self.assertEqual(vtk_vectors.GetNumberOfComponents(), 3)

        recovered = numpy_support.vtk_to_numpy(vtk_vectors)
        expected = expected_array.reshape(-1, 3, order="C")

        np.testing.assert_array_equal(recovered, expected)

    def test_binary_mask_surface(self):
        mask_array = np.zeros((8, 9, 10), dtype=np.uint8)
        mask_array[2:6, 2:7, 3:8] = 1

        mask = sitk.GetImageFromArray(mask_array)
        mask.SetSpacing((1.0, 1.5, 2.0))
        mask.SetOrigin((10.0, 20.0, 30.0))

        surface = binary_mask_to_polydata(
            mask,
            space="lps",
        )

        self.assertIsInstance(surface, vtk.vtkPolyData)
        self.assertGreater(surface.GetNumberOfPoints(), 0)
        self.assertGreater(surface.GetNumberOfCells(), 0)

    def test_empty_mask_surface_raises(self):
        mask = sitk.GetImageFromArray(
            np.zeros((5, 5, 5), dtype=np.uint8)
        )

        with self.assertRaises(ValueError):
            binary_mask_to_polydata(mask)

    def test_label_map_multiblock(self):
        labels_array = np.zeros((10, 10, 10), dtype=np.uint8)
        labels_array[1:4, 1:4, 1:4] = 1
        labels_array[6:9, 6:9, 6:9] = 2

        labels = sitk.GetImageFromArray(labels_array)

        multiblock = label_map_to_multiblock(
            labels,
            space="lps",
        )

        self.assertIsInstance(multiblock, vtk.vtkMultiBlockDataSet)
        self.assertEqual(multiblock.GetNumberOfBlocks(), 2)

        for block_index in range(2):
            block = multiblock.GetBlock(block_index)
            self.assertIsInstance(block, vtk.vtkPolyData)
            self.assertGreater(block.GetNumberOfPoints(), 0)

    def test_write_and_read_vts(self):
        image, _ = make_scalar_image()
        grid = sitk_to_structured_grid(
            image,
            space="lps",
            array_name="intensity",
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "image.vts"

            result = write_vtk_dataset(grid, output_path)

            self.assertEqual(result, str(output_path))
            self.assertTrue(output_path.exists())
            self.assertGreater(output_path.stat().st_size, 0)

            reader = vtk.vtkXMLStructuredGridReader()
            reader.SetFileName(str(output_path))
            reader.Update()

            recovered_grid = reader.GetOutput()

            self.assertEqual(
                get_vtk_dimensions(recovered_grid),
                get_vtk_dimensions(grid),
            )
            self.assertEqual(
                recovered_grid.GetNumberOfPoints(),
                grid.GetNumberOfPoints(),
            )

    def test_write_and_read_vtp(self):
        mask_array = np.zeros((8, 8, 8), dtype=np.uint8)
        mask_array[2:6, 2:6, 2:6] = 1

        mask = sitk.GetImageFromArray(mask_array)
        surface = binary_mask_to_polydata(mask)

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "mask.vtp"

            write_vtk_dataset(surface, output_path)

            self.assertTrue(output_path.exists())
            self.assertGreater(output_path.stat().st_size, 0)

            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(str(output_path))
            reader.Update()

            recovered_surface = reader.GetOutput()

            self.assertGreater(recovered_surface.GetNumberOfPoints(), 0)
            self.assertGreater(recovered_surface.GetNumberOfCells(), 0)

    def test_writer_rejects_wrong_extension_for_dataset(self):
        image, _ = make_scalar_image()
        grid = sitk_to_structured_grid(image)

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "wrong.vtp"

            with self.assertRaises(TypeError):
                write_vtk_dataset(grid, output_path)


if __name__ == "__main__":
    unittest.main()
