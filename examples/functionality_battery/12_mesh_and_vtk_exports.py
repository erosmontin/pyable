from _common import make_image, make_labelmap, make_roi, output_dir, print_created
from pyable import Imaginable, meshable


image = make_image()
roi = make_roi()
labelmap = make_labelmap()
out = output_dir("12_mesh_and_vtk_exports")

vtk_image = meshable.sitk2vtk(image.getImage(), array_name="intensity")
roundtrip_image = meshable.vtk2sitk(vtk_image)

image_volume = image.writeParaView(out / "image_structured_grid.vts", array_name="intensity")
roi_surface = roi.writeParaView(
    out / "roi_smoothed_surface.vtp",
    smooth_iterations=20,
    relaxation_factor=0.1,
)
label_surfaces = labelmap.writeParaView(
    out / "label_surfaces",
    prefix="synthetic",
    smooth_iterations=10,
    relaxation_factor=0.1,
)
label_multiblock = labelmap.writeParaView(
    out / "labelmap_multiblock.vtm",
    smooth_iterations=10,
    relaxation_factor=0.1,
)
roundtrip_output = Imaginable(image=roundtrip_image).writeImageAs(
    str(out / "sitk_vtk_sitk_roundtrip.nii.gz")
)

print_created(image_volume, roi_surface, *label_surfaces.values(), label_multiblock, roundtrip_output)