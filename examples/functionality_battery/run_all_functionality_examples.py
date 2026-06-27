from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "00_generate_synthetic_dataset.py",
    "01_imaginable_core_and_io.py",
    "02_geometry_resampling_and_transforms.py",
    "03_image_math_filters_and_edges.py",
    "04_segmentation_methods.py",
    "05_roiable_morphology_and_shape.py",
    "06_roiable_comparison_and_change_maps.py",
    "07_labelmapable_workflows.py",
    "08_deformations_and_transforms.py",
    "09_vectorable_examples.py",
    "10_timeseriesable_examples.py",
    "11_plotting_and_overlay_exports.py",
    "12_mesh_and_vtk_exports.py",
    "13_batch_processing_directory.py",
]

for script in SCRIPTS:
    print(f"\n=== {script} ===", flush=True)
    subprocess.run(
        [sys.executable, str(HERE / script)],
        check=True,
        cwd=str(HERE),
    )