from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "00_check_inputs.py",
    "01_image_volume_vts.py",
    "02_roi_surface_vtp.py",
    "03_roi_surface_smoothing_levels.py",
    "04_roi_volume_vts.py",
    "05_labelmap_per_label_vtp.py",
    "06_labelmap_smoothing_levels.py",
    "07_labelmap_multiblock_vtm.py",
    "08_coordinate_spaces.py",
]

for script in SCRIPTS:
    print(f"\n=== {script} ===", flush=True)
    subprocess.run(
        [sys.executable, str(HERE / script)],
        check=True,
        cwd=str(HERE),
    )