#!/usr/bin/env python
"""
Example: load the specified ROI and render it in VTK using `renderIsosurface()`.

Usage:
    python example_show_roi_vtk.py

The script defaults to interactive display (VTK window). If you want a non-interactive
run (for CI or headless), pass `--no-show` and the script will create the isosurface
without opening a window and print summary statistics.
"""
import os
import sys
import argparse

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

from pyable.imaginable import Roiable, SITKImaginable
import SimpleITK as sitk

ROI_PATH = "/media/erosm/DATA/aging20/DATA/HCP/HCA6002236_V1_MR/aparc+aseg.nii.gz"


def main():
    parser = argparse.ArgumentParser(description="Load ROI and render with VTK")
    parser.add_argument("--no-show", action="store_true", help="Do not open interactive VTK window; run non-interactive and print stats")
    parser.add_argument("--path", type=str, default=ROI_PATH, help="Path to ROI NIfTI file")
    parser.add_argument("--color", type=float, nargs=3, default=(0.0, 1.0, 0.0), help="RGB color (0-1) for surface")
    parser.add_argument("--opacity", type=float, default=0.9, help="Opacity 0-1 for surface")
    args = parser.parse_args()

    path = args.path

    if not os.path.exists(path):
        print(f"ERROR: ROI file not found: {path}")
        print("Please verify the path or copy the file to the expected location.")
        sys.exit(2)

    # Prefer to create Roiable by filename so metadata is read automatically
    try:
        roi = Roiable(filename=path)
    except Exception:
        # Fallback: read with SimpleITK then wrap
        sitk_img = sitk.ReadImage(path)
        roi = Roiable(image=sitk_img)

    print(f"Loaded ROI: {path}")
    print(f"  Pixel type: {roi.getImagePixelTypeAsString()[1]}")
    print(f"  Size: {roi.getImageSize()}")
    print(f"  Non-zero voxels: {roi.getNumberOfNonZeroVoxels()}")

    # Non-interactive: create actor and print polygon/point counts
    if args.no_show:
        actor, renderer, window = roi.renderIsosurface(show=False, color=tuple(args.color), opacity=args.opacity)
        polydata = actor.GetMapper().GetInput()
        print("Isosurface created (non-interactive):")
        print(f"  Polygons: {polydata.GetNumberOfCells()}")
        print(f"  Points: {polydata.GetNumberOfPoints()}")
        sys.exit(0)

    # Interactive: open VTK window (blocks until closed)
    print("Opening interactive VTK window. Close the window to exit.")
    try:
        roi.renderIsosurface(show=True, color=tuple(args.color), opacity=args.opacity)
    except Exception as e:
        print(f"Failed to open interactive VTK window: {e}")
        print("You can run with --no-show to generate geometry without display.")
        sys.exit(1)


if __name__ == '__main__':
    main()
