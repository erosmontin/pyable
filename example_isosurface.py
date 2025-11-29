#!/usr/bin/env python
"""
Quick start examples for VTK isosurface rendering.
Run this to see the basic functionality in action.
"""

import sys
import numpy as np
import SimpleITK as sitk

sys.path.insert(0, '/home/erosm/pyable')

from pyable.imaginable import SITKImaginable, Roiable


def example_1_continuous_image():
    """Example 1: Render a continuous image at mean intensity."""
    print("\n" + "="*60)
    print("EXAMPLE 1: Continuous Image Rendering")
    print("="*60)
    
    # Create synthetic 3D Gaussian volume
    size = (100, 100, 100)
    image = sitk.Image(size, sitk.sitkFloat32)
    center = [s // 2 for s in size]
    radius = 30
    
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                dist = np.sqrt((i - center[0])**2 + (j - center[1])**2 + (k - center[2])**2)
                value = np.exp(-(dist**2) / (2 * radius**2)) * 255
                image.SetPixel([i, j, k], value)
    
    # Wrap in Imaginable
    img = SITKImaginable(image=image)
    
    print("Created synthetic Gaussian volume (100×100×100)")
    print(f"Intensity range: [{img.getMinimumValue():.1f}, {img.getMaximumValue():.1f}]")
    print(f"Mean intensity: {img.getMeanValue():.1f}")
    print("\nUsage:")
    print("  img.renderIsosurface()                    # Auto (mean intensity)")
    print("  img.renderIsosurface(isosurface_value=100)  # Custom value")
    print("  img.renderIsosurface(color=(0,1,0))       # Green color")
    
    return img


def example_2_roi():
    """Example 2: Render an ROI boundary."""
    print("\n" + "="*60)
    print("EXAMPLE 2: ROI/Segmentation Boundary")
    print("="*60)
    
    # Create synthetic binary ROI
    size = (100, 100, 100)
    image = sitk.Image(size, sitk.sitkUInt8)
    center = [s // 2 for s in size]
    radius = 25
    
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                dist = np.sqrt((i - center[0])**2 + (j - center[1])**2 + (k - center[2])**2)
                if dist <= radius:
                    image.SetPixel([i, j, k], 1)
                else:
                    image.SetPixel([i, j, k], 0)
    
    # Wrap in Roiable
    roi = Roiable(image=image)
    
    print("Created synthetic binary ROI (100×100×100)")
    print(f"Non-zero voxels: {roi.getNumberOfNonZeroVoxels()}")
    print("\nUsage:")
    print("  roi.renderIsosurface()                    # Auto (boundary at 0.5)")
    print("  roi.renderIsosurface(color=(0,1,0))       # Green boundary")
    print("  roi.renderIsosurface(opacity=0.7)         # Semi-transparent")
    
    return roi


def example_3_batch_processing():
    """Example 3: Generate isosurfaces without displaying."""
    print("\n" + "="*60)
    print("EXAMPLE 3: Batch Processing (No Display)")
    print("="*60)
    
    # Create synthetic volume
    size = (100, 100, 100)
    image = sitk.Image(size, sitk.sitkFloat32)
    center = [s // 2 for s in size]
    radius = 30
    
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                dist = np.sqrt((i - center[0])**2 + (j - center[1])**2 + (k - center[2])**2)
                value = np.exp(-(dist**2) / (2 * radius**2)) * 255
                image.SetPixel([i, j, k], value)
    
    img = SITKImaginable(image=image)
    
    print("Creating isosurfaces without displaying...")
    print("\nUsage:")
    print("  actor, _, _ = img.renderIsosurface(show=False)")
    print("  # Returns actor, renderer, window tuple")
    print("\nExample code:")
    print("  actor, _, _ = img.renderIsosurface(isosurface_value=100, show=False)")
    print("  polydata = actor.GetMapper().GetInput()")
    print("  print(f'Polygons: {polydata.GetNumberOfCells()}')")
    print("  print(f'Points: {polydata.GetNumberOfPoints()}')")
    
    # Actually generate one
    actor, _, _ = img.renderIsosurface(isosurface_value=100, show=False)
    polydata = actor.GetMapper().GetInput()
    
    print("\nGenerated isosurface statistics:")
    print(f"  Polygons: {polydata.GetNumberOfCells()}")
    print(f"  Points: {polydata.GetNumberOfPoints()}")
    
    return actor


def example_4_colors():
    """Example 4: Different colors for isosurfaces."""
    print("\n" + "="*60)
    print("EXAMPLE 4: Color Options")
    print("="*60)
    
    size = (100, 100, 100)
    image = sitk.Image(size, sitk.sitkFloat32)
    center = [s // 2 for s in size]
    
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                dist = np.sqrt((i - center[0])**2 + (j - center[1])**2 + (k - center[2])**2)
                value = np.exp(-(dist**2) / (2 * 30**2)) * 255
                image.SetPixel([i, j, k], value)
    
    img = SITKImaginable(image=image)
    
    print("Common colors (RGB, 0-1 range):")
    colors = [
        ((1.0, 0.0, 0.0), "Red"),
        ((0.0, 1.0, 0.0), "Green"),
        ((0.0, 0.0, 1.0), "Blue"),
        ((1.0, 1.0, 0.0), "Yellow"),
        ((1.0, 0.0, 1.0), "Magenta"),
        ((0.0, 1.0, 1.0), "Cyan"),
    ]
    
    print("\nUsage examples:")
    for color, name in colors:
        print(f"  img.renderIsosurface(color={color})  # {name}")
    
    # Generate all colors
    print("\nGenerating all color variations...")
    for color, name in colors:
        actor, _, _ = img.renderIsosurface(isosurface_value=100, show=False, color=color)
        print(f"  ✓ {name}")
    
    return img


def example_5_transparency():
    """Example 5: Transparency/opacity control."""
    print("\n" + "="*60)
    print("EXAMPLE 5: Transparency Control")
    print("="*60)
    
    size = (100, 100, 100)
    image = sitk.Image(size, sitk.sitkFloat32)
    center = [s // 2 for s in size]
    
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                dist = np.sqrt((i - center[0])**2 + (j - center[1])**2 + (k - center[2])**2)
                value = np.exp(-(dist**2) / (2 * 30**2)) * 255
                image.SetPixel([i, j, k], value)
    
    img = SITKImaginable(image=image)
    
    print("Transparency examples:")
    opacity_values = [0.2, 0.5, 0.8, 1.0]
    
    print("\nUsage examples:")
    for opacity in opacity_values:
        print(f"  img.renderIsosurface(opacity={opacity})  # {int(opacity*100)}% opaque")
    
    print("\nGenerating isosurfaces with different opacities...")
    for opacity in opacity_values:
        actor, _, _ = img.renderIsosurface(isosurface_value=100, show=False, opacity=opacity)
        print(f"  ✓ opacity={opacity}")
    
    return img


def main():
    """Run all examples."""
    print("\n" + "="*70)
    print("ISOSURFACE RENDERING - QUICK START EXAMPLES")
    print("="*70)
    
    # Create examples
    img = example_1_continuous_image()
    roi = example_2_roi()
    actor = example_3_batch_processing()
    img2 = example_4_colors()
    img3 = example_5_transparency()
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    print("""
✓ All examples executed successfully!

Key Features:
  • Continuous images: render at any intensity value
  • ROIs: automatically renders boundary at 0.5
  • Batch processing: generate without displaying (show=False)
  • Colors: full RGB color support (0-1 range)
  • Transparency: control opacity for 3D visualization
  • Multi-component: extract components before rendering
  • 4D images: extract time frame before rendering

To display interactively, call:
  img.renderIsosurface()

To render without display (for automation):
  actor, _, _ = img.renderIsosurface(show=False)

For more details, see: docs/ISOSURFACE_RENDERING_GUIDE.md
    """)
    
    print("="*70)


if __name__ == "__main__":
    main()
