"""
Pyable - Medical Image Processing Library

A powerful toolkit for working with SimpleITK images, ROIs, and segmentations.
Supports registration transforms, deformations, and multi-label operations.

Classes:
    - Imaginable: Base class for medical images
    - SITKImaginable: SimpleITK image wrapper with chainable methods
    - Roiable: Region of Interest/mask handling
    - LabelMapable: Multi-label segmentation support
    - Fieldable: Vector field operations

Modules:
    - imaginable: Core image classes
    - meshable: Mesh and VTK operations
    - utilizers: ROI comparison metrics
    - deformations: Registration and warping utilities
    - utils: Helper functions

Example:
    >>> from pyable_eros_montin import SITKImaginable, Roiable, LabelMapable
    >>> 
    >>> # Load and manipulate image
    >>> img = SITKImaginable('image.nii.gz')
    >>> img.rotateImage([10, 0, 0]).scaleImage([2, 2, 2]).write('result.nii.gz')
    >>> 
    >>> # Load ROI and apply deformation
    >>> roi = Roiable('segmentation.nii.gz')
    >>> roi.warpROI('deformation.mha').write('warped_roi.nii.gz')
    >>> 
    >>> # Work with multi-label maps
    >>> labels = LabelMapable('anatomy.nii.gz')
    >>> labels.applyTransformToLabelMap('atlas_to_patient.tfm')
    >>> labels.write('patient_anatomy.nii.gz')
"""

from .imaginable import Imaginable, SITKImaginable, Roiable, LabelMapable, Fieldable, LabelMapableROI
from .vectorable import Vectorable, TimeSeriesable
from .plotable import PlotViewer, ScalarPlotter, VectorPlotter, TimeSeriesPlotter, GridPlotter, plotOverlay
from .meshable import vtk2sitk, sitk2vtk
from .utils import processImageDirectory
from . import deformations
from . import utils

# Backward compatibility alias
ROIable = Roiable

__version__ = "3.1.0"
__author__ = "Eros Montin"
__all__ = [
    'Imaginable',
    'SITKImaginable',
    'Roiable',
    'ROIable',  # Backward compatibility
    'LabelMapable',
    'Fieldable',
    'LabelMapableROI',
    'Vectorable',
    'TimeSeriesable',
    'PlotViewer',
    'ScalarPlotter',
    'VectorPlotter',
    'TimeSeriesPlotter',
    'GridPlotter',
    'plotOverlay',
    'processImageDirectory',
    'vtk2sitk',
    'sitk2vtk',
    'deformations',
    'utils',
]
