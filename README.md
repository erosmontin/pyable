# pyable
My collection of image and mesh functions
based on SimpleITK

1. imaginable
1. meshable

# Classes
    - Imaginable:
        image data 
    - Roiable:
        Mask with balue 1
    - LabelaMapable
    - LabelMapableROI:
        Roi with multiple values (DEV)
# versions:
- 0.1.1 (Sept, 24)
    - 
- 0.1.0.6 (May, 24)
    - LabelMap are Imaginable with interpolator = sitkNearestNeighbor and dfltuseNearestNeighborExtrapolator=True 

- 0.0.4 pre release
    - updated the concept of change and set
    - dflt interpolation and deflt usenearest..
    - divide and multiply are casted to float and then cast back to their original pixeltype
    - getWavelets
    - left and right functions for HF (tested with Rview)
    - WIP rigid_transform_3D resampleoncanonicalDirections() using [this git repo](https://github.com/nghiaho12/rigid_transform_3D/blob/master/test_rigid_transform_3D.py)
    
[*Dr. Eros Montin, PhD*](http://me.biodimensional.com)
**46&2 just ahead of me!**

