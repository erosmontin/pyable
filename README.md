# pyable
My collection of image and mesh functions
based on SimpleITK

1. imaginable
1. meshable


## Installation

```
python3 -m venv able
source able/bin/activate
pip install git+https://github.com/erosmontin/pyable.git

```
## Cite Us

1. Montin, E., Belfatto, A., Bologna, M., Meroni, S., Cavatorta, C., Pecori, E., Diletto, B., Massimino, M., Oprandi, M. C., Poggi, G., Arrigoni, F., Peruzzo, D., Pignoli, E., Gandola, L., Cerveri, P., & Mainardi, L. (2020). A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. In Medical &amp; Biological Engineering &amp; Computing (Vol. 58, Issue 4, pp. 843–855). Springer Science and Business Media LLC. https://doi.org/10.1007/s11517-019-02109-4

1. Cavatorta, C., Meroni, S., Montin, E., Oprandi, M. C., Pecori, E., Lecchi, M., Diletto, B., Alessandro, O., Peruzzo, D., Biassoni, V., Schiavello, E., Bologna, M., Massimino, M., Poggi, G., Mainardi, L., Arrigoni, F., Spreafico, F., Verderio, P., Pignoli, E., & Gandola, L. (2021). Retrospective study of late radiation-induced damages after focal radiotherapy for childhood brain tumors. In S. D. Ginsberg (Ed.), PLOS ONE (Vol. 16, Issue 2, p. e0247748). Public Library of Science (PLoS). https://doi.org/10.1371/journal.pone.0247748

## Classes
    - Imaginable:
        image data 
    - Roiable:
        Mask with balue 1
    - LabelaMapable
    - LabelMapableROI:
        Roi with multiple values (DEV)
## versions:
- 0.2.0 (Oct, 24)
    - Fieldable
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

