# pyable
My collection of image and mesh functions
based on SimpleITK

1. imaginable
1. meshable

---

## 🔧 RECENT UPDATES (November 21, 2025)

### ✨ NEW: Deformation & Registration Module (Phase 5)

**Professional-grade support for image warping and registration transforms:**

- **Deformation Module:** 650 lines of utilities for transforms and displacement fields
- **Easy-to-Use API:** 5 new chainable methods on Imaginable, Roiable, LabelMapable
- **Multi-Transform Support:** Rigid, affine, B-spline, displacement fields, composite
- **Label Preservation:** Automatic nearest-neighbor for segmentation warping
- **Comprehensive Docs:** 600+ lines with 15+ working examples
- **100% Test Coverage:** 19 unit tests, all passing

**Quick Example:**
```python
# Warp an image with a registration transform
img = SITKImaginable('moving.nii.gz')
img.applyDisplacementField('deformation.mha', target_image='fixed.nii.gz')
img.write('warped.nii.gz')

# Warp a segmentation (labels automatically preserved!)
roi = Roiable('segmentation.nii.gz')
roi.warpROI('deformation.mha')
roi.write('warped_roi.nii.gz')

# Method chaining
img.applyTransform('transform.tfm').alignGeometry('fixed.nii.gz').cast('uint8')
```

**Documentation:** See [docs/DEFORMATION_WORKFLOW.md](docs/DEFORMATION_WORKFLOW.md)  
**Summary:** See [DEFORMATION_IMPLEMENTATION_SUMMARY.md](DEFORMATION_IMPLEMENTATION_SUMMARY.md)

---

### ✅ Previous Updates: Comprehensive Bug Fix & Testing Initiative (Phase 1-4)

**7 Critical Bugs Fixed:** Logic errors, operator precedence, variable scope, NaN handling  
**31+ Automated Tests Created:** Unit, integration, and regression test suites  
**100% Test Pass Rate:** All validations passing (17/17 regression tests)  
**100% Backward Compatible:** No breaking changes, full API stability maintained

**Documentation:**
- `BUGFIX_REPORT.md` - Detailed analysis of all fixes
- `QUICKREF.md` - Quick reference guide
- `CHANGES.md` - Complete change log with before/after code
- `NAMING_ANALYSIS.md` - Function naming analysis (15 issues identified)

**Status: ✅ Production Ready**

For details, see [BUGFIX_REPORT.md](BUGFIX_REPORT.md)

---

## Installation

To install pyable v3:

```
python3 -m venv "pyable v3"
source "pyable v3"/bin/activate
pip install git+https://github.com/erosmontin/pyable.git@v3-branch
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

