# FUNCTION NAMING ANALYSIS - PYABLE PACKAGE

**Date:** November 21, 2025  
**Status:** Detailed Code Review

---

## EXECUTIVE SUMMARY

Your function naming has **both strengths and issues**:

✅ **Good Practices:**
- Consistent Get/Set patterns for properties
- Clear descriptive names for most operations
- Good use of camelCase

⚠️ **Issues Found:** 13 naming problems
- Inconsistent naming conventions
- Typos in function/class names
- Unclear abbreviations
- Missing method aliases
- Grammatical errors

---

## CRITICAL ISSUES

### 1. **Duplicate "getCenterOfGravity" Bug** (Already Fixed) ✅
**Location:** `imaginable.py` Line 1411  
**Status:** FIXED in Phase 1

```python
# WAS: def getCenterOfGravitygetCenterOfGravityIndex(self):  ❌
# NOW: def getCenterOfGravityIndex(self):                     ✅
```

---

## ISSUES BY CATEGORY

### A. TYPOS & MISSPELLINGS (3 issues)

#### 1. **`getMaskedNunmpyArray`** - Typo in function name
**File:** `imaginable.py` Line 105  
**Issue:** "Nunmpy" should be "Numpy"

```python
# CURRENT (WRONG):
def getMaskedNunmpyArray(IM,ROI):
    
# SHOULD BE:
def getMaskedNumpyArray(IM,ROI):
```

**Usage:** Called in `imaginable.py` line 315  
**Fix Priority:** 🔴 HIGH (used internally)

---

#### 2. **`copythethreeinfosonandsetthemtimage`** - Multiple typos
**File:** `imaginable.py` Line 260  
**Issue:** "threeinfoson" should be "threeinfo" + missing caps

```python
# CURRENT (WRONG):
def copythethreeinfosonandsetthemtimage(source,reference):
    sp,o,d=getSITKImageInfo(reference)
    return setSITKImageInfo(source,sp,o,d)

# SHOULD BE:
def copyTheThreeInfosAndSetThemToImage(source,reference):
    sp,o,d=getSITKImageInfo(reference)
    return setSITKImageInfo(source,sp,o,d)
```

**Better Alternative:** `copyImageMetadata()` or `copyImageInfo()`

**Fix Priority:** 🔴 HIGH (exported function)

---

#### 3. **`setSITKImageInforFromImage`** - Typo
**File:** `imaginable.py` Line 288  
**Issue:** "Infor" should be "Info"

```python
# CURRENT (WRONG):
def setSITKImageInforFromImage(nda,ima):

# SHOULD BE:
def setSITKImageInfoFromImage(nda,ima):
```

**Fix Priority:** 🟡 MEDIUM (appears to be unused)

---

#### 4. **README.md - Typo in class name**
**File:** `README.md` Line 30  
**Issue:** "LabelaMapable" should be "LabelMapable"

```markdown
# CURRENT (WRONG):
- LabelaMapable

# SHOULD BE:
- LabelMapable
```

Also: "balue" should be "value"

```markdown
# CURRENT (WRONG):
- Roiable: Mask with balue 1

# SHOULD BE:
- Roiable: Mask with value 1
```

**Fix Priority:** 🟡 MEDIUM (documentation)

---

### B. INCONSISTENT NAMING CONVENTIONS (5 issues)

#### 5. **Inconsistent Method Pair Names**
**File:** `imaginable.py`  
**Issue:** Some pairs follow `getX` / `setX`, but some don't

```python
# GOOD PAIRS:
getImage() / setImage()
getImageSpacing() / setImageSpacing()
getImageOrigin() / setImageOrigin()
getImageDirection() / setImageDirection()

# INCONSISTENT:
getImageSize() / changeImageSize()  # Should be setImageSize() or resize()
getImageSpacing() / changeImageSpacing()  # Should be setImageSpacing()
getImageOrigin() / changeImageOrigin()  # Should be setImageOrigin()
getImageDirection() / changeImageDirection()  # Should be setImageDirection()
```

**Problem:** Mixed terminology: `set` vs `change`  
**Recommendation:**
- Use `set*` for direct property assignment
- Use `*AndResample()` or `*WithInterpolation()` for operations with side effects

**Fix Priority:** 🟠 MEDIUM-HIGH

---

#### 6. **Arithmetic Operations - Inconsistent Naming**
**File:** `imaginable.py` Lines 932-952  
**Issue:** Some use `add()`, some use `*Image()` suffix

```python
# CURRENT (INCONSISTENT):
add(self, toadd)              # Generic name
multiply(self, toadd)         # Generic name
subtract(self, toadd)         # Generic name
divide(self, toadd)           # Generic name

# SHOULD BE (one of these approaches):
addImage(self, image)         # Consistent suffix
multiplyImage(self, image)
subtractImage(self, image)
divideImage(self, image)

# OR add Image Variants (aliases):
addImage() / add()
multiplyImage() / multiply()  # With appropriate parameter naming
```

**Current Parameter Names:** `toadd` is unclear (could be image or scalar)

**Fix Priority:** 🟡 MEDIUM

---

#### 7. **Get Bounding Box - Confusing Name**
**File:** `imaginable.py` Line 1019  
**Issue:** Parameter `exclude` is confusing

```python
# CURRENT (UNCLEAR):
def getBoundingBox(self, exclude=[0]):
    """Exclude values... or exclude axes?"""

# SHOULD BE:
def getBoundingBox(self, excludeValues=None):
    """Clearer intent"""
    if excludeValues is None:
        excludeValues = [0]

# OR BETTER:
def getBoundingBox(self, background_value=0):
    """More explicit"""
```

**Fix Priority:** 🟡 MEDIUM

---

#### 8. **Private Methods - Inconsistent Naming**
**File:** `imaginable.py`  
**Issue:** Mix of `__name__()` and `_name()` styles

```python
# CURRENT (INCONSISTENT):
__readImage__()              # Double underscore
__filterSelfAndImage__()     # Double underscore
__filterSelfAndImageMat__()  # Double underscore
__getOrientation__()         # Double underscore
__derodeRadius__()           # Double underscore
__tellme__()                 # Double underscore (non-standard)
__RegionExtractor__()        # Double underscore + CamelCase

# STANDARD PYTHON:
_readImage()                 # Single underscore for "protected"
_filterSelfAndImage()
# Double underscore (__) only for name-mangling (very rare)
```

**Problem:** Double underscores are meant for name mangling, not "private" marking  
**Fix Priority:** 🟠 MEDIUM-HIGH (affects API clarity)

---

#### 9. **Filter Methods - Poor Naming**
**File:** `imaginable.py` Line 945  
**Issue:** `filterValues()` has unclear semantics

```python
# CURRENT (UNCLEAR):
def filterValues(self, values):
    """What does this do? Remove? Keep? """

# SHOULD BE (one of):
def keepValues(self, values)        # More explicit
def retainValues(self, values)      # More explicit
def filterByValues(self, values)    # Or add documentation
```

**Fix Priority:** 🟡 MEDIUM

---

### C. UNCLEAR ABBREVIATIONS (2 issues)

#### 10. **`dflt` Prefix - Unclear Abbreviation**
**File:** `imaginable.py` Multiple lines  
**Issue:** `dflt` is not immediately clear

```python
# CURRENT:
self.dfltInterpolator = sitk.sitkLinear
self.dfltuseNearestNeighborExtrapolator = False

# SHOULD BE:
self.defaultInterpolator = sitk.sitkLinear
self.defaultUseNearestNeighborExtrapolator = False

# OR (if conciseness important):
self.default_interpolator = sitk.sitkLinear
self.default_use_nearest_neighbor_extrapolator = False
```

**Impact:** Affects readability for new users  
**Fix Priority:** 🟡 MEDIUM

---

#### 11. **`UD` / `LR` - Cryptic Variable Names**
**File:** `imaginable.py` Line 1175  
**Issue:** What do UD/LR mean? (Presumably Up-Down / Left-Right?)

```python
# CURRENT (UNCLEAR):
# if self.UD:
#     o=np.flipud(o)
# if self.LR:
#     o=np.fliplr(o)

# SHOULD BE:
# if self.flip_up_down:
#     o = np.flipud(o)
# if self.flip_left_right:
#     o = np.fliplr(o)
```

**Fix Priority:** 🟡 MEDIUM (commented out, low priority)

---

### D. INCONSISTENT NAMING SCHEMES (2 issues)

#### 12. **"Imaginable" vs "SITK" Prefix Inconsistency**
**File:** `imaginable.py`  
**Issue:** Some functions use "Imaginable" prefix, some don't

```python
# INCONSISTENT:
numpyToImaginable()        # Noun → Imaginable
getImaginableSlice()       # Imaginable → Noun
overlayAble()              # Just "Able" (incomplete name!)
instantiateAnotherAble()   # "Able" instead of "Imaginable"

# SHOULD BE (one approach):
numpyToImaginable()        # Keep
imaginableFromNumpy()      # More consistent
getImaginableSlice()       # Keep
createFromNumpy()          # Alternative
```

**Fix Priority:** 🟡 MEDIUM

---

#### 13. **`viewSagittal()` vs `getSliceNormalK()` - Inconsistent**
**File:** `imaginable.py`  
**Issue:** Different naming schemes for similar operations

```python
# CURRENT (MIXED):
viewSagittal()             # Anatomical term
view2D()                   # Dimension term
getSliceNormalK()          # Index-based term
getSliceNormalJ()          
getSliceNormalI()          

# SHOULD BE (one approach):
view_sagittal()            / getSliceSagittal()
view_coronal()             / getSliceCoronal()
view_axial()               / getSliceAxial()

# OR (index-based throughout):
getSliceAlongAxis0()
getSliceAlongAxis1()
getSliceAlongAxis2()
```

**Fix Priority:** 🟡 MEDIUM

---

### E. GRAMMATICAL ISSUES (1 issue)

#### 14. **"whathappened()" - Grammatically Incorrect**
**File:** `imaginable.py` Line 593  
**Issue:** Should be "What happened" with camelCase

```python
# CURRENT (WRONG):
def whathappened(self):
    self.log.getWhatHappened()

# SHOULD BE:
def whatHappened(self):  # or getHistory()
    self.log.getWhatHappened()
```

**Fix Priority:** 🟡 LOW (internal method, uncommon use)

---

## SUMMARY TABLE

| # | Issue | File | Severity | Type |
|---|-------|------|----------|------|
| 1 | `getCenterOfGravitygetCenterOfGravityIndex` | imaginable.py | ✅ FIXED | Duplicate |
| 2 | `getMaskedNunmpyArray` | imaginable.py | 🔴 HIGH | Typo |
| 3 | `copythethreeinfosonandsetthemtimage` | imaginable.py | 🔴 HIGH | Typo + Case |
| 4 | `setSITKImageInforFromImage` | imaginable.py | 🟡 MEDIUM | Typo |
| 5 | README `LabelaMapable` + `balue` | README.md | 🟡 MEDIUM | Typo |
| 6 | `getImageSize()` / `changeImageSize()` | imaginable.py | 🟠 MED-HIGH | Inconsistency |
| 7 | `add()` / `multiply()` naming | imaginable.py | 🟡 MEDIUM | Inconsistency |
| 8 | `getBoundingBox(exclude=[0])` | imaginable.py | 🟡 MEDIUM | Unclear param |
| 9 | Double-underscore private methods | imaginable.py | 🟠 MED-HIGH | Convention |
| 10 | `filterValues()` semantics | imaginable.py | 🟡 MEDIUM | Unclear name |
| 11 | `dflt` abbreviation | imaginable.py | 🟡 MEDIUM | Abbreviation |
| 12 | `UD` / `LR` variables | imaginable.py | 🟡 MEDIUM | Abbreviation |
| 13 | "Imaginable" prefix inconsistency | imaginable.py | 🟡 MEDIUM | Inconsistency |
| 14 | View method naming | imaginable.py | 🟡 MEDIUM | Inconsistency |
| 15 | `whathappened()` capitalization | imaginable.py | 🟡 LOW | Grammar |

---

## RECOMMENDATIONS (Prioritized)

### PHASE 1 (Critical - Fix Immediately)
1. ✅ Fix `getCenterOfGravitygetCenterOfGravityIndex` → `getCenterOfGravityIndex` (DONE)
2. Fix `getMaskedNunmpyArray` → `getMaskedNumpyArray`
3. Fix `copythethreeinfosonandsetthemtimage` → `copyImageMetadata`
4. Fix `setSITKImageInforFromImage` → `setSITKImageInfoFromImage`
5. Update README: `LabelaMapable` → `LabelMapable` + `balue` → `value`

### PHASE 2 (Important - Fix Soon)
6. Rename private methods: `__name__()` → `_name()` (Python convention)
7. Standardize Get/Set naming: resolve `get*/change*` inconsistency
8. Add method aliases for backward compatibility:
   - `add()` → `addImage()` (or vice versa)
   - `multiply()` → `multiplyImage()`

### PHASE 3 (Nice to Have - Consider for v4)
9. Replace `dflt` → `default`
10. Clarify parameter names: `exclude` → `excludeValues` or `backgroundValue`
11. Document viewing methods (anatomical vs index-based)
12. Consider: `whathappened()` → `getHistory()`

---

## CODE STYLE OBSERVATIONS

### Positive ✅
- Consistent camelCase for methods
- Good use of verb-noun pairs (`get*`, `set*`)
- Descriptive parameter names in most cases
- Clear class hierarchy

### Needs Work ⚠️
- Inconsistent private method convention (use `_` not `__`)
- Mixed terminology (`set` vs `change` vs `add`)
- Some unclear abbreviations (`dflt`, `UD`, `LR`)
- Occasional typos in exported function names

---

## QUICK FIX SCRIPT

```python
# Recommended fixes for backward compatibility:

# In imaginable.py, add these aliases:

# Rename with aliases to maintain backward compatibility
getMaskedNumpyArray = getMaskedNunmpyArray  # NEW (correct)
copyImageMetadata = copythethreeinfosonandsetthemtimage  # NEW (clear)
setSITKImageInfoFromImage = setSITKImageInforFromImage  # NEW (fixed)

# Fix capitalization
whatHappened = whathappened  # NEW

# Add method aliases for consistency
def addImage(self, toadd):
    """Add image. See add() for details."""
    return self.add(toadd)

# Deprecation helpers
def setImageSize(self, newSize, *args, **kwargs):
    """Deprecated: use changeImageSize(). This is an alias."""
    return self.changeImageSize(newSize, *args, **kwargs)
```

---

## CONCLUSION

Your naming is **generally good with room for improvement**. The issues are:

1. **2 critical typos** that affect API clarity
2. **5 inconsistencies** in naming patterns
3. **Multiple opportunities** for standardization

**Recommendation:** Address critical typos now, standardize naming conventions in v4.0.

---

*Analysis Complete - November 21, 2025*
