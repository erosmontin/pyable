# Interesting Utilities from d2/ Directory

## Summary

The `d2/` directory contains several useful utilities that could enhance the pyable classes:

---

## 1. **Image Slicing Utilities** (`common.py`)

### `sliceImages(nifti_file)` - Extract Multiple Slices from Different Planes

**Purpose**: Extract representative 2D slices from a 3D volume around the center-of-gravity.

**Functionality**:
- Orients image to standard LPS
- Resamples to isotropic spacing (1.0, 1.0, 1.0)
- Computes center-of-gravity using label statistics
- Extracts slices at ±10, 0, +10 mm offsets
- Returns slices in sagittal, coronal, and axial planes

**Potential Integration**:
```python
# Add to Imaginable class
def getRepresentativeSlices(self, offsets=[-10, 0, 10]):
    """Extract representative slices from 3 orthogonal planes"""
    # Could use for quick preview or batch processing
```

**Current Usage**:
```python
images = sliceImages('mri_scan.nii.gz')
# Returns list of 9 2D numpy arrays (3 planes × 3 offsets)
```

---

## 2. **DICOM JSON Classification** (`common.py`)

### `decode_dcm_to_nii_json_and_nii()` - AI-Powered Sequence Classification

**Purpose**: Classify MRI sequences using both DICOM JSON metadata AND vision models (AWS Bedrock Claude).

**Features**:
- Reads DICOM dcm2niix JSON files
- Extracts slices from 3 planes with context
- Encodes images as PNG for vision model
- Uses Claude 3.5 Sonnet for intelligent classification
- Combines JSON metadata + visual inspection

**AI Integration**:
- AWS Bedrock API for vision model inference
- Base64 image encoding
- Temperature/parameter control

**Potential Integration**:
```python
# Could add to Imaginable as:
def classifySequence(self, dicom_json_file, use_ai=True):
    """Classify MRI sequence using metadata + optional AI vision model"""
    # Vision model could identify:
    # - Sequence type (T1, T2, PD, FLAIR, etc.)
    # - Anatomical coverage
    # - Image quality
    # - Artifacts
```

---

## 3. **Batch Processing** (`read_dir.py`)

### `read_category_files(directory)` - Process Multiple Studies

**Purpose**: Batch process directories of DICOM JSON + NIfTI pairs.

**Features**:
- Walks directory structure
- Matches JSON with corresponding NIfTI files (.nii or .nii.gz)
- Applies classification to each
- Aggregates results to CSV

**Workflow**:
```python
# Process entire dataset
D = glob.glob("/data/MYDATA/hip_mri/nifti/*")
out = pd.DataFrame()

for directory in D:
    results = read_category_files(directory)
    out = pd.concat([out, pd.json_normalize(results)])

out.to_csv("classified_sequences.csv")
```

**Potential Integration**:
```python
# Add to Imaginable/utils:
def batchProcessDirectory(directory, processor_func, output_format='csv'):
    """Apply processor to all images in directory, aggregate results"""
```

---

## 4. **Interactive Jupyter Widgets** (`subdivide_sequences_data.py`)

### Interactive Sequence Review UI

**Purpose**: Manual verification/correction of automated classifications.

**Features**:
- Display 3×3 grid of representative slices
- Show automated classifications from JSON and vision model
- Interactive form for manual correction
- Side-by-side comparison (JSON vs Vision model predictions)

**UI Components**:
```python
# 3×3 image grid with matplotlib
fig, axes = plt.subplots(3, 3, figsize=(10, 10))

# Interactive widgets
- Text fields for sequence names
- Checkboxes for boolean properties
- Dropdowns for categorical choices
- Side-by-side labels showing automated predictions
```

**Potential Integration**:
```python
# Could add to plotting system:
def interactiveReview(self, metadata_json=None, predictions=None):
    """Display interactive review UI for QC"""
    # Show slices + form for corrections
    # Export corrected metadata
```

---

## 5. **Streamlit Web App** (`app.py`)

### Batch Classification UI

**Purpose**: Web-based UI for batch processing and manual review.

**Features**:
- Streamlit-based web interface
- Patient-level batch processing
- Dynamic form generation from metadata
- Mirroring directory structure (nifti → corrected_jsons)
- Session state management
- Summary statistics per patient

**Architecture**:
```
Patient Directory Structure:
├── patient_001/
│   ├── sequence_001.nii.gz
│   ├── sequence_001.json
│   └── sequence_002.nii.gz
├── patient_002/
│   └── ...

Corrected Output:
├── patient_001/
│   ├── corrected_jsons/
│   │   ├── sequence_001.json (corrected)
│   │   └── sequence_002.json (corrected)
```

**Potential Integration**:
```python
# Could create a Streamlit app for:
def streamlitViewer(base_directory, output_directory):
    """Web UI for batch image review and correction"""
    # Display images in grid
    # Interactive form for metadata correction
    # Export corrected data
```

---

## Recommended Integration Strategies

### **Strategy 1: Minimal - Add Utility Methods**
Add key functions to existing classes:
```python
# Add to Imaginable
def getRepresentativeSlices(self, num_offsets=3)
def classifyFromDICOMJSON(self, json_file)
```

### **Strategy 2: Moderate - Create New Utility Module**
Create `pyable/batch.py`:
```python
class BatchProcessor:
    def processDirectory(self, directory, processor_func)
    def aggregate_results(self, results_list)
    def export_to_csv(self, data, output_path)

class InteractiveReviewer:
    def show_slices_grid(self)
    def create_review_form(self)
    def export_corrections(self)
```

### **Strategy 3: Full - Add Web Interface**
Create optional Streamlit app:
```
pyable/apps/streamlit_viewer.py
```

---

## Key Insights from d2/

| Component | Value | Integration Level |
|-----------|-------|-------------------|
| **Slice extraction** | Representative sampling | Easy ✅ |
| **AI classification** | Vision model integration | Medium 🟡 |
| **Batch processing** | Workflow automation | Easy ✅ |
| **Interactive QC** | Manual verification UI | Medium 🟡 |
| **Web app** | Production interface | Hard 🔴 |

---

## Recommended Quick Wins

### 1. **Slice Extraction** (10 min)
```python
# Add to Imaginable
def extractRepresentativeSlices(self, planes=['sagittal', 'coronal', 'axial'], 
                                offsets=[-10, 0, 10]):
    """Get representative slices from center of gravity"""
    # From sliceImages() function
```

### 2. **Batch Processing** (20 min)
```python
# Add to utils
def processImageDirectory(directory, processor, file_pattern='*.nii.gz',
                         output_format='csv'):
    """Apply processor to all images in directory"""
    # From read_dir.py
```

### 3. **Interactive Grid Viewer** (15 min)
```python
# Extend plotable.py
class GridPlotter:
    def show_slices_grid(self, slices, titles=None)
    # From subdivide_sequences_data.py
```

---

## Summary

The d2/ directory contains **production-ready utilities** for:

✅ **Multi-plane slice extraction** - Quick preview of 3D volumes  
✅ **Batch processing workflows** - Handle large image datasets  
✅ **AI-assisted classification** - Vision model integration  
✅ **Interactive review UI** - Quality control and corrections  
✅ **Web interface** - Production deployment  

**Recommended approach**: Extract the most useful functions (slicing, batch processing) and integrate them into pyable core classes as optional utilities.

---

**Status**: Analysis Complete  
**Recommendation**: Start with slice extraction + batch processing (Easy wins)  
**Follow-up**: Consider web app for production deployment
