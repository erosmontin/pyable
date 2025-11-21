import json
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from common import read_json, sliceImages
import os, sys, copy, pandas as pd
from glob import glob
st.markdown("""
    <style>
        .stMainBlockContainer {
            max-width: 2000px !important;
        }
        .block-container {
            max-width: 2000px !important;
        }
                /* Increase sidebar width */
        [data-testid="stSidebar"] {
            width: 800px !important;
            min-width: 800px !important;
        }
        [data-testid="stSidebar"] > div:first-child {
            width: 800px !important;
        }
        /* Make form text larger and bold */
        .stRadio > label {
            font-size: 18px !important;
            font-weight: bold !important;
        }
        .stRadio > div {
            font-size: 20px !important;
        }
        .stColumns {
            font-size: 20px !important;
            font-weight: bold !important;
        }
        /* Increase form correction section font size */
        [data-testid="stColumn"] {
            font-size: 20px !important;
        }
        [data-testid="stColumn"] > div > div > p {
            font-size: 20px !important;
            font-weight: bold !important;
        }
    </style>
    """, unsafe_allow_html=True)
# ===========================================================
# INPUT ROOT DIRECTORY
# ===========================================================
# Usage: streamlit run app.py /data/MYDATA/hip_mri/nifti
# ROOT contains patient directories, each with NIfTI + JSONs
# ===========================================================
ROOT = sys.argv[1]  # e.g. `/data/MYDATA/hip_mri/nifti/`
ROOT = ROOT.rstrip("/")

patient_dirs = sorted(
    [d for d in glob(os.path.join(ROOT, "*")) if os.path.isdir(d)]
)

st.title("MRI Batch Sequence Classification")

# Session state for skipped (faulty) JSONs
if "skip_jsons" not in st.session_state:
    st.session_state.skip_jsons = set()


# ===========================================================
# HELPERS
# ===========================================================
def get_output_dir(base_dir: str) -> str:
    """
    Mirror '.../nifti/...' into '.../corrected_jsons/...'.
    """
    return base_dir.replace("nifti", "corrected_jsons")


def _fs_label(seq_name: str, fat_value) -> str:
    """
    Build a column label from (SequenceName, FatSuppression).

    Examples:
        ('PD', True)   -> 'PD_FS'
        ('PD', False)  -> 'PD'
        ('T2w', True)  -> 'T2w_FS'
    """
    seq = str(seq_name) if seq_name is not None else "Unknown"
    fat_str = str(fat_value).strip().lower()

    is_fs = fat_str in ("true", "1", "yes", "fs", "fat", "fatsat", "fat sat")
    return f"{seq}_FS" if is_fs else seq


def compute_patient_summary(base_dir: str) -> pd.DataFrame:
    """
    For a given patient directory (base_dir), read ALL corrected JSONs
    in its OUTPUT_DIR and compute:

        # usable series per (SequenceName, FatSuppression)

    Returns a **wide** DataFrame with columns like:
        ['Patient', 'PD', 'PD_FS', 'T2w', 'T2w_FS', ...]

    Also writes per-patient CSV:
        OUTPUT_DIR/series_summary.csv
    """
    output_dir = get_output_dir(base_dir)
    os.makedirs(output_dir, exist_ok=True)

    patient_id = os.path.basename(base_dir.rstrip("/"))
    corrected_jsons = sorted(glob(os.path.join(output_dir, "*.json")))

    rows = []
    for cjson in corrected_jsons:
        try:
            info = read_json(cjson)
        except Exception:
            continue

        # Prefer corrected classification; fallback if needed
        cls = info.get("Classification")
        if cls is None:
            cls = info.get("Classification_json", {})

        # ONLY count usable series
        usable = bool(cls.get("UsableForDiagnosis", False))
        if not usable:
            continue

        seq_name = cls.get("SequenceName", "Unknown")
        fat = cls.get("FatSuppression", "Unknown")

        rows.append(
            {
                "Patient": patient_id,
                "SequenceName": str(seq_name),
                "FatSuppression": fat,
            }
        )

    if rows:
        df_raw = pd.DataFrame(rows)

        # long format: counts per (SequenceName, FatSuppression)
        df_long = (
            df_raw
            .groupby(["Patient", "SequenceName", "FatSuppression"])
            .size()
            .reset_index(name="UsableCount")
        )

        # wide format: one row per patient, columns per (SequenceName, FS)
        df_wide = df_long.pivot_table(
            index="Patient",
            columns=["SequenceName", "FatSuppression"],
            values="UsableCount",
            aggfunc="sum",
            fill_value=0,
        )

        # flatten MultiIndex columns -> 'PD', 'PD_FS', 'T2w_FS', ...
        df_wide.columns = [
            _fs_label(seq, fat)
            for (seq, fat) in df_wide.columns
        ]
        df_wide = df_wide.reset_index()

    else:
        # patient with no usable series yet
        df_wide = pd.DataFrame({"Patient": [patient_id]})

    # Save per-patient wide summary
    csv_path = os.path.join(output_dir, "series_summary.csv")
    df_wide.to_csv(csv_path, index=False)

    return df_wide


def update_global_summary(root: str, patient_dirs_list: list) -> pd.DataFrame:
    """
    Build/refresh a global wide CSV with ALL patients' statistics, by
    concatenating all per-patient `series_summary.csv` files
    (one row per patient).

    File:
        ROOT/global_series_summary.csv
    """
    all_rows = []
    for pdir in patient_dirs_list:
        out_dir = get_output_dir(pdir)
        csv_path = os.path.join(out_dir, "series_summary.csv")
        if os.path.exists(csv_path):
            try:
                dfp = pd.read_csv(csv_path)
                all_rows.append(dfp)
            except Exception:
                continue

    if all_rows:
        gdf = pd.concat(all_rows, ignore_index=True)

        # Ensure Patient is first column
        cols = ["Patient"] + [c for c in gdf.columns if c != "Patient"]
        gdf = gdf[cols]

        # Any missing combinations -> 0
        gdf = gdf.fillna(0)
    else:
        gdf = pd.DataFrame(columns=["Patient"])

    global_csv = os.path.join(root, "global_series_summary.csv")
    gdf.to_csv(global_csv, index=False)
    return gdf


@st.cache_data(show_spinner=False)
def load_slices(nii_path: str):
    """
    Cached wrapper around sliceImages to avoid re-loading
    the same NIfTI on every rerun.
    """
    return sliceImages(nii_path)


def is_corrected(json_path):
    try:
        data = read_json(json_path)
        return "Classification" in data
    except Exception:
        return False


# ===========================================================
# SELECT PATIENT DIRECTORY
# ===========================================================
st.sidebar.header("📁 Select Patient")

patient = st.sidebar.selectbox(
    "Choose a patient folder:",
    patient_dirs,
)

BASE_DIR = patient
OUTPUT_DIR = get_output_dir(BASE_DIR)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ===========================================================
# SIDEBAR: MODE SELECTION (NEW/REVIEW)
# ===========================================================
st.sidebar.markdown("---")
mode = st.sidebar.radio(
    "Select Mode:",
    ["Correct New", "Review & Edit"],
    horizontal=False
)

# ===========================================================
# FIND JSONS TO CORRECT FOR THIS PATIENT
# ===========================================================
json_files = sorted(glob(os.path.join(BASE_DIR, "*.json")))

# Corrected JSONs already saved for this patient
corrected_files = sorted(glob(os.path.join(OUTPUT_DIR, "*.json")))
corrected_names = {os.path.basename(p) for p in corrected_files}

# SKIP list for faulty NIfTIs
if "skip_jsons" not in st.session_state:
    st.session_state.skip_jsons = set()

if mode == "Correct New":
    # Final list of JSONs that still need correction
    remaining = []
    for f in json_files:
        fname = os.path.basename(f)
        if fname in st.session_state.skip_jsons:
            continue

        # skip if corrected JSON exists AND contains valid Classification
        corrected_path = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(corrected_path) and is_corrected(corrected_path):
            continue

        remaining.append(f)

    # Progress info
    total_jsons = len(json_files)
    remaining_count = len(remaining)
    corrected_count = total_jsons - remaining_count

    st.markdown(
        f"**Patient:** `{os.path.basename(BASE_DIR)}`  "
        f"— Corrected: **{corrected_count}/{total_jsons}**"
    )
    st.progress(corrected_count / total_jsons if total_jsons > 0 else 0.0)

    if remaining_count == 0:
        st.success(
            f"✅ All JSON files corrected for patient: **{os.path.basename(BASE_DIR)}**"
        )
        st.stop()

    # ===========================================================
    # SELECT JSON FILE TO CORRECT
    # ===========================================================
    JSON = st.selectbox("Select JSON to correct:", remaining)
    save_path = os.path.join(OUTPUT_DIR, os.path.basename(JSON))

    INFO = read_json(JSON)
    L = INFO.get("Classification_json", {})
    V = INFO.get("Classification_vision", {})

else:  # Review & Edit mode
    st.markdown(f"**Patient:** `{os.path.basename(BASE_DIR)}` — **Review Mode**")
    
    if len(corrected_files) == 0:
        st.warning("No corrected JSON files available for this patient yet.")
        st.stop()

    # Select from corrected files
    corrected_basenames = [os.path.basename(f) for f in corrected_files]
    selected_corrected = st.selectbox("Select corrected JSON to review:", corrected_basenames)
    JSON = os.path.join(OUTPUT_DIR, selected_corrected)
    save_path = JSON

    INFO = read_json(JSON)
    L = INFO.get("Classification", {})
    V = {}  # No Vision classification in review mode

# ===========================================================
# LOAD IMAGES
# ===========================================================
# NIfTI files are always in BASE_DIR, not in OUTPUT_DIR
nii = os.path.join(BASE_DIR, os.path.basename(JSON).replace(".json", ".nii"))
if not os.path.exists(nii):
    nii = os.path.join(BASE_DIR, os.path.basename(JSON).replace(".json", ".nii.gz"))

try:
    images = load_slices(nii)
except Exception as e:
    st.error(f"⚠️ Error loading NIfTI file:\n{nii}\n{e}")
    st.warning("Skipping this file due to error.")
    if mode == "Correct New":
        st.session_state.skip_jsons.add(os.path.basename(JSON))
    st.rerun()

# ===========================================================
# KEYS TO CORRECT
# ===========================================================
all_keys = sorted(set(L.keys()).union(V.keys()))
SKIP_KEYS = [
    "model_id",
    "version",
    "ImageSpacingITK",
    "ImageSizeITK",
    "ParallelImaging",
    "Classification_json",
    "Classification_vision",
    "singleslice",
    "AcquisitionDimension",
    "ContrastUsed",
]
all_keys = [k for k in all_keys if k not in SKIP_KEYS]

# ===========================================================
# SIDEBAR: IMAGE PREVIEW + JSON PREVIEW
# ===========================================================
with st.sidebar:
    st.header("🖼️ Preview Slices")
    fig, axes = plt.subplots(3, 3, figsize=(8, 8))
    for ax, img in zip(axes.flatten(), images):
        ax.imshow(img, cmap="gray")
        ax.axis("off")
    plt.tight_layout()
    st.pyplot(fig)

    st.header("📄 Live JSON Preview")
    json_preview_box = st.empty()

# ===========================================================
# MAIN FORM: INTERACTIVE CORRECTION
# ===========================================================
st.subheader(f"{'Reviewing' if mode == 'Review & Edit' else 'Correcting'}: **{os.path.basename(JSON)}**")
final = {}
for key in all_keys:
    val_L = L.get(key, "")
    val_V = V.get(key, "")

    col1, col2, col3, col4, col5 = st.columns([2, 2, 2, 3, 2])

    with col1:
        st.markdown(f"**{key}**")
    with col2:
        st.write(f"L: `{val_L}`")
    with col3:
        if mode == "Correct New":
            st.write(f"V: `{val_V}`")

    is_bool = isinstance(val_L, bool) or isinstance(val_V, bool)

    if is_bool:
        if mode == "Correct New":
            options = [f"L ({val_L})", f"V ({val_V})", "True", "False"]
        else:
            options = ["True", "False"]
    else:
        if mode == "Correct New":
            options = [f"L ({val_L})", f"V ({val_V})", "Custom"]
        else:
            options = [f"Current ({val_L})", "Custom"]

    # Remove empty options
    options = [opt for opt in options if not (opt.endswith("()") or opt.endswith("('')"))]

    # Ensure we have at least one option
    if not options:
        options = ["Custom"]

    # Determine default index
    default_idx = 0
    if mode == "Correct New":
        if val_L != val_V:
            for i, opt in enumerate(options):
                if opt.startswith("V"):
                    default_idx = i
                    break
    else:
        default_idx = 0

    # Adjust default index if it exceeds available options
    if default_idx >= len(options):
        default_idx = 0

    with col4:
        choice = st.radio(
            label=f"choice_for_{key}",
            options=options,
            index=default_idx,
            horizontal=True,
            key=f"radio_{key}",
            label_visibility="collapsed",
        )

    with col5:
        if choice.startswith("L") or choice.startswith("Current"):
            final[key] = val_L
        elif choice.startswith("V"):
            final[key] = val_V
        elif is_bool:
            final[key] = True if choice == "True" else False
        else:
            final[key] = st.text_input("", key=f"{key}_custom")

# ===========================================================
# UPDATE PREVIEW JSON (LIVE)
# ===========================================================
preview = copy.deepcopy(INFO)

# Process final values
final_classification = {}
for key in all_keys:
    val_L = L.get(key, "")
    val_V = V.get(key, "")
    
    if mode == "Correct New":
        # If L and V are identical, automatically use V
        if val_L == val_V:
            final_classification[key] = val_V
        else:
            # Otherwise use the user's choice
            final_classification[key] = final[key]
    else:
        # In review mode, use user's choice
        final_classification[key] = final[key]

preview["Classification"] = final_classification
preview.pop("Classification_json", None)
preview.pop("Classification_vision", None)
with st.sidebar:
    json_preview_box.json(preview)

# ===========================================================
# SAVE + UPDATE PER-PATIENT AND GLOBAL STATISTICS
# ===========================================================
button_label = "💾 Save Changes" if mode == "Review & Edit" else "💾 Save and Continue"
if st.button(button_label):
    # Save corrected JSON
    with open(save_path, "w") as f:
        json.dump(preview, f, indent=2)

    success_msg = "Saved changes to JSON!" if mode == "Review & Edit" else "Saved corrected JSON"
    st.success(f"{success_msg}: {save_path}")

    # Recompute this patient's summary from ALL corrected JSONs (wide)
    patient_df = compute_patient_summary(BASE_DIR)

    # Update global summary
    global_df = update_global_summary(ROOT, patient_dirs)

    st.success("📊 Statistics updated (patient + global).")
    if mode == "Correct New":
        st.info("Moving to next JSON...")
    st.rerun()

# ===========================================================
# DISPLAY STATISTICS TABLES
# ===========================================================
# st.markdown("---")
# st.header("📊 Statistics")

# col1, col2 = st.columns(2)

# with col1:
#     st.subheader("👤 Current Patient Summary")
#     patient_df = compute_patient_summary(BASE_DIR)
#     st.dataframe(patient_df, use_container_width=True)

# with col2:
#     st.subheader("🏥 Global Summary (All Patients)")
#     global_df = update_global_summary(ROOT, patient_dirs)
#     st.dataframe(global_df, use_container_width=True)
