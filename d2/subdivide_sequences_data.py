

import json
from common import read_json,sliceImages
import sys
import matplotlib.pyplot as plt
from ipywidgets import VBox, HBox, Text, Dropdown, Checkbox, Label, Output, Layout
import ipywidgets as widgets
import numpy as np
import os
from IPython.display import display, HTML, clear_output

JSON='/data/MYDATA/hip_mri/nifti/preop_217_dei/preop_217_dei_AX_PD_OBL_UNILAT_20150629081560_8.json'
dcm2nii_json = JSON

INFO=read_json(dcm2nii_json)

if "Classification" in INFO.keys():
    parsed = INFO["Classification"]
else:
    L=INFO["Classification_json"]
    V=INFO["Classification_vision"]


nifti_file = dcm2nii_json.replace('.json', '.nii')
if not os.path.exists(nifti_file):
    nifti_file = dcm2nii_json.replace('.json', '.nii.gz')
    if not os.path.exists(nifti_file):
        nifti_file = None

images = sliceImages(nifti_file)
# Merge all keys from both
all_keys = sorted(set(L.keys()).union(V.keys()))

# -----------------------------------------------------------
# 3. DISPLAY THE 3×3 GRID OF IMAGES
# -----------------------------------------------------------
# ---- 3. DISPLAY IMAGES INSIDE A WIDGET OUTPUT ----
img_out = Output()

with img_out:
    fig, axes = plt.subplots(3, 3, figsize=(10, 10))
    for ax, img in zip(axes.flatten(), images):
        ax.imshow(img, cmap="gray")
        ax.axis("off")
    plt.tight_layout()
    plt.show()


# ---- 4. BUILD THE FORM WIDGETS ----
def make_row(key):
    val_L = L.get(key, "")
    val_V = V.get(key, "")

    # placeholder logic
    placeholder = str(val_L) if val_L == val_V else "choose..."

    # widget type logic
    if isinstance(val_L, bool) or isinstance(val_V, bool):
        user_widget = Checkbox(value=(val_L == val_V and isinstance(val_L, bool)))
    elif isinstance(val_L, list) or isinstance(val_V, list):
        user_widget = Text(value="", placeholder=placeholder)
    else:
        user_widget = Text(value="", placeholder=placeholder)

    row = HBox([
        Label(value=key, layout=Layout(width="200px")),
        user_widget,
        Label(value="L:", layout=Layout(width="30px")),
        Label(value=str(val_L), layout=Layout(width="200px")),
        Label(value="V:", layout=Layout(width="30px")),
        Label(value=str(val_V), layout=Layout(width="200px")),
    ], layout=Layout(padding="4px"))

    return row


form_rows = [make_row(k) for k in all_keys]
form_box = VBox(form_rows)


# ---- 5. COMBINE IMAGES + FORM IN ONE DISPLAY ----
final_ui = VBox([
    img_out,        # top: the 3×3 grid
    Label(value="Review and correct the sequence classification:"), 
    form_box        # bottom: the L/V form
])

display(final_ui)