# make a reader f the json files to subvide by category
import json
import glob
import boto3, json, sys
import base64
import numpy as np
import matplotlib.image as mpimg
import pyable_eros_montin.imaginable as ima
import pyable_eros_montin.utils as utu
import uuid
# Connect to Bedrock Runtime in your AWS region
session = boto3.Session(profile_name="nyu", region_name="us-east-1")
bedrock = session.client("bedrock-runtime")
import SimpleITK as sitk
# us.anthropic.claude-3-5-sonnet-20241022-v2:0


def sliceImages(nifti_file):
    A = ima.Imaginable(nifti_file)    
    A.dicomOrient('LPS')
    A.changeImageSpacing((1.0, 1.0, 1.0))
    S = A.getImageSize()
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(A.getImage())


    center_of_gravity = stats.GetCentroid(1)
    S=A.getIndexFromCoordinates(center_of_gravity)
    imagesavailable=[]    
    for u in range(3):
        for offset in [-10,0,10]:
            try:
                slice_index = S[u] + offset
                if slice_index < 0 or slice_index >= A.getImageSize()[u]:
                    continue
                im = utu.getImaginableSliceNumpy(A, u, slice_index)
            except:
                im=np.zeros((10,10))
                continue
            imagesavailable.append(im)
    return imagesavailable


def decode_dcm_to_nii_json_and_nii(json_file, nifti_file, model_id="us.anthropic.claude-3-5-sonnet-20241022-v2:0", verbose=False, temperature=0):
    """
    Classify MRI sequences using DICOM JSON metadata and optional images.
    Returns classification data and original DICOM metadata.
    """
    if "b1pre" in nifti_file.lower() or "b1map" in nifti_file.lower():
        if verbose:
            print("B1 map detected based on filename; skipping vision model.")
        return {}, json.load(open(json_file))
    # Load DICOM JSON metadata
    data = json.load(open(json_file))
    try:
        del data["Classification_vision"]
    except:
        pass
    addinfo={}
    processed_images = []
    A = ima.Imaginable(nifti_file)
    addinfo["ImageSizeITK"]=A.getImageSize()
    addinfo["ImageSpacingITK"]=A.getImageSpacing()
    
    A.dicomOrient('LPS')
    A.changeImageSpacing((1.0, 1.0, 1.0))
    S = A.getImageSize()
    # get the baricenter slices in each plane
    
    
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(A.getImage())

    # If your mask uses label 1
    center_of_gravity = stats.GetCentroid(1)
    S=A.getIndexFromCoordinates(center_of_gravity)
    orderslice=["saggital","coronal","axial"]
    imagesavailable=[]    
    for u in range(3):
        for offset in [-10,5,0,5,10]:
            try:
                slice_index = S[u] + offset
                if slice_index < 0 or slice_index >= A.getImageSize()[u]:
                    continue
                im = utu.getImaginableSliceNumpy(A, u, slice_index)
                if im.shape[0] < 50 or im.shape[1] < 50:
                    continue    
                arr = np.array(im)
                if arr.dtype.kind in ("U", "S", "O"):
                    arr = arr.astype(float)
                output_path = f"/tmp/{str(uuid.uuid4())}.png"
                mpimg.imsave(output_path, arr, cmap='gray')
                processed_images.append(output_path)
                imagesavailable
            except:
                continue
            
    if len(processed_images)==0:
        if verbose:
            print("No valid images extracted from NIfTI for vision model.")
            return {},data
    

    # Build textual prompt
    prompt = f"""
    You are an MRI expert.

    Task:
    Given this DICOM JSON metadata from dcm2nii and the associated images (if provided),
    classify the MRI characteristics and the anatomical region.

    Important rules:
    - Ignore any "bodypart "'Classification' field if present in the JSON.
    - Use BOTH the metadata and the visual content.
    - First think step-by-step internally about the anatomy, but DO NOT output your reasoning.
    - BodyRegion must be exactly one of: "HIP", "KNEE", "PELVIS", "OTHER", "UNKNOWN".
    - Do NOT guess. If confidence < 80%, output "UNKNOWN".

    Visual cues:
    - HIP/PELVIS: acetabulum, femoral head/neck, pelvic ring, gluteal muscles.
    - KNEE: femoral condyles, tibial plateau, patella, menisci usually two different structures separated by the void.
    - If anatomy is ambiguous or mostly soft tissue without bone landmarks → "UNKNOWN".
    - knees may show parts of the distal femur and proximal tibia only.
    
    For KNEE identification, look specifically for these anatomical landmarks:
    - Femoral condyles (rounded ends of the thigh bone)
    - Tibial plateau (flat top surface of the shin bone)
    - Patella (kneecap - small round/oval bone structure)
    - Menisci (two crescent-shaped structures with void/dark space between them)
    - Joint space between femur and tibia
    - May show only distal femur and/or proximal tibia portions

    For HIP/PELVIS identification, look for:
    - Acetabulum (hip socket)
    - Femoral head/neck (ball-shaped top of thigh bone)
    - Pelvic ring structures
    - Gluteal muscle masses

    Return ONLY a valid JSON object with:
    - SequenceType
    - SequenceName (based on classical naming, e.g., T1w, T2w, FLAIR, DWI, etc.)
    - AcquisitionDimension (2D/3D)
    - Plane
    - FatSuppression (True/False)
    - ParallelImaging (True/False)
    - BodyRegion: (KNEE/PELVIS/OTHER/UNKNOWN)
    - ContrastUsed (True/False)
    - UsableForDiagnosis (True/False)
    - ImageQuality (Poor/Fair/Good/Excellent)

    Input JSON ({json_file}):
    {json.dumps(data, indent=2)}
    images Provided separately as base64-encoded images ordered as {', '.join(orderslice)} planes with offsets -10,-5 0,+5, +10 from center of gravity slice.
    
"""



    # Prepare message content
    content = [{"type": "text", "text": prompt}]

    if processed_images:
        for image_file in processed_images:
            with open(image_file, "rb") as f:
                img_b64 = base64.b64encode(f.read()).decode("utf-8")
            content.append({
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": "image/png",
                    "data": img_b64
                }
            })

    # Build proper request body (Anthropic format)
    body = json.dumps({
        "anthropic_version": "bedrock-2023-05-31",
        "max_tokens": 400,
        "temperature": temperature,
        "messages": [{"role": "user", "content": content}]
    })

    # Call model
    response = bedrock.invoke_model(
        modelId=model_id,
        body=body
    )

    result = json.loads(response["body"].read())
    output_text = result["content"][0]["text"].strip()

    try:
        parsed = json.loads(output_text)
    except json.JSONDecodeError:
        # Recover if model added commentary
        start = output_text.find("{")
        end = output_text.rfind("}")
        parsed = json.loads(output_text[start:end+1])
    
    if verbose:
        print(json.dumps(parsed, indent=2))
    
    parsed["model_id"] = model_id
    parsed["version"] = "v1"
    parsed.update(addinfo)
    return parsed, data


def decode_dcm_to_nii_json(json_file,model_id="us.anthropic.claude-3-5-haiku-20241022-v1:0",verbose=False):
    # Load DICOM JSON metadata
    data = json.load(open(json_file))

    # Build the prompt
    prompt = f"""
You are an MRI expert. Based only on this dcm2nii JSON output metadata:
{json.dumps(data, indent=2)}
First think step-by-step internally about the anatomy, but DO NOT output your reasoning.

Return **only** a single valid JSON object — no text, no explanations, no comments — 
containing these keys:
- SequenceType
- SequenceName (based on classical naming, e.g., T1w, T2w, FLAIR, DWI, etc.)
- AcquisitionDimension (2D/3D)
- Plane
- FatSuppression (True/False)
- ParallelImaging (True/False)
- HASTE or setting image
- singleslice (True/False)
"""

    # Correct Bedrock request format for Anthropic models
    body = json.dumps({
        "anthropic_version": "bedrock-2023-05-31",
        "max_tokens": 400,
        "temperature": 0,
        "messages": [
            {"role": "user", "content": [{"type": "text", "text": prompt}]}
        ]
    })

    # Use the inference profile ID
    response = bedrock.invoke_model(
        modelId=model_id,
        body=body
    )

    result = json.loads(response["body"].read())
    output_text = result["content"][0]["text"].strip()

    try:
        parsed = json.loads(output_text)
    except json.JSONDecodeError:
        # recover if model added commentary
        start = output_text.find("{")
        end = output_text.rfind("}")
        parsed = json.loads(output_text[start:end+1])
    if verbose:
        print(json.dumps(parsed, indent=2))
    parsed["model_id"] = model_id
    parsed["version"] = "v0"
    return parsed,data

def read_json(json_file):
    parsed = json.load(open(json_file))
    return parsed

def get_sequence_enhanced_info(json_file, output_json=None, classification_json=None,nifti_file=None,FORCE=False):
    data = json.load(open(json_file))
    parsed =data.get("Classification", None)
    
    # if a classification already exists, return the full json
    if parsed and not FORCE:
        return data
    
    if classification_json and not parsed:
        parsed = read_json(classification_json)
        data["Classification"] = parsed
    else:
        if nifti_file:
            parsed, data = decode_dcm_to_nii_json_and_nii(json_file, nifti_file)
        else:
            parsed, data = decode_dcm_to_nii_json(json_file)
        data["Classification"] = parsed
    if output_json:
        with open(output_json, "w") as f:
            json.dump(data, f, indent=2)
    
    return data


def get_sequence_enhanced_info_2(json_file, output_json=None, classification_json=None,nifti_file=None,FORCE=False):
    data = json.load(open(json_file))
    parsed =data.get("Classification", None)
    
    
    # if a classification already exists, return the full json
    if parsed:
        del data["Classification"]
    with open(json_file, 'w') as file:
        json.dump(data, file, indent=2)
    
    if not 'Classification_json' in data.keys():
        parsed1, data = decode_dcm_to_nii_json(json_file)
        data["Classification_json"] = parsed1
        with open(json_file, 'w') as file:
            json.dump(data, file, indent=2)
    if not "Classification_vision" in data.keys():
        parsed, data = decode_dcm_to_nii_json_and_nii(json_file, nifti_file)
        data["Classification_vision"] = parsed
    
    if output_json:
        with open(output_json, "w") as f:
            json.dump(data, f, indent=2)
    
    return data


