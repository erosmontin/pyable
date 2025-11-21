import os
import json
import glob
import tqdm
from common import get_sequence_enhanced_info_2
def read_category_files(directory):
    json_files = glob.glob(f"{directory}/*.json")
    OUT=[]
    for file_path in tqdm.tqdm(json_files):
        try:
            with open(file_path, 'r') as file:
                fn=file_path.replace('.json', '.nii')
                if not os.path.exists(fn):
                    fn=file_path.replace('.json', '.nii.gz')
                    if not os.path.exists(fn):
                        fn=None
                classified_data = get_sequence_enhanced_info_2(file_path,file_path,nifti_file=fn,FORCE=True)                
                classified_data["file"] = {"json": file_path, "nifti": fn}
                OUT.append(classified_data)
                
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    return OUT

            
                
                
            



if __name__ == "__main__":
    import pandas as pd
    out =pd.DataFrame()
    
    import glob
    
    D=glob.glob("/data/MYDATA/hip_mri/nifti/*")
    out=pd.DataFrame()    
    for directory in D:
        print(f"Processing directory: {directory}")
        o=[]
        for c in read_category_files(directory):
            togo=c["Classification_vision"]
            togo["Classification_json"]=c["Classification_json"]
            togo["file"]=c["file"]
            o.append(togo) 
        out = pd.concat([out,pd.json_normalize(o)], ignore_index=True, axis=0)
    
    out.to_csv("classified_sequences.csv", index=False)    
    
