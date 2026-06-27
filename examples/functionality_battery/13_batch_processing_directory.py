from _common import DATA_DIR, output_dir, save_synthetic_dataset, write_json, print_created
from pyable import processImageDirectory


save_synthetic_dataset()
out = output_dir("13_batch_processing_directory")


def summarize_image(image):
    return {
        "size": image.getImageSize(),
        "spacing": image.getImageSpacing(),
        "mean": image.getMeanValue(),
        "max": image.getMaximumValue(),
    }


csv_path = out / "batch_summary.csv"
dataframe = processImageDirectory(
    str(DATA_DIR),
    summarize_image,
    file_pattern="synthetic_[ilr]*.nii.gz",
    output_csv=str(csv_path),
    recursive=False,
    verbose=True,
)

json_path = write_json(
    out / "batch_summary_preview.json",
    dataframe.head().to_dict(orient="records"),
)

print_created(csv_path, json_path)