from _common import DATA_DIR, print_created, save_synthetic_dataset, write_json


paths = save_synthetic_dataset()
summary = write_json(
    DATA_DIR / "dataset_manifest.json",
    {key: str(path) for key, path in paths.items()},
)

print_created(*paths.values(), summary)