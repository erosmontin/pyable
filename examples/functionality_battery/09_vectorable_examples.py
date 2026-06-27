from _common import make_vector_field, output_dir, print_created, write_json


vector = make_vector_field()
out = output_dir("09_vectorable_examples")

magnitude = vector.getMagnitude()
x_component = vector.getComponent(0)
scaled = vector.getDuplicate().scaleVector([2.0, 1.0, 0.5])
normalized = vector.getDuplicate().normalize(target_magnitude=1.0)
smoothed = vector.getDuplicate().applyGaussianSmoothing(sigma=1.0)

paths = [
    magnitude.writeImageAs(str(out / "vector_magnitude.nii.gz")),
    x_component.writeImageAs(str(out / "vector_x_component.nii.gz")),
    scaled.writeImageAs(str(out / "vector_scaled.mha")),
    normalized.writeImageAs(str(out / "vector_normalized.mha")),
    smoothed.writeImageAs(str(out / "vector_smoothed.mha")),
]

summary_path = write_json(
    out / "vector_summary.json",
    {
        "number_of_components": vector.getNumberOfComponents(),
        "mean_vector": vector.getMeanVector(),
        "magnitude_statistics": vector.getVectorStatistics(),
        "scaled_mean_vector": scaled.getMeanVector(),
        "normalized_magnitude_statistics": normalized.getVectorStatistics(),
    },
)
paths.append(summary_path)

print_created(*paths)