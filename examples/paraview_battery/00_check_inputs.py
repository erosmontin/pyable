from _common import IMAGE_PATH, ROI_PATH, LABELMAP_PATH
from _common import ensure_inputs, load_image, load_roi, load_labelmap


ensure_inputs()

image = load_image()
roi = load_roi()
labelmap = load_labelmap()

print("Inputs are available.")
print("Image:", IMAGE_PATH)
print("ROI:", ROI_PATH)
print("Label map:", LABELMAP_PATH)
print("Image size:", image.getImageSize())
print("Image spacing:", image.getImageSpacing())
print("ROI size:", roi.getImageSize())
print("Label values:", sorted(int(v) for v in labelmap.getImageUniqueValues(exclude=[0])))