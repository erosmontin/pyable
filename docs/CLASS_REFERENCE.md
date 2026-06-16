# pyable Class Reference

This document is the package-level reference for the public classes in `pyable`.

## Conventions

- Most mutating methods return `self`.
- `Imaginable.getImageAsNumpy()` returns arrays in `(z, y, x)` order in v3.
- `Roiable` and `LabelMapable` default to nearest-neighbor interpolation for resampling.
- Many high-level methods delegate to the functional modules `pyable.deformations`, `pyable.segmentation`, and `pyable.metrics`.

## Imaginable

**Constructor**

- `Imaginable(filename=None, image=None, verbose=False)`

**Purpose**

- Base mutable wrapper for a `SimpleITK.Image`.
- Keeps an image stack for undo/reset behavior.
- Central place for I/O, geometry edits, transforms, array conversion, math, and quick visualization.

**Key attributes**

- `imageStack`: internal history stack used by `setImage()`, `undo()`, and `reset()`.
- `log`: operation log used by `whathappened()`.
- `dfltInterpolator`: default resampling interpolator. Scalar images start with `sitkLinear`.
- `dfltuseNearestNeighborExtrapolator`: default extrapolator flag for resampling.
- `verbose`: prints operations when enabled.
- `InputFileName`: remembered source path when created from disk.
- `settings["spacingMinSize"]`: rounding precision used when spacing is recomputed.

**Lifecycle and I/O**

- `getImage() -> sitk.Image`: return the current image from the stack.
- `setImage(p, w=None) -> self`: push a new image onto the stack and optionally log a message.
- `reset() -> None`: undo back to the first image in the stack.
- `undo() -> self`: pop one image from the stack.
- `isImageSet() -> bool`: report whether a valid image is available.
- `setVerbose(v)`, `getVerbose()`: configure/read logging verbosity.
- `getInputFileName() -> str | None`, `setInputFileName(fn)`: manage the source filename.
- `writeImageAs(filename, force=True) -> str`: write the current image to disk.
- `write(filename, force=True) -> str`: alias for `writeImageAs`.
- `whathappened() -> None`: print the stored operation log.
- `describe() -> dict`: print and return a compact image summary.
- `printImageInfo() -> dict`: print image metadata fields.
- `forkDuplicate() -> Imaginable`: deep copy of the object.
- `getDuplicate() -> Imaginable`: same class, same current image, fresh object.

**Array and image conversion**

- `getImageAsNumpy() -> np.ndarray`: return the image in `(z, y, x)` order.
- `getImageAsNumpyZYX() -> np.ndarray`: explicit alias for `(z, y, x)`.
- `getImageAsNumpyXYZ() -> np.ndarray`: return the array in legacy `(x, y, z)` order.
- `getImageAsNumpyForPyTorch() -> np.ndarray`: return a tensor-friendly layout.
- `setImageFromNumpy(nparray, refimage=None, vector=False, spacing=None, origin=None, direction=None) -> self`: load from NumPy in `(z, y, x)` order.
- `setImageFromNumpyZYX(...) -> self`: explicit alias for `setImageFromNumpy`.
- `setImageFromNumpyXYZ(...) -> self`: load from `(x, y, z)` order.
- `getITKImage() -> sitk.Image`: return the wrapped SimpleITK object.
- `getVTKImage() -> vtk.vtkImageData`: convert to VTK image data.

**Geometry, coordinates, and orientation**

- `getImageDirection()`, `setImageDirection(direction)`: raw direction matrix access.
- `getDirectionCosines()`, `setDirectionCosines(direction)`: semantic aliases for direction access.
- `getImageSpacing()`, `setImageSpacing(spacing)`: raw spacing access.
- `getImageOrigin()`, `setImageOrigin(origin)`: raw origin access.
- `getImageDimension() -> int`: image dimensionality.
- `getImageNumberOfComponentsPerPixel() -> int`: scalar/vector channel count.
- `getImagePixelTypeAsString() -> str`, `getImagePixelTypeAsID() -> int`: pixel type inspection.
- `getImageSize(index=None) -> tuple | int`: full size or a single dimension length.
- `getImageCenterIndex() -> list[int]`, `getImageCenterCoordinate() -> list[float]`: center in index or physical space.
- `getCoordinatesFromIndex(P) -> tuple`: legacy physical-point lookup from ITK index.
- `getIndexFromCoordinates(I) -> tuple`: legacy ITK index lookup from physical point.
- `getPhysicalPointFromArrayIndex(kji_index) -> tuple`: `(z, y, x)` NumPy index to physical point.
- `getArrayIndexFromPhysicalPoint(xyz_point) -> tuple`: physical point to `(z, y, x)` array index.
- `getPhysicalPointFromITKIndex(ijk_index) -> tuple`: ITK index to physical point.
- `getITKIndexFromPhysicalPoint(xyz_point) -> tuple`: physical point to ITK index.
- `getCornersCoordinates() -> list[tuple]`: physical coordinates of image corners.
- `isInsidePoint(P) -> bool`, `isInsideIndex(V) -> bool`: inclusion tests.
- `getOrientationCode() -> str`: current anatomical orientation code.
- `dicomOrient(orientation="LPS") -> self`: reorient via DICOM orientation rules.
- `reorientToLPS()`, `reorientToRAS()`, `reorientToRPI() -> self`: named orientation helpers.
- `isAxisAligned(tolerance=1e-6) -> bool`: check for near-identity direction cosines.

**Resampling and transforms**

- `resampleOnCanonicalSpace(interpolator=None, useNearestNeighborExtrapolator=None, bgvalue=0.0) -> self`: move the image to canonical axis-aligned space.
- `changeImageSpacing(spacing, interpolator=None, useNearestNeighborExtrapolator=None, bgvalue=0.0) -> self`: resample to a new spacing.
- `changeImageDirection(direction, interpolator=None, useNearestNeighborExtrapolator=None, bgvalue=0.0) -> self`: resample to a new direction matrix.
- `resampleToAxisAligned(interpolator=None, useNearestNeighborExtrapolator=None, bgvalue=0.0) -> self`: axis-align oblique images.
- `resampleOnTargetImage(target, interpolator=None, default_value=0, useNearestNeighborExtrapolator=None) -> self`: match another image geometry.
- `changeImageSize(newSize, interpolator=None, bgvalue=0.0, useNearestNeighborExtrapolator=None) -> self`: resample to a new voxel grid size.
- `padImage(lower_padding, upper_padding, padding_value=0) -> self`: constant padding.
- `getPaddedImage(padding, padding_value=0) -> Imaginable`: symmetric padded copy.
- `cropImage(lowerB, upperB, coordinates=None) -> self`: extract a subvolume.
- `cropToBoundingBox() -> self`: crop to the non-background extent.
- `translateImage(T, interpolator=sitk.sitkLinear, reference_image=None, default_value=0, useNearestNeighborExtrapolator=False) -> self`: rigid translation.
- `scaleImage(S, center=None, centerindex=False, interpolator=None, reference_image=None, default_value=0, useNearestNeighborExtrapolator=None) -> self`: scale about a center.
- `rotateImage(rotation=None, center=None, centerindex=False, translation=None, interpolator=None, reference_image=None, default_value=0.0, useNearestNeighborExtrapolator=None, angle=None) -> self`: Euler rotation. `angle` is a backward-compatible in-plane shortcut.
- `transform(T, interpolator=None, reference_image=None, default_value=0, useNearestNeighborExtrapolator=None) -> self`: generic transform application using inverse resampling logic.
- `transformImageAffine(A, translation=None, center=[], centerindex=False, interpolator=None, reference_image=None, default_value=0.0, useNearestNeighborExtrapolator=None) -> self`: explicit affine matrix transform.
- `applyTransform(transform, target_image=None, interpolator=None, default_value=0) -> self`: registration transform with automatic label-preserving behavior for integer images.
- `applyDisplacementField(displacement_field, target_image=None, interpolator=None, default_value=0) -> self`: warp by vector field.
- `warpImage(displacement_field, **kwargs) -> self`: alias for `applyDisplacementField`.
- `alignGeometry(reference_image) -> self`: copy/resample to a reference geometry.
- `invertDisplacementField(max_iterations=100, mean_error_tolerance=0.001) -> Fieldable`: invert a deformation field.
- `convertTransformToField(transform) -> Fieldable`: sample a transform on the current grid.
- `composeTransforms(transforms, inverse_flags=None) -> self`: build and apply a composite transform.
- `transformFromRegitration(...) -> self`: legacy transform helper kept for compatibility.

**Pixel type, math, and statistics**

- `changePixelType(dtype) -> self`, `cast(dtype) -> self`: cast image type.
- `getPossiblePixelTypes() -> list`: report accepted type names.
- `applyAbs() -> self`, `applyModulus() -> self`: absolute value / complex modulus.
- `add(toadd)`, `subtract(toadd)`, `multiply(toadd)`, `divide(toadd) -> self`: image or scalar arithmetic.
- `addImage(toadd)`, `subtractImage(toadd)`, `multiplyImage(toadd)`, `divideImage(toadd) -> self`: backward-compatible aliases.
- `filterValues(values) -> self`: keep/filter selected scalar values.
- `getImageUniqueValues(exclude=[]) -> list`: unique voxel values.
- `getMaximumValue()`, `getMinimumValue()`, `getMeanValue()`, `getVariance()`, `getStdValue()`, `getSum()`: scalar summary statistics.
- `getNumberOfNonZeroVoxels() -> int`: count of non-zero voxels.
- `getVoxelVolume() -> float`, `getNumberOfVoxels() -> int`, `getNumberofVoxels() -> int`, `getVolume() -> float`: voxel count and physical volume.
- `getRoiableValuesUpper(th) -> list`: convenience thresholded value extraction.
- `getValuesInRoi(roi) -> np.ndarray`: values inside a mask.
- `getBoundingBox(exclude=[0]) -> tuple`: bounding box of non-excluded voxels.

**Slices and quick viewers**

- `getSliceNormalKAsNumpy(slice)`, `getSliceNormalJAsNumpy(slice)`, `getSliceNormalIAsNumpy(slice) -> np.ndarray`: extract orthogonal slices as arrays.
- `getSliceNormalK(slice)`, `getSliceNormalJ(slice)`, `getSliceNormalI(slice) -> Imaginable`: same slice extraction as wrapped images.
- `viewK(km=[True, True])`, `viewJ(km=[True, True])`, `viewI(km=[True, True])`: legacy plane viewers.
- `viewAxial()`, `viewCoronal()`, `viewSagittal()`, `view2D()`: quick inspection viewers.
- `overlayAble(secondimaginable, axis, index, ..., titles=None, figsize=None, colorbar=None, index_mode="auto") -> matplotlib object`: 2D overlay helper. Pass list/array `axis` and/or `index` values to render a compact grid; multi-axis plus a 3D `index` point uses `index[axis]` by default. Use `index_mode="cartesian"` to combine every axis with every index.
- `overlayAbleImage(secondimaginable, axis, index, ..., title=None, titles=None, ncols=None, slice_offsets=None, as_base64=False, data_uri=False, save=None, index_mode="auto") -> np.ndarray | str`: image+overlay raster only, returned as RGBA pixels or PNG base64. Pass multiple axes and a 3D `index` point for orthogonal planes, or use `slice_offsets` for a tight 2.5D montage.
- `overlayReport(overlay, spacing=None, orientation="LPS", views="all", slice_offsets=None, fill_alpha=0.2, contour_alpha=0.95, overlay_color=(1, 0, 0), contour_iterations=2, image_cmap="gray", figsize=None, title=None, show=False, save=None, dpi=150, stats=True) -> matplotlib figure`: publication-style report.
- `plotOverlay(overlay=None, alpha=0.5, title=None, slice_idx=None, **kwargs) -> PlotViewer`: plotting wrapper.
- `viewInteractive(overlays=None, orientation=2, slice_idx=None, title=None, figsize=(14, 10), cmap="gray") -> InteractiveViewer`: interactive viewer entry point.
- `extractRepresentativeSlices(planes="all", offsets=[-10, 0, 10], verbose=False) -> dict`: standard orthogonal slice extraction around the center.
- `renderIsosurface(isosurface_value=None, component_index=0, time_index=0, color=(1, 0, 0), opacity=1.0, show=True, title=None) -> vtk actor`: VTK isosurface rendering.

**Utility identity checks**

- `instantiateAnotherAble() -> object`: create another wrapper of the same family.
- `isImaginable() -> bool`, `isSITKImaginable() -> bool`: legacy type helpers.
- `isImaginableInTheSameSpace(image) -> bool`: strict geometry comparison.

**Algorithm notes**

- Geometry changes and transforms use SimpleITK resampling.
- `overlayReport()` builds contour overlays and summary panels with Matplotlib.
- `renderIsosurface()` delegates to VTK isosurface extraction.

## SITKImaginable

- `SITKImaginable(filename=None, image=None, verbose=False)`
- Thin subclass of `Imaginable`.
- Use it when you want an explicit scalar-image wrapper name without changing behavior.

## Roiable

**Constructor**

- `Roiable(filename=None, image=None, verbose=False, roivalue=None)`

**Purpose**

- Binary mask/ROI wrapper built on `Imaginable`.
- Defaults to nearest-neighbor interpolation and nearest-neighbor extrapolation.
- If `roivalue` is passed, the input is binarized to that label.

**Core methods**

- `getCenterOfGravityCoordinates() -> tuple`, `getCenterOfGravityIndex() -> tuple`: intensity-weighted center of gravity.
- `getCentroidCoordinates() -> tuple`, `getCentroidIndex() -> tuple`: geometric centroid.
- `dilateRadius(radius=2) -> self`, `erodeRadius(radius=2) -> self`: binary morphology.
- `removeSmallObj(voxel_threshold=50, connectivity=26) -> self`: remove small connected components.
- `removeHoles(voxel_threshold=50, connectivity=26) -> self`: fill small holes.
- `keepBiggestObj(connectivity=26) -> self`: keep only the largest component.
- `fillBinaryHoles() -> self`: fill all enclosed holes.
- `smoothBinary(radius=1, mode="closing") -> self`: binary opening/closing smoothing.

**Transform and segmentation methods**

- `warpROI(displacement_field, target_image=None) -> self`: warp by deformation field with label-safe interpolation.
- `applyTransformToROI(transform, target_image=None) -> self`: backward-compatible alias for `applyTransform`.
- `refineWatershed(image, height_map=None, erosion_iters=3, dilation_iters=5, min_voxels=50) -> self`: marker-based watershed refinement.
- `refineRegionGrowing(image, multiplier=2.5, neighborhood_radius=1, n_iterations=3, max_distance_mm=10.0, n_seeds=100, prob_map=None, prob_threshold=0.1, min_voxels=50) -> self`: confidence-connected region growing.
- `refineGeodesicActiveContour(image, propagation=0.5, curvature=0.3, advection=1.0, iterations=50, rms_tolerance=0.001, allow_shrink=False, speed_image=None) -> self`: level-set refinement.
- `expandByProbability(prob_map, threshold=0.25, max_layers=5, max_distance_mm=None, fill_holes=True, min_voxels=50) -> self`: probability-constrained expansion.
- `shrinkByProbability(prob_map, threshold=0.15, max_layers=5, min_preserve_fraction=0.5, fill_holes=True) -> self`: probability-constrained shrinkage.
- `fuseSTAPLE(segmentations, confidence_threshold=0.5, min_voxels=50) -> Roiable`: combine several ROIs with STAPLE.

**Metrics, surfaces, and morphology**

- `compareTo(other) -> dict`: overlap and surface metrics.
- `getSurfaceDistances(other) -> dict`: mean/RMS/HD95/surface-Dice metrics.
- `exportSurfaceWithError(other, out_dir="/tmp/surface_error", highlight_percentile=95) -> None`: export meshes with error scalars.
- `getEdgeAlignmentScore(image, band_width=2, sigma_mm=1.0) -> float`: ROI-to-edge agreement score.
- `getCompactnessScore() -> float`: compactness ratio.
- `getConnectedComponentCount() -> int`: connected-component count.
- `getMorphometrics() -> dict`: PCA extents, Feret diameter, inscribed thickness, local thickness stats.
- `getPrincipalExtents() -> dict`, `getMaxFeret() -> float`, `getMaxInscribedThickness() -> float`, `getThicknessStats() -> dict`: focused morphometry helpers.
- `getDistanceMap() -> Imaginable`: Euclidean distance transform.
- `getShell(width_mm=5.0) -> Roiable`: annular band around the ROI.
- `getSurfaceMask() -> Roiable`: 1-voxel-thick surface.
- `getSignedDistanceMap() -> Imaginable`: signed Maurer distance map.
- `getSurfaceDistanceMap(other) -> Imaginable`: per-surface-voxel distance-to-target image.
- `getChangeMaps(other) -> dict[str, Roiable]`: removed/added/changed regions.
- `splitByConnectedComponents(sort_by="size") -> list[Roiable]`: separate components.
- `splitLeftRight(axis=0) -> tuple[Roiable, Roiable]`: component-based left/right split.
- `exportSurfaceVTP(out_path, distance_to=None, scalar_name="dist_mm") -> None`: VTK surface export.
- `describe() -> dict`: ROI-oriented summary.

**Algorithm notes**

- Refinement methods use watershed, confidence-connected region growing, and geodesic active contour.
- Morphometrics rely on Euclidean distance transforms, PCA extents, and Feret-style diameter estimation.
- Surface export methods use marching cubes / VTK mesh pipelines.

## LabelMapable

**Constructor**

- `LabelMapable(filename=None, image=None, verbose=False)`

**Purpose**

- Multi-label segmentation wrapper with label-preserving transforms and per-label analysis.

**Attributes**

- `dfltInterpolator = sitk.sitkNearestNeighbor`
- `dfltuseNearestNeighborExtrapolator = True`

**Label statistics and extraction**

- `noNearestNeighborExtrapolator() -> self`: disable NN extrapolation when desired.
- `getLabels(exclude_background=True) -> list[int]`: discover label values.
- `getCenterOfGravityCoordinatesPerLabel() -> dict`, `getCentroidCoordinatesPerLabel() -> dict`: per-label centers.
- `getCenterOfGravityCoordinates() -> tuple`, `getCenterOfGravityIndex() -> tuple`: aggregate center across labels.
- `getCentroidCoordinates() -> tuple`, `getCentroidIndex() -> tuple`: aggregate centroid across labels.
- `extractLabel(label_value) -> Roiable`: convert one label into a binary ROI.
- `setLabel(label_value, roi) -> self`: overwrite one label from a binary mask.
- `describe() -> dict`: label distribution summary.

**Comparison, priors, and refinement**

- `compareToByLabel(other, labels=None) -> dict[int, dict]`: per-label overlap/surface comparison.
- `buildPriors(tau=0.8, blur_sigma_mm=0.6, classes=None) -> tuple[Imaginable, list[int]]`: build soft priors from signed distances.
- `combineBinaryMasks(mask_list, priority_order=None) -> LabelMapable`: merge binary ROIs into a label map.
- `refineLabel(label_value, image, method="watershed", **kwargs) -> self`: refine a single label through ROI workflows.
- `refineAllLabels(image, method="watershed", resolve_overlaps=True, **kwargs) -> self`: batch refinement.
- `resolveOverlaps() -> self`: assign overlapping voxels to the nearest centroid label.

**Registration and longitudinal analysis**

- `warpLabelMap(displacement_field, target_image=None) -> self`: warp a label map with preserved labels.
- `applyTransformToLabelMap(transform, target_image=None) -> self`: backward-compatible alias for `applyTransform`.
- `registerTo(other, method="rigid", roi_values=None, iterations=200) -> tuple[LabelMapable, sitk.Transform]`: register using union-mask distance maps.
- `getChangeMapsByLabel(other, labels=None) -> dict[int, dict[str, Roiable]]`: per-label added/removed/changed regions.

**Algorithm notes**

- Priors are derived from per-label signed distance maps.
- `registerTo()` performs rigid or affine alignment on union-mask distance maps, then resamples the label map with nearest-neighbor interpolation.

## LabelMapableROI

- `LabelMapableROI(filename=None, image=None, verbose=False, labelsvalues=None)`
- Legacy compatibility wrapper around a list of `Roiable` objects.
- `removeSmallObj(voxel_threshold=50, connectivity=26) -> self`, `removeHoles(voxel_threshold=50, connectivity=26) -> self`, `keepBiggestObj(connectivity=26) -> self`: apply binary cleanup independently to each label ROI.
- `mergeLabels() -> self`: rebuild a label map from the internal ROI list.
- Prefer `LabelMapable` for new code.

## Fieldable

- `Fieldable(filename=None, image=None, verbose=False)`
- Vector/displacement-field flavored `Imaginable`.
- Uses nearest-neighbor interpolation by default but disables nearest-neighbor extrapolation.

**Methods**

- `setImageFromNumpy(nparray, refimage=None, spacing=None, origin=None, direction=None) -> self`: load vector data from `(z, y, x, components)` order.
- `setImageFromNumpyZYX(...) -> self`: explicit alias.
- `setImageFromNumpyXYZ(...) -> self`: backward-compatible loader for `(x, y, z, components)` order.
- `describe() -> dict`: report vector component count and inherited image information.

## Vectorable

- `Vectorable(filename=None, image=None, verbose=False)`
- Explicit vector-field analysis class for displacement fields, velocity fields, and similar multi-component images.

**Attributes**

- Inherits `Imaginable` state plus vector-specific validation that the image has more than one component.

**Methods**

- `getNumberOfComponents() -> int`: number of vector channels.
- `getMagnitude() -> Imaginable`: scalar magnitude image.
- `getComponent(component) -> Imaginable`: extract one vector component.
- `setComponent(component, scalar_image) -> self`: replace one component.
- `scaleVector(factors=1.0) -> self`: scale the vector field component-wise or uniformly.
- `normalize(target_magnitude=1.0) -> self`: normalize vectors to a target magnitude.
- `applyGaussianSmoothing(sigma=1.0, use_spacing=True) -> self`: Gaussian smoothing.
- `getVectorStatistics() -> dict`, `getStatistics() -> dict`: magnitude mean/std/min/max.
- `getMeanVector() -> np.ndarray`: global mean vector.
- `describe() -> dict`: vector-field summary with magnitude stats.
- `getDuplicate() -> Vectorable`: copy of the vector field.
- `plotOverlay(overlay=None, alpha=0.5, title=None, component=None, slice_idx=None, **kwargs) -> VectorPlotter`: static viewer.
- `viewInteractive(overlays=None, orientation=2, slice_idx=None, component=None, title=None, figsize=(14, 10), cmap="gray") -> InteractiveViewer`: interactive viewer.

**Algorithm notes**

- Magnitude is computed with `SimpleITK.VectorMagnitudeImageFilter`.
- Smoothing uses recursive Gaussian filtering.

## TimeSeriesable

- `TimeSeriesable(filename=None, image=None, verbose=False)`
- 4D scalar time-series wrapper where the 4th axis is time.

**Attributes**

- Inherits `Imaginable` and validates that the underlying image is 4D.

**Methods**

- `getNumberOfFrames() -> int`: number of time frames.
- `getFrame(frame_index) -> Imaginable`: extract one 3D frame.
- `setFrame(frame_index, frame_image) -> self`: replace one frame while preserving 4D metadata.
- `getFrameRange(start, end) -> TimeSeriesable`: contiguous temporal subset.
- `getTemporalMean() -> Imaginable`: mean over the time axis.
- `getTemporalVariance() -> Imaginable`: variance over the time axis.
- `getTemporalStandardDeviation() -> Imaginable`: standard deviation over the time axis.
- `applyFilterToAllFrames(filter_name, **kwargs) -> self`: apply the same scalar filter to every frame.
- `transformAllFrames(transform, interpolator="linear") -> self`: warp every frame with the same transform.
- `extractPhase(phase_number) -> Imaginable`: semantic alias for `getFrame`.
- `describe() -> dict`: time-series summary.
- `getDuplicate() -> TimeSeriesable`: copy of the 4D series.
- `plotOverlay(overlay=None, alpha=0.5, title=None, frame=None, slice_idx=None, **kwargs) -> TimeSeriesPlotter`: static viewer.
- `viewInteractive(overlays=None, orientation=2, slice_idx=None, frame=None, title=None, figsize=(14, 10), cmap="gray") -> InteractiveViewer`: interactive viewer.

**Algorithm notes**

- Temporal statistics are plain reductions along the frame axis.
- Frame extraction uses SimpleITK slicing/extraction and preserves spatial metadata in the 3D output.

## RoiComparison

- `RoiComparison(ref=None, test=None)`
- Classic overlap/shape comparison helper around two `Roiable` objects.

**Attributes**

- `Reference`: reference ROI.
- `Test`: ROI being evaluated.
- `Overlap`: cached `LabelOverlapMeasuresImageFilter`.

**Methods**

- `getReference()`, `setReference(reference)`: access/update the reference ROI.
- `getTest()`, `setTest(test)`: access/update the test ROI.
- `setOverlapFilter(overlap)`, `resetOverlap()`, `getOverlapFilter()`: cache and reuse the overlap filter.
- `getJaccard() -> float`, `getDice() -> float`, `getMeanOverlap() -> float`: overlap coefficients.
- `getVolumeSimilarity() -> float`, `getSimilarity() -> float`, `getVolumeSimilarityFromOverlap() -> float`, `getVolmeSimilarity() -> float`: same volume-similarity metric with compatibility aliases.
- `getFalseNegativeError() -> float`, `getFalsePositiveError() -> float`, `getFalsePostiveError() -> float`: classification errors.
- `getHausdorff() -> float`, `getHahusdorf() -> float`: Hausdorff distance and typo alias.
- `getSimilarityIndex() -> float`: `SimilarityIndexImageFilter` metric.
- `getOverlappedVoxels() -> float`, `getNonOverlappedVoxels() -> float`: voxel counts.
- `getAllMetrics() -> dict`: aggregate report.
- `printAllMetrics(json=None) -> None`: print or optionally save metrics.

**Algorithm notes**

- Uses SimpleITK overlap, Hausdorff, and similarity-index filters.
- Automatically resamples the reference ROI if the two masks are not in the same physical space.

## PlotViewer and derived plotters

### PlotViewer

- `PlotViewer(image, overlay=None, title="Image Viewer", cmap="gray", cmap_overlay="hot", figsize=(10, 8))`
- Base Matplotlib viewer for 2D/3D scalar images.

**Attributes**

- `image`, `overlay`, `title`, `cmap`, `cmap_overlay`, `alpha`, `current_slice`, `fig`, `ax`, `slider`.

**Methods**

- `show(slice_idx=None, alpha=0.5) -> None`: display the current image and optional overlay.
- `saveFigure(output_path, dpi=150) -> None`: save the current figure.

### ScalarPlotter

- `ScalarPlotter(image, overlay=None, title="Scalar Image Viewer", **kwargs)`
- No extra public methods beyond `PlotViewer`; it is the scalar-specific specialization.

### VectorPlotter

- `VectorPlotter(vector_image, overlay=None, title="Vector Field Viewer", show_vectors=False, **kwargs)`
- `show(component=0, alpha=0.5, slice_idx=None) -> None`: render a selected vector component.

### TimeSeriesPlotter

- `TimeSeriesPlotter(timeseries_image, overlay=None, title="Time Series Viewer", **kwargs)`
- `show(frame=0, alpha=0.5, slice_idx=None) -> None`: render one time frame.

### GridPlotter

- `GridPlotter(figsize=(12, 10))`
- `show_grid(slices, rows=None, cols=None, titles=None, cmap="gray", vmin=None, vmax=None, overlays=None, cmap_overlay="hot", alpha_overlay=0.5) -> None`: draw many 2D slices in a grid.

**Algorithm notes**

- Plot classes normalize arrays before display and resample overlays to the primary image when needed.

## OverlayManager

- `OverlayManager()`
- Lightweight registry for named overlay layers.

**Attributes**

- `overlays`: mapping of overlay name to image/array metadata.
- `order`: rendering order.

**Methods**

- `add(name, image, cmap="hot", opacity=0.5, visible=True) -> None`: store a SimpleITK overlay.
- `add_array(name, array, cmap="hot", opacity=0.5, visible=True) -> None`: store a NumPy overlay.
- `set_opacity(name, opacity) -> None`, `get_opacity(name) -> float`: opacity control.
- `set_visible(name, visible) -> None`, `is_visible(name) -> bool`: visibility control.
- `get_all() -> dict`, `get_count() -> int`: inspection.
- `remove(name) -> None`, `clear() -> None`: deletion helpers.

## InteractiveViewer

- `InteractiveViewer(image, title="Interactive Viewer", figsize=(14, 10), cmap="gray")`
- Full interactive browser for scalar, vector, and 4D images.

**Attributes**

- `image`: primary image.
- `overlays`: `OverlayManager` instance.
- `current_orientation`: `0=axial`, `1=sagittal`, `2=coronal`.
- `current_slice`: active slice.
- `current_frame`: active frame for 4D images.
- `current_component`: active vector component.
- `fig`, `ax_main`, `sliders`, `checkboxes`, `radio_component`: UI state.

**Methods**

- `add_overlay(overlay, name=None, cmap="hot", opacity=0.5) -> None`: add one overlay layer.
- `add_overlays(overlays, names=None, cmap="hot") -> None`: add several overlays.
- `show() -> None`: launch the GUI.

**Algorithm notes**

- Uses Matplotlib widgets for slice, orientation, frame, component, and overlay controls.
- Extracts 2D slices from 3D/4D inputs on demand and caches repeated requests.

## Related function modules

- `pyable.deformations`: deformation-field creation, transform conversion, multi-step transform application.
- `pyable.segmentation`: watershed, region growing, geodesic active contour, STAPLE, overlap resolution.
- `pyable.metrics`: ROI overlap metrics, surface distances, morphometrics, priors, change maps.
- `pyable.meshable`: `sitk <-> vtk` conversions.
- `pyable.utils`: directory processing, overlay utilities, rigid transform helper, slice export.
