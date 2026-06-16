try:
    import imaginable as ima
except:
    import pyable.imaginable as ima

import base64
from io import BytesIO

import numpy as np

# Input: expects 3xN matrix of points
# Returns R,t
# R = 3x3 rotation matrix
# t = 3x1 column vector

def rigid_transform_3D(A, B):
    assert A.shape == B.shape

    num_rows, num_cols = A.shape
    if num_rows != 3:
        raise Exception(f"matrix A is not 3xN, it is {num_rows}x{num_cols}")

    num_rows, num_cols = B.shape
    if num_rows != 3:
        raise Exception(f"matrix B is not 3xN, it is {num_rows}x{num_cols}")

    # find mean column wise
    centroid_A = np.mean(A, axis=1)
    centroid_B = np.mean(B, axis=1)

    # ensure centroids are 3x1
    centroid_A = centroid_A.reshape(-1, 1)
    centroid_B = centroid_B.reshape(-1, 1)

    # subtract mean
    Am = A - centroid_A
    Bm = B - centroid_B

    H = Am @ np.transpose(Bm)

    # sanity check
    #if linalg.matrix_rank(H) < 3:
    #    raise ValueError("rank of H = {}, expecting 3".format(linalg.matrix_rank(H)))

    # find rotation
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # special reflection case
    if np.linalg.det(R) < 0:
        print("det(R) < R, reflection detected!, correcting for it ...")
        Vt[2,:] *= -1
        R = Vt.T @ U.T

    t = -R @ centroid_A + centroid_B

    return R, t

def getImaginableSlice(I,axis,index):
    if axis==0:
        slicer=I.getSliceNormalI
        # ratio=k[2]/k[1]
    elif axis==1:
        slicer=I.getSliceNormalJ
        # self.ratio=k[2]/k[0]
    elif axis==2:
        slicer=I.getSliceNormalK
        # self.ratio=k[1]/k[0]
    o=type(I)(image=slicer(index))
    return o

def getImaginableSliceNumpy(I,axis,index):
    return getImaginableSlice(I,axis,index).getImageAsNumpy()


from PIL import Image, ImageDraw, ImageFont

import matplotlib
def saveSliceToImage(I,axis,index,fn,spacing=None):
    if spacing:
        I.changeImageSpacing(spacing)
    f=getImaginableSlice(I,axis=axis,index=index)
    matplotlib.image.imsave(fn, f,vmin=0.2,vmax=1,cmap='jet'  )


import matplotlib.pyplot as plt

def isListLikeIndex(value):
    return (
        (isinstance(value, np.ndarray) and value.ndim > 0)
        or isinstance(value, (list, tuple, range))
    )


def asIntList(value, name):
    if isinstance(value, np.ndarray):
        values = value.ravel()
    else:
        values = np.asarray(list(value) if isListLikeIndex(value) else [value]).ravel()
    if values.size == 0:
        raise ValueError(f"{name} must contain at least one value.")
    return [int(v) for v in values]


def makeAxisIndexPairs(axis, index, slice_offsets=None, index_mode='auto'):
    axes = asIntList(axis, "axis")
    multi_axis = isListLikeIndex(axis)
    multi_index = isListLikeIndex(index)
    index_values = asIntList(index, "index")

    if index_mode not in ('auto', 'point', 'cartesian'):
        raise ValueError("index_mode must be 'auto', 'point', or 'cartesian'.")

    if index_mode == 'auto':
        if multi_axis and multi_index and min(axes) >= 0 and max(axes) < len(index_values):
            index_mode = 'point'
        else:
            index_mode = 'cartesian'

    offsets = None
    if slice_offsets is not None:
        offsets = asIntList(slice_offsets, "slice_offsets")

    if index_mode == 'point':
        if min(axes) < 0 or max(axes) >= len(index_values):
            raise ValueError("index point must contain one coordinate for each requested axis.")
        axis_centers = [(axis_value, index_values[axis_value]) for axis_value in axes]
        if offsets is None:
            pairs = axis_centers
        else:
            pairs = [
                (axis_value, center_index + offset)
                for axis_value, center_index in axis_centers
                for offset in offsets
            ]
    else:
        if offsets is not None:
            if multi_index:
                raise ValueError(
                    "slice_offsets with multiple index values is ambiguous; "
                    "pass an index point with index_mode='point' or a scalar center index."
                )
            indices = [index_values[0] + offset for offset in offsets]
        else:
            indices = index_values
        pairs = [(axis_value, index_value) for axis_value in axes for index_value in indices]

    return pairs, (multi_axis or multi_index or slice_offsets is not None)


def _rgba_tuple(value):
    if len(value) == 3:
        return tuple(value) + (255,)
    return tuple(value)


def _get_default_font(size):
    for font_name in ("arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(font_name, size=size)
        except OSError:
            pass
    return ImageFont.load_default()


def _text_bbox(text, font):
    probe = Image.new('RGBA', (1, 1))
    draw = ImageDraw.Draw(probe)
    try:
        return draw.textbbox((0, 0), text, font=font)
    except AttributeError:
        width, height = draw.textsize(text, font=font)
        return (0, 0, width, height)


def _add_tight_title_band(rgba_uint8, title, font_size=12, padding=2, title_color=(255, 255, 255, 255), background=(0, 0, 0, 255)):
    if title is None or title == "":
        return rgba_uint8

    title = str(title)
    image = Image.fromarray(rgba_uint8, mode='RGBA')
    font = _get_default_font(font_size)
    bbox = _text_bbox(title, font)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    pad = max(0, int(padding))
    band_h = text_h + 2 * pad
    out_w = max(image.width, text_w + 2 * pad)
    out_h = image.height + band_h

    canvas = Image.new('RGBA', (out_w, out_h), _rgba_tuple(background))
    draw = ImageDraw.Draw(canvas)
    text_x = (out_w - text_w) // 2 - bbox[0]
    text_y = pad - bbox[1]
    draw.text((text_x, text_y), title, font=font, fill=_rgba_tuple(title_color))
    canvas.alpha_composite(image, ((out_w - image.width) // 2, band_h))
    return np.asarray(canvas)


def _save_encode_or_return_rgba(rgba_uint8, as_base64=False, data_uri=False, save=None):
    if save or as_base64:
        pil_image = Image.fromarray(rgba_uint8, mode='RGBA')
        if save:
            pil_image.save(save)
        if as_base64:
            buffer = BytesIO()
            pil_image.save(buffer, format='PNG')
            encoded = base64.b64encode(buffer.getvalue()).decode('ascii')
            if data_uri:
                encoded = f"data:image/png;base64,{encoded}"
            return encoded

    return rgba_uint8


def overlayNumpyImageAndNumpyLabelmapToImage(image, labelmap, image_cmap='gray', labelmap_cmap='jet', alpha_value=0.5, image_vmin=None, image_vmax=None, labelmap_vmin=None, labelmap_vmax=None, as_base64=False, data_uri=False, save=None, origin='lower', title=None, title_font_size=12, title_padding=2, title_color=(255, 255, 255, 255), background=(0, 0, 0, 255)):
    """Return only the composited image+overlay raster.

    By default this returns an ``(H, W, 4)`` uint8 RGBA NumPy array. If
    ``as_base64`` is True, it returns a PNG-encoded base64 string instead.
    """
    image = np.asarray(image)
    labelmap = np.asarray(labelmap)
    if image.shape != labelmap.shape:
        raise ValueError("image and labelmap must have the same shape.")
    if origin not in ('lower', 'upper'):
        raise ValueError("origin must be 'lower' or 'upper'.")

    image_norm = plt.Normalize(image_vmin, image_vmax)
    labelmap_norm = plt.Normalize(labelmap_vmin, labelmap_vmax)

    image_rgba = plt.get_cmap(image_cmap)(image_norm(image))
    labelmap_rgba = plt.get_cmap(labelmap_cmap)(labelmap_norm(labelmap))

    alpha = np.where(labelmap == 0, 0.0, alpha_value).astype(float)
    alpha = np.clip(alpha, 0.0, 1.0)[..., np.newaxis]

    rgb = image_rgba[..., :3] * (1.0 - alpha) + labelmap_rgba[..., :3] * alpha
    rgba = np.concatenate([rgb, np.ones_like(alpha)], axis=-1)
    if origin == 'lower':
        rgba = np.flipud(rgba)

    rgba_uint8 = np.round(np.clip(rgba, 0.0, 1.0) * 255).astype(np.uint8)
    rgba_uint8 = _add_tight_title_band(
        rgba_uint8,
        title,
        font_size=title_font_size,
        padding=title_padding,
        title_color=title_color,
        background=background,
    )

    return _save_encode_or_return_rgba(
        rgba_uint8,
        as_base64=as_base64,
        data_uri=data_uri,
        save=save,
    )


def overlayNumpyImageAndNumpyLabelmapGridToImage(images, labelmaps, image_cmap='gray', labelmap_cmap='jet', alpha_value=0.5, image_vmin=None, image_vmax=None, labelmap_vmin=None, labelmap_vmax=None, as_base64=False, data_uri=False, save=None, origin='lower', title=None, titles=None, ncols=None, tile_gap=0, title_font_size=12, title_padding=2, title_color=(255, 255, 255, 255), background=(0, 0, 0, 255)):
    """Return a tight raster montage of image+overlay slices."""
    images = list(images)
    labelmaps = list(labelmaps)
    n_items = len(images)

    if n_items == 0:
        raise ValueError("At least one image/labelmap pair is required.")
    if n_items != len(labelmaps):
        raise ValueError("images and labelmaps must contain the same number of items.")

    if titles is None:
        panel_titles = [None] * n_items
    else:
        panel_titles = np.asarray(titles, dtype=object).ravel().tolist()
        if len(panel_titles) != n_items:
            raise ValueError("titles must contain one title per image/labelmap pair.")

    if ncols is None:
        ncols = int(np.ceil(np.sqrt(n_items)))
    ncols = int(ncols)
    if ncols < 1:
        raise ValueError("ncols must be at least 1.")
    ncols = min(ncols, n_items)
    nrows = int(np.ceil(n_items / ncols))
    gap = max(0, int(tile_gap))

    tiles = []
    for image, labelmap, panel_title in zip(images, labelmaps, panel_titles):
        tiles.append(
            overlayNumpyImageAndNumpyLabelmapToImage(
                image,
                labelmap,
                image_cmap=image_cmap,
                labelmap_cmap=labelmap_cmap,
                alpha_value=alpha_value,
                image_vmin=image_vmin,
                image_vmax=image_vmax,
                labelmap_vmin=labelmap_vmin,
                labelmap_vmax=labelmap_vmax,
                as_base64=False,
                data_uri=False,
                save=None,
                origin=origin,
                title=panel_title,
                title_font_size=title_font_size,
                title_padding=title_padding,
                title_color=title_color,
                background=background,
            )
        )

    tile_h = max(tile.shape[0] for tile in tiles)
    tile_w = max(tile.shape[1] for tile in tiles)
    out_w = ncols * tile_w + (ncols - 1) * gap
    out_h = nrows * tile_h + (nrows - 1) * gap

    canvas = Image.new('RGBA', (out_w, out_h), _rgba_tuple(background))
    for i, tile in enumerate(tiles):
        row = i // ncols
        col = i % ncols
        tile_image = Image.fromarray(tile, mode='RGBA')
        x = col * (tile_w + gap) + (tile_w - tile_image.width) // 2
        y = row * (tile_h + gap) + (tile_h - tile_image.height) // 2
        canvas.alpha_composite(tile_image, (x, y))

    rgba_uint8 = np.asarray(canvas)
    rgba_uint8 = _add_tight_title_band(
        rgba_uint8,
        title,
        font_size=title_font_size,
        padding=title_padding,
        title_color=title_color,
        background=background,
    )

    return _save_encode_or_return_rgba(
        rgba_uint8,
        as_base64=as_base64,
        data_uri=data_uri,
        save=save,
    )

def overlayNumpyImageAndNumpyLabelmap(image, labelmap, image_cmap='gray', labelmap_cmap='jet', alpha_value=0.5, image_vmin=None, image_vmax=None, labelmap_vmin=None, labelmap_vmax=None,show=False,save=None,title=None,labelmap_name=None, ax=None, colorbar=True):
    if ax is None:
        ax = plt.gca()

    # Display the image as it is
    im_handle = ax.imshow(image, cmap=image_cmap, vmin=image_vmin, vmax=image_vmax,origin='lower')

    # Clip the labelmap values to the desired range
    # labelmap_clipped = labelmap
    # if (labelmap_vmin is None) or (labelmap_vmax is None):
    #     if labelmap_vmin is None:
    #         labelmap_vmin = np.min(labelmap)
    #     if labelmap_vmax is None:
    #         labelmap_vmax = np.max(labelmap)
            
    #     labelmap_clipped = np.clip(labelmap, labelmap_vmin, labelmap_vmax)

    # Normalize the labelmap
    labelmap_norm = plt.Normalize(labelmap_vmin, labelmap_vmax)

    # Apply colormap to the normalized labelmap
    labelmap_colored = plt.get_cmap(labelmap_cmap)(labelmap_norm(labelmap))


    # Create an alpha channel based on the labelmap
    alpha_channel = np.where(labelmap == 0, 0, alpha_value)

    # Replace the alpha channel in the colored labelmap
    labelmap_colored[..., 3] = alpha_channel

    # Overlay the labelmap on top of the image
    lbl=ax.imshow(labelmap_colored,origin='lower')
    # Create a ScalarMappable object for the colorbar
    sm = plt.cm.ScalarMappable(cmap=labelmap_cmap, norm=labelmap_norm)
    sm.set_array([])

    # Add the colorbar
    colorbar_handle = None
    if colorbar:
        colorbar_handle = ax.figure.colorbar(sm, ax=ax, label=labelmap_name)
    

    if title:
        ax.set_title(title)
    
    if save:
        ax.figure.savefig(save,dpi=300)
    if show:
        plt.show()

    return {
        'figure': ax.figure,
        'axis': ax,
        'image': im_handle,
        'overlay': lbl,
        'colorbar': colorbar_handle,
    }


def overlayNumpyImageAndNumpyLabelmapGrid(images, labelmaps, image_cmap='gray', labelmap_cmap='jet', alpha_value=0.5, image_vmin=None, image_vmax=None, labelmap_vmin=None, labelmap_vmax=None,show=False,save=None,title=None,titles=None,labelmap_name=None, colorbar=False, figsize=None):
    images = list(images)
    labelmaps = list(labelmaps)
    n_items = len(images)

    if n_items == 0:
        raise ValueError("At least one image/labelmap pair is required.")
    if n_items != len(labelmaps):
        raise ValueError("images and labelmaps must contain the same number of items.")

    if titles is None:
        panel_titles = [None] * n_items
    else:
        panel_titles = np.asarray(titles, dtype=object).ravel().tolist()
        if len(panel_titles) != n_items:
            raise ValueError("titles must contain one title per image/labelmap pair.")

    n_grid = int(np.ceil(np.sqrt(n_items)))
    if figsize is None:
        figsize = (3.0 * n_grid, 3.0 * n_grid)

    fig, axes = plt.subplots(n_grid, n_grid, figsize=figsize, squeeze=False)
    axes_flat = axes.ravel()
    panels = []

    for i, (image, labelmap) in enumerate(zip(images, labelmaps)):
        ax = axes_flat[i]
        panel = overlayNumpyImageAndNumpyLabelmap(
            image,
            labelmap,
            image_cmap=image_cmap,
            labelmap_cmap=labelmap_cmap,
            alpha_value=alpha_value,
            image_vmin=image_vmin,
            image_vmax=image_vmax,
            labelmap_vmin=labelmap_vmin,
            labelmap_vmax=labelmap_vmax,
            show=False,
            save=None,
            title=None,
            labelmap_name=labelmap_name,
            ax=ax,
            colorbar=colorbar,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        if panel_titles[i]:
            ax.set_title(str(panel_titles[i]), fontsize=9, pad=2)
        panels.append(panel)

    for ax in axes_flat[n_items:]:
        ax.axis('off')

    if title:
        fig.suptitle(title)

    has_panel_titles = any(panel_titles)
    top = 0.92 if title else (0.98 if has_panel_titles else 1.0)
    hspace = 0.16 if has_panel_titles else 0.02
    fig.subplots_adjust(left=0, right=1, bottom=0, top=top, wspace=0.02, hspace=hspace)

    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight', pad_inches=0.02)
    if show:
        plt.show()

    return {
        'figure': fig,
        'axes': axes,
        'panels': panels,
    }
    
    



if __name__=="__main__":
    IM=ima.Imaginable('/data/MYDATA/fulldixon-images/C-1/data/wo.nii')
    R=ima.LabelMapable('/data/MYDATA/fulldixon-images/C-1/data/fo.nii')
    P=ima.LabelMapable('/data/MYDATA/fulldixon-images/C-1/data/roi.nii.gz')
    ORIENTATION='RPI'
    # IM.dicomOrient(ORIENTATION)
    # R.resampleOnTargetImage(IM)
    # R.dicomOrient(ORIENTATION)
    # P.resampleOnTargetImage(IM)
    # P.dicomOrient(ORIENTATION)
    # for SL in range(IM.getImageSize(0)):

    #     im=getImaginableSliceNumpy(IM,0,SL)
    #     im2=getImaginableSliceNumpy(R,0,SL)
    #     r=getImaginableSliceNumpy(P,0,SL)
    #     im2[r==0]=0
    #     overlayNumpyImageAndNumpyLabelmap(im.T,im2.T)
    #     plt.savefig(f'/g/{SL}.png')
    #     plt.close()
    
    
    NP=P.getImageAsNumpy()
    NP[NP>0]=1
    P.setImageFromNumpy(NP)
    
    R.multiply(P)
    plt.subplot(1,2,1)
    IM.overlayAble(R,1,160,show=False,title='Overlayed',labelmap_name='Femur',alpha_value=0.3)
    plt.subplot(1,2,2)
    IM.overlayAble(R,2,40,show=False,title='Overlayed',labelmap_name='Femur',alpha_value=0.3)
    plt.show()
    
    
    

# if __name__=="__main__":

#     # P=pn.Pathable('/g/a/f.nii.gz')
#     # for p in P.getFilesInPathByExtension():
#     #     I=Imaginable(filename=p)
#     #     l=pn.Pathable(p)
#     #     l.changeExtension('png')
#     #     l.appendPath('f')
#     #     l.ensureDirectoryExistence()
#     #     for c,v in zip(range(3),[40,80,30]):
#     #         l.addPrefix(f'{c}_')
#     #         if not l.exists():
#     #             saveSliceToImage(I,axis=c,index=v,fn=l.getPosition())
#     #         l.undo()
        


   

#     # NEW=np.array([-0.9993099234645667, -6.385608658999428e-08, -0.037144005107034125, 1.898319864628949e-05, -0.9999998702797991, -0.0005089990447097432, -0.03714399688886282, -0.0005093529042066276, 0.9993097937099291]).reshape((3,3))
#     # ORIGINAL=np.eye(3)

#     # n=10
#     # t=np.array([[0],[0],[0]])
#     # A = np.random.rand(3, n)
#     # B = NEW@A + t

    #     # ret_R, ret_t = rigid_transform_3D(A, B)
    #     # print(ret_R)
    #     # print(ret_t)


import glob
import pandas as pd
import os

def processImageDirectory(directory, processor_func, file_pattern='*.nii.gz', 
                         output_csv=None, recursive=True, verbose=False):
    """
    Batch process images in a directory with a processor function.
    
    Walks through directory, applies processor_func to each matching image,
    aggregates results to pandas DataFrame and optionally exports to CSV.
    
    Args:
        directory (str): Root directory to search
        processor_func (callable): Function that takes (Imaginable) -> dict
                                  The dict should contain results to aggregate
        file_pattern (str): Glob pattern for files. Default: '*.nii.gz'
                           Examples: '*.nii.gz', '*.nii', '*.img'
        output_csv (str, optional): Path to save results as CSV. If None, no export.
        recursive (bool): Search recursively in subdirectories. Default: True
        verbose (bool): Print debug information. Default: False
        
    Returns:
        pd.DataFrame: Results with one row per processed image
                     'filepath' column contains the file path
                     Other columns from processor_func output
                     
    Example:
        >>> def process_img(img):
        ...     # img is an Imaginable object
        ...     return {
        ...         'size': img.getImageSize(),
        ...         'spacing': img.getImageSpacing(),
        ...         'num_pixels': np.prod(img.getImageSize())
        ...     }
        >>> 
        >>> results = processImageDirectory(
        ...     '/data/images', 
        ...     process_img,
        ...     file_pattern='*.nii.gz',
        ...     output_csv='/output/results.csv',
        ...     verbose=True
        ... )
        >>> print(results.head())
    """
    if not os.path.isdir(directory):
        raise ValueError(f"Directory does not exist: {directory}")
    
    # Find all matching files
    if recursive:
        search_pattern = os.path.join(directory, '**', file_pattern)
        files = glob.glob(search_pattern, recursive=True)
    else:
        search_pattern = os.path.join(directory, file_pattern)
        files = glob.glob(search_pattern, recursive=False)
    
    if len(files) == 0:
        if verbose:
            print(f"No files matching '{file_pattern}' found in {directory}")
        return pd.DataFrame()
    
    if verbose:
        print(f"Found {len(files)} files matching '{file_pattern}'")
    
    # Process each file
    results = []
    for i, filepath in enumerate(files):
        try:
            if verbose:
                print(f"[{i+1}/{len(files)}] Processing: {filepath}")
            
            # Load image as Imaginable
            img = ima.Imaginable(filepath)
            
            # Apply processor function
            result = processor_func(img)
            
            # Add filepath to result
            if isinstance(result, dict):
                result['filepath'] = filepath
                results.append(result)
            else:
                if verbose:
                    print(f"  Warning: processor returned non-dict type {type(result)}")
                continue
                
        except Exception as e:
            if verbose:
                print(f"  Error processing {filepath}: {e}")
            continue
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    if verbose:
        print(f"Successfully processed {len(df)} files")
        print(f"Columns: {list(df.columns)}")
    
    # Export to CSV if requested
    if output_csv:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
        df.to_csv(output_csv, index=False)
        if verbose:
            print(f"Results saved to: {output_csv}")
    
    return df
