import numbers
from pynico import pynico as pn
import SimpleITK as sitk
import numpy as np
import copy
import matplotlib.pyplot as plt
import os

try:
    from .utils import *
except:
    try:
        from utils import *
    except:
        from pyable.utils import *
        
from skimage import data, filters, measure, morphology
import itertools
def get_image_corners_nd(image):
    """
    Get the corner pixel values of an n-dimensional SimpleITK image.

    Parameters:
    - image: A SimpleITK image of any dimension.

    Returns:
    - A dictionary with the corner coordinates as keys and pixel values as values.
    """
    # Ensure the image is not null
    if image is None:
        raise ValueError("The input image is None.")
    
    # Get the size of the image (size per dimension)
    size = image.GetSize()
    
    # Generate all corner coordinates by creating combinations of (0 or max index) for each dimension
    corners = list(itertools.product(*[(0, s-1) for s in size]))
  
    return corners

def get_image_corners_coordinates(image):
    """
    Get the corner coordinates of an n-dimensional SimpleITK image.

    Parameters:
    - image: A SimpleITK image of any dimension.

    Returns:
    - A list of corner coordinates.
    """
    # Ensure the image is not null
    if image is None:
        raise ValueError("The input image is None.")
    corners=[]
    for c in get_image_corners_nd(image):
       corners.append(image.TransformIndexToPhysicalPoint(c))
    return corners

def dcm2niixFieldsToJson(fn,field_list=["RepetitionTime","FlipAngle","MagneticFieldStrength","ScanningSequence","NonlinearGradientCorrection","SliceThickness","SpacingBetweenSlices","SAR","EchoTime","RepetitionTime","SpoilingState","FlipAngle","PartialFourier","TxRefAmp","PixelBandwidth","PatientPosition","MRAcquisitionType","ImagingFrequency","ScanOptions"]):
    """
    This function convert the fields of dcm2niix file to a json file with the fields specified in field_list

    Args:
        fn (str): filename
        field_list (list, optional): _description_. Defaults to ["RepetitionTime","FlipAngle","MagneticFieldStrength","ScanningSequence","NonlinearGradientCorrection","SliceThickness","SpacingBetweenSlices","SAR","EchoTime","RepetitionTime","SpoilingState","FlipAngle","PartialFourier","TxRefAmp","PixelBandwidth","PatientPosition","MRAcquisitionType","ImagingFrequency","ScanOptions"].
    """    
    import json
    with open(fn) as f:
        data = json.load(f)
    o={}
    for f in field_list:
        try:
            o[f]=data[f]
        except:
            o[f]=None

    return o
    
def create_affine_matrix(rotation=[0,0,0], scaling=[1,1,1]):
    """Creates a 3D affine matrix given three rotations, three scalings, and three translations.

    Args:
        rotation: The angle of rotations 3d in degrees
        scaling: The scaling factor along the x-axis.

    Returns:
        A 3x3 NumPy array representing the affine matrix.
    """
    rx,ry,rz=[np.deg2rad(r) for r in rotation]
    sx,sy,sz=scaling
    rotation_matrix= np.array([[np.cos(rx), -np.sin(rx), 0],
                                [np.sin(rx), np.cos(rx), 0],
                                [0, 0, 1]]) @ np.array([[np.cos(ry), 0, np.sin(ry)],
                                                    [0, 1, 0],
                                                    [-np.sin(ry), 0, np.cos(ry)]]) @ np.array([[np.cos(rz), np.sin(rz), 0],
                                                                                                    [-np.sin(rz), np.cos(rz), 0],
                                                                                                    [0, 0, 1]])
    
    scaling_matrix = np.array([[sx, 0, 0],
                              [0, sy, 0],
                              [0, 0, sz]])


    return rotation_matrix @ scaling_matrix 

def getMaskedNunmpyArray(IM,ROI):
    """Return the values inside a region of interest
    Args:
      IM:
        Imaginable
      ROI:
        ROIable

    Returns:
      Numpy arrays aof values
    """
    o=IM.getImageAsNumpy()
    m=ROI.getImageAsNumpy()
    return o[np.where(m>0)]

def transform_point(P,transform):
    return transform.TransformPoint(P)


def getMatrixToPatientOrientationHF(L):
    #input is an imaginable
    
    transform=sitk.AffineTransform(L.getImageDimension())
    transform.SetMatrix(L.getImageDirection())
    O0= transform_point([0]*L.getImageDimension(),transform)
    O1= transform_point(L.getImageSize(),transform)
    
    
    L=[['Left','Right'],['Anterior','Posterior'],['Down','Up']]
    if O0[0]<O1[0]:
        L[0].reverse()   
    if O0[1]>O1[1]:
        L[1].reverse()   
    if O0[2]>O1[2]:
        L[2].reverse()   
    return L

def getMatrixToPatientOrientation(IM):
    """_summary_
    The function help you understand the real world directions and the matrix one
    an output that goes like:
    [['Left','Right'],['Anterior','Posterior'],['Down','Up']]
    means that the first direction goes from left to right, the second from Anterior to posterior and so on
    Args:
        IM (able): input image

    Returns:
        _type_: _description_
    """    
    #input is an imaginable  
    L=[['Right','Left'],['Anterior','Posterior'],['Down','Up']]
    O0= IM.getCoordinatesFromIndex([0]*IM.getImageDimension())
    O1= IM.getCoordinatesFromIndex(IM.getImageSize())    
    for a in range(IM.getImageDimension()):
        if O0[a]<O1[a]:
            L[a].reverse()   
    return L


def numpyToImaginable(x,ref=None,vector=False):
    T=Imaginable()
    T.setImageFromNumpy(x,refimage=ref,vector=vector)
    return T


def saveNumpy(x,fn,ref=None):
    T=numpyToImaginable(x,ref)
    T.writeImageAs(fn)
def getTransformFromFile(h):
    if not pn.isCollection(h):
        h=[h]
    o=[]
    for t in h:
        if isinstance(h,pn.Pathable):
            st=h.getPosition()
        else:
            st=t
        o.append(sitk.ReadTransform(st))
    return o    

def getSITKMetaInfo(image):
    o={}
    for key in image.GetMetaDataKeys():
            o[key]=image.GetMetaData(key)
    return o
def setSITKMetaInfo(image,m):
    for key,value in m.items():
            image.SetMetaData(key,value)
    return image
        
class IndexViewer(object):
        def __init__(self, ax, Ima,km=[True,True],normal=2,ind=None):
            self.ax = ax

            # ax.set_title('use scroll wheel to navigate the image size(' + str(X.shape) +')')
            self.slices=Ima.getImageSize(normal)
            self.ind=ind
            self.Ima=Ima
            k=list(Ima.getImageSpacing())
            if self.ind==None:
                self.ind = self.slices//2
                if normal==0:
                    self.slicer=self.Ima.getSliceNormalI
                    self.ratio=k[2]/k[1]
                elif normal==1:
                    self.slicer=self.Ima.getSliceNormalJ
                    self.ratio=k[2]/k[0]
                elif normal==2:
                    self.slicer=self.Ima.getSliceNormalK
                    self.ratio=k[1]/k[0]
            self.UD=not km[0]
            self.LR=not km[1]
            self.X=self.getTheX()
            self.im = ax.imshow(self.X)
            ax.set_aspect(self.ratio)
            self.update()
        def getTheX(self):
            o=SITKImaginable(image=self.slicer(self.ind)).getImageAsNumpyZYX()
            if self.UD:
                o=np.flipud(o)
            if self.LR:
                o=np.fliplr(o)
            return o

        def onscroll(self, event):
            if event.button == 'up':
                self.ind = (self.ind + 1) % self.slices
            else:
                self.ind = (self.ind - 1) % self.slices
            self.update()

        def update(self):
            self.X=self.getTheX()
            self.im.set_data(self.X)
            self.ax.set_ylabel('slice ' + str(self.ind)  + "/" +str(self.slices))
            self.im.axes.figure.canvas.draw()       
        def onclick(self,event):           
            if event.button == 2:
               print('you pressed', event.button, event.xdata, event.ydata)
               print('value %s' % self.X[int(np.floor(event.ydata)),int(np.floor(event.xdata))])



def getmeTheSimpleITKImage(x):
    if(issubclass(type(x),sitk.Image)):
        return x
    elif (issubclass(type(x),Imaginable)):
        return x.getImage()
    elif (isinstance(x,str)):
        return Imaginable(filename=x).getImage()
    else:
        raise Exception("I don't know this image type!!! what should i do?? ask Eros eros.montin@gmail.com")
def copythethreeinfosonandsetthemtimage(source,reference):
    sp,o,d=getSITKImageInfo(reference)
    return setSITKImageInfo(source,sp,o,d)


def setSITKImageInforFromImage(nda,ima):
    # ORSP=nda.GetSpacing()
    REF=getmeTheSimpleITKImage(ima)
    nda=copythethreeinfosonandsetthemtimage(nda,REF)
    #if the resolution is identical no prob otherwise we need to fix the poition of the origin
    # if not np.array_equiv(ORSP,REF.GetSpacing()):
    #     nda=fixOriginOnResamping(nda,ima)
    return nda



def setSITKImageInfo(nda,spacing,origin,direction):
    if spacing:
        osp=nda.GetSpacing()
        Onda=copy.deepcopy(nda)
        nda.SetSpacing(spacing)
    if direction:
        nda.SetDirection(direction)
    if origin:
        nda.SetOrigin(origin)
        # if np.array_equiv(spacing,osp):
        #     nda=fixOriginOnResamping(nda,Onda)
    return nda
def getSITKImageInfo(nda):
    return nda.GetSpacing(),nda.GetOrigin(), nda.GetDirection()




class Imaginable:
    """
    Base mutable wrapper around a ``SimpleITK.Image``.

    The class stores image history in an internal stack so geometry edits,
    arithmetic, transforms, and filtering can be chained and undone.
    """
    def __init__(self,filename=None,image=None,verbose=False):
        self.verbose=verbose
        self.dfltInterpolator=sitk.sitkLinear
        self.dfltuseNearestNeighborExtrapolator=False
        self.imageStack =pn.Stack()    
        self.log=pn.Log()
        self.settings={
            "spacingMinSize":2 #max number of digit 
            }

        if ((filename) and (image)):
            if(pn.Pathable(filename).exists()):
                raise Exception("you can't put an existing filename and an image")
        if filename:
            self.setInputFileName(filename)
        if image:
            self.setImage(getmeTheSimpleITKImage(image),'image set as simpeitk image at class initialization')
            if not self.getInputFileName():
                self.InputFileName=pn.createRandomTemporaryPathableFromFileName('a.nii.gz').getPosition()


        
    def getValuesInRoi(self,roi):
        return getMaskedNunmpyArray(self,roi)
    
    def getImage(self):
        return self.imageStack.peek()

    def setImage(self,p,w=None):
        self.imageStack.push(p)
        if w:
            self.__tellme__(w)
        return self
    
    def reset(self):
        while self.imageStack.size()>1:
            self.undo()
    def isImaginableInTheSameSpace(self,image):
        return ((self.getImageSize() == image.getImageSize()) and
             (self.getImageDirection() == image.getImageDirection()) and
             (self.getImageOrigin() == image.getImageOrigin()) and
               (self.getImageSpacing() == image.getImageSpacing() ))
    
    def undo(self):
        if self.imageStack.size()>1:
            self.imageStack.pop()
            self.__tellme__('image popped')
        return self

    def isImageSet(self):
        try:
            if self.getImage() is None:
                return False
            else:
                return True
        except:
            return False

    def setVerbose(self,v):
        self.verbose=v
    def getVerbose(self):
        return self.verbose

    def getInputFileName(self):
        try:
            return self.InputFileName
        except:
            return None
    def setInputFileName(self,fn):
        self.InputFileName = fn
        
        if pn.Pathable(fn).exists():
            self.setImage(self.__readImage__(fn),'image read from filename')
            

    
    def __readImage__(self,f=None):
        if not f:
            f=self.getInputFileName()
        return sitk.ReadImage(f)
    
    def writeImageAs(self,filename,force=True):
        if ((not pn.Pathable(filename).exists()) or (force)):
            if os.path.dirname(filename)!='':
                pn.Pathable(filename).ensureDirectoryExistence()
            try:
                sitk.WriteImage(self.getImage(), filename)
                return filename
            except:
                raise Exception(f" file {filename} can't be written!!")

    # Convenience alias
    write = writeImageAs
    
    def printImageInfo(self):
        image= self.getImage()
        o={}
        for key in image.GetMetaDataKeys():
            print("\"{0}\":\"{1}\"".format(key, image.GetMetaData(key)))
            o[key]=image.GetMetaData(key)
        return o
    def getImageAsNumpy(self):
        """
        Returns the image as a numpy array in standard (Z, Y, X) ordering.
        
        This is the standard convention for numpy arrays and PyTorch tensors in medical imaging.
        For 3D images: (Z, Y, X) = (depth/slices, height/rows, width/cols)
        For 2D images: (Y, X) = (height/rows, width/cols)
        
        Returns:
            numpy.ndarray: Image array in (Z, Y, X) order for 3D, (Y, X) for 2D
        
        Note: This changed in v3! Previously returned (X, Y, Z). 
              Use getImageAsNumpyXYZ() if you need the old behavior.
        """
        image = self.getImage() 
        return sitk.GetArrayFromImage(image)
    
    def getImageAsNumpyZYX(self):
        """
        Returns the image as a numpy array in (Z, Y, X) ordering.
        Alias for getImageAsNumpy() for explicit clarity.
        
        Returns:
            numpy.ndarray: Image array in (Z, Y, X) order
        """
        return self.getImageAsNumpy()
    
    def getImageAsNumpyXYZ(self):
        """
        Returns the image as a numpy array in (X, Y, Z) ordering.
        
        DEPRECATED: This is non-standard for numpy. Provided for backward compatibility only.
        The old getImageAsNumpy() returned this ordering in v2.
        
        Returns:
            numpy.ndarray: Image array in (X, Y, Z) order - NON-STANDARD
        """
        # Transpose from (Z,Y,X) to (X,Y,Z)
        L = list(range(self.getImageDimension()))
        L.reverse()
        o = np.transpose(self.getImageAsNumpy(), L)
        return o
    
    def getImageAsNumpyForPyTorch(self):
        """
        Returns the image as a numpy array in PyTorch-compatible format.
        Alias for getImageAsNumpy() since v3 uses standard (Z,Y,X) ordering.
        
        For 3D: returns (D, H, W) = (Z, Y, X)
        For 2D: returns (H, W) = (Y, X)
        
        Returns:
            numpy.ndarray: Image array ready for PyTorch tensors
        
        Example:
            >>> arr = img.getImageAsNumpyForPyTorch()
            >>> tensor = torch.from_numpy(arr).float()
            >>> # Shape is (D, H, W) for 3D or (H, W) for 2D
        """
        return self.getImageAsNumpy()
    
    def getITKImage(self):
        """
        Returns the underlying SimpleITK Image object.
        Explicit alias for getImage() for clarity.
        
        The ITK image uses (X, Y, Z) indexing and physical coordinates in millimeters.
        
        Returns:
            SimpleITK.Image: The SimpleITK image object
        """
        return self.getImage()
    
    def getVTKImage(self):
        """
        Convert and return the image as a VTK object.
        Useful for 3D rendering and visualization with VTK-based tools.
        
        Returns:
            vtk.vtkImageData: The VTK image object
        """
        from .meshable import sitk2vtk
        return sitk2vtk(self.getImage())
    
    
    
    ###
    def getParaViewData(
        self,
        space="lps",
        stride=1,
        array_name="values",
    ):
        """
        Return the image as a coordinate-aware vtkStructuredGrid.
        """
        from .meshable import sitk_to_structured_grid

        return sitk_to_structured_grid(
            self.getImage(),
            space=space,
            stride=stride,
            array_name=array_name,
        )


    def writeParaView(
        self,
        output_path,
        space="lps",
        stride=1,
        array_name="values",
    ):
        """
        Write the image as a ParaView-compatible VTS file.
        """
        from .meshable import write_vtk_dataset

        dataset = self.getParaViewData(
            space=space,
            stride=stride,
            array_name=array_name,
        )

        return write_vtk_dataset(dataset, output_path)
    ###
    
    def overlayReport(
        self,
        overlay,
        spacing=None,
        orientation='LPS',
        views='all',
        slice_offsets=None,
        fill_alpha=0.20,
        contour_alpha=0.95,
        overlay_color=(1.0, 0.0, 0.0),
        contour_iterations=2,
        image_cmap='gray',
        figsize=None,
        title=None,
        show=False,
        save=None,
        dpi=150,
        stats=True,
    ):
        """Generate a publication-ready overlay report figure.

        Creates a multi-view (axial/coronal/sagittal) overlay of a
        Roiable, LabelMapable, or any Imaginable on top of *self*.
        Both images are oriented to a canonical orientation and optionally
        resampled to isotropic resolution before slicing.

        Parameters
        ----------
        overlay : Imaginable | Roiable | LabelMapable
            The overlay image.  Binary masks get a filled region + bold
            contour; multi-label maps use the provided ``overlay_color``
            for non-zero voxels.
        spacing : list[float] | float | None
            Target isotropic spacing in mm (e.g. ``1.5`` or ``[1.5, 1.5, 1.5]``).
            ``None`` keeps the original spacing.
        orientation : str
            Three-letter DICOM orientation code applied to **both**
            images before slicing (default ``'LPS'``).
        views : str | list[str]
            ``'all'`` (default) → ``['axial', 'coronal', 'sagittal']``,
            or a subset list, or a single view name.
        slice_offsets : list[int] | None
            Offsets (in voxels) from the overlay center-of-mass.
            ``None`` → ``[0]`` (center slice only).
        fill_alpha : float
            Opacity of the filled overlay region (0–1).
        contour_alpha : float
            Opacity of the bold boundary contour (0–1).
        overlay_color : tuple[float, float, float]
            RGB colour for overlay fill and contour (0–1 range).
        contour_iterations : int
            Dilation iterations for contour thickness (0 = no contour).
        image_cmap : str
            Matplotlib colormap name for the background image.
        figsize : tuple | None
            ``(width, height)`` in inches.  ``None`` auto-computes.
        title : str | None
            Figure suptitle.  ``None`` auto-generates.
        show : bool
            Call ``plt.show()`` after rendering.
        save : str | None
            File path to save the figure (PNG, PDF, …).
        dpi : int
            DPI used when saving.
        stats : bool
            If ``True`` add an extra column with volume / mean / std
            statistics for the overlay region.

        Returns
        -------
        dict
            ``{'figure': fig, 'axes': axes_array, 'stats': stats_dict}``

        Examples
        --------
        >>> from pyable import Imaginable, Roiable
        >>> img = Imaginable('t1.nii.gz')
        >>> roi = Roiable('tumor.nii.gz')
        >>> img.overlayReport(roi, spacing=1.5, save='report.png')

        >>> # specific views, custom colour
        >>> img.overlayReport(roi, views=['axial','coronal'],
        ...                   overlay_color=(0, 1, 0), fill_alpha=0.3)

        >>> # LabelMapable overlay
        >>> from pyable import LabelMapable
        >>> seg = LabelMapable('seg.nii.gz')
        >>> img.overlayReport(seg, spacing=[1,1,1], orientation='RAS')
        """
        import copy
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.cm as mcm
        from scipy import ndimage

        # ── Resolve parameters ──────────────────────────────────────
        if isinstance(spacing, (int, float)):
            spacing = [float(spacing)] * 3
        if isinstance(views, str):
            views = ['axial', 'coronal', 'sagittal'] if views == 'all' else [views]
        views = [v.lower() for v in views]
        if slice_offsets is None:
            slice_offsets = [0]
        overlay_color = np.array(overlay_color, dtype=np.float64)

        view_axis_map = {'axial': 0, 'coronal': 1, 'sagittal': 2}

        # ── Prepare working copies ──────────────────────────────────
        img_work = copy.deepcopy(self)
        ovl_work = copy.deepcopy(overlay)

        # Orient both to the same canonical orientation
        img_work.dicomOrient(orientation)
        ovl_work.dicomOrient(orientation)

        # Resample to target spacing
        if spacing is not None:
            img_work.changeImageSpacing(spacing)
            ovl_work.changeImageSpacing(spacing)

        # Resample overlay onto image grid for guaranteed shape match
        ovl_sitk = sitk.Resample(
            ovl_work.getImage(),
            img_work.getImage(),
            sitk.Transform(),
            sitk.sitkNearestNeighbor,
            0.0,
            ovl_work.getImage().GetPixelID(),
        )
        ovl_work.setImage(ovl_sitk)

        img_arr = img_work.getImageAsNumpy().astype(np.float64)
        ovl_arr = ovl_work.getImageAsNumpy().astype(np.float64)
        ovl_mask = (ovl_arr > 0).astype(np.uint8)

        # ── Find center of mass of overlay ──────────────────────────
        if ovl_mask.max() > 0:
            com = ndimage.center_of_mass(ovl_mask)
            center = tuple(int(round(c)) for c in com)
        else:
            s = ovl_mask.shape
            center = (s[0] // 2, s[1] // 2, s[2] // 2)

        # ── Helper: normalise a 2D slice ────────────────────────────
        def _norm(slc):
            valid = slc[slc > 0]
            if valid.size > 0:
                p1, p99 = np.percentile(valid, [1, 99])
            else:
                p1, p99 = 0.0, 1.0
            if p99 <= p1:
                p99 = p1 + 1.0
            return np.clip((slc - p1) / (p99 - p1), 0.0, 1.0)

        # ── Helper: contour of a 2D mask ────────────────────────────
        def _contour(m2d):
            if m2d.max() == 0 or contour_iterations <= 0:
                return np.zeros_like(m2d, dtype=bool)
            struct = ndimage.generate_binary_structure(2, 2)
            dilated = ndimage.binary_dilation(m2d > 0, structure=struct,
                                              iterations=contour_iterations)
            eroded = ndimage.binary_erosion(m2d > 0, structure=struct,
                                            iterations=max(1, contour_iterations - 1))
            return dilated & ~eroded

        # ── Helper: overlay onto RGB ────────────────────────────────
        def _apply_overlay(rgb, mask2d):
            m = mask2d > 0
            if not m.any():
                return rgb
            for c in range(3):
                rgb[:, :, c] = np.where(
                    m,
                    rgb[:, :, c] * (1.0 - fill_alpha) + overlay_color[c] * fill_alpha,
                    rgb[:, :, c],
                )
            ct = _contour(mask2d)
            if ct.any():
                for c in range(3):
                    rgb[:, :, c] = np.where(
                        ct,
                        rgb[:, :, c] * (1.0 - contour_alpha) + overlay_color[c] * contour_alpha,
                        rgb[:, :, c],
                    )
            return rgb

        # ── Helper: extract & flip 2D slice ─────────────────────────
        def _slice(vol, axis, idx):
            idx = int(np.clip(idx, 0, vol.shape[axis] - 1))
            slc = np.take(vol, idx, axis=axis).astype(np.float64)
            return np.flipud(slc)

        # ── Compute ROI stats ───────────────────────────────────────
        spx = spacing if spacing is not None else list(img_work.getImageSpacing())
        vox_vol_ml = float(np.prod(spx[:3])) / 1000.0
        n_vox = int(ovl_mask.sum())
        vol_ml = round(n_vox * vox_vol_ml, 2)
        if n_vox > 0:
            vals = img_arr[ovl_mask > 0]
            mean_v, std_v = round(float(np.mean(vals)), 2), round(float(np.std(vals)), 2)
        else:
            mean_v = std_v = 0.0
        stats_dict = {'n_voxels': n_vox, 'volume_ml': vol_ml,
                      'mean': mean_v, 'std': std_v}

        # ── Figure layout ───────────────────────────────────────────
        n_offsets = len(slice_offsets)
        n_views = len(views)
        n_cols = n_views * n_offsets + (1 if stats else 0)
        width_ratios = [1.0] * (n_views * n_offsets)
        if stats:
            width_ratios.append(0.45)

        if figsize is None:
            figsize = (n_cols * 4.0, 4.5)

        fig, axes = plt.subplots(
            1, n_cols, figsize=figsize, facecolor='black',
            gridspec_kw={'width_ratios': width_ratios},
        )
        if n_cols == 1:
            axes = np.array([axes])

        if title is None:
            sp_txt = f'{spacing[0]} mm iso' if spacing else 'native'
            title = f'Overlay Report  ({orientation}, {sp_txt})'
        fig.suptitle(title, color='white', fontsize=14, fontweight='bold', y=1.02)

        # ── Render panels ───────────────────────────────────────────
        col = 0
        for view_name in views:
            axis_num = view_axis_map[view_name]
            for offset in slice_offsets:
                idx = center[axis_num] + offset
                slc = _slice(img_arr, axis_num, idx)
                msk = _slice(ovl_mask, axis_num, idx)

                slc_n = _norm(slc)
                # Apply cmap
                cmap_fn = matplotlib.colormaps.get_cmap(image_cmap)
                rgb = cmap_fn(slc_n)[:, :, :3].copy()
                rgb = _apply_overlay(rgb, msk)
                rgb = np.clip(rgb, 0, 1)

                ax = axes[col]
                ax.imshow(rgb, aspect='equal')
                ax.set_facecolor('black')
                ax.set_xticks([])
                ax.set_yticks([])

                lbl = view_name.capitalize()
                if len(slice_offsets) > 1:
                    lbl += f'\n(offset {offset:+d})'
                ax.set_title(lbl, color='white', fontsize=11, fontweight='bold')
                col += 1

        # ── Stats column ────────────────────────────────────────────
        if stats:
            ax_s = axes[col]
            ax_s.set_facecolor('black')
            ax_s.set_xticks([])
            ax_s.set_yticks([])
            for sp in ax_s.spines.values():
                sp.set_visible(False)
            txt = (
                f"Volume\n{stats_dict['volume_ml']:.2f} mL\n"
                f"({stats_dict['n_voxels']} vox)\n\n"
                f"Mean\n{stats_dict['mean']:.1f}\n\n"
                f"Std\n{stats_dict['std']:.1f}"
            )
            ax_s.text(0.5, 0.5, txt, color='white', fontsize=11,
                      fontweight='bold', ha='center', va='center',
                      transform=ax_s.transAxes, family='monospace')
            ax_s.set_title('Stats', color='white', fontsize=11,
                           fontweight='bold')

        plt.tight_layout()

        if save:
            fig.savefig(save, dpi=dpi, facecolor='black', bbox_inches='tight')
        if show:
            plt.show()

        return {'figure': fig, 'axes': axes, 'stats': stats_dict}


    def overlayAble(self,secondimaginable, axis,index,image_cmap='gray', labelmap_cmap='jet', alpha_value=0.5, image_vmin=None, image_vmax=None, labelmap_vmin=None, labelmap_vmax=None,show=False,save=None,title=None,labelmap_name=None,titles=None,figsize=None,colorbar=None,index_mode='auto'):
        """Overlay one or more slices from another Imaginable/LabelMapable.

        If ``axis`` and ``index`` are scalars, a single overlay is drawn on the
        current axes. If either is a list, tuple, range, or NumPy array, each
        requested slice is drawn in a compact ceil(sqrt(n)) by ceil(sqrt(n))
        grid. In ``index_mode='auto'``, multi-axis plus a 3D index point uses
        each axis coordinate from that point. Use ``index_mode='cartesian'`` to
        combine every requested axis with every requested index.

        For multi-slice grids, use ``titles`` for per-panel titles. Passing a
        list/array to ``title`` is also accepted as shorthand for ``titles``;
        passing a scalar string to ``title`` sets the figure title.
        """
        IM=self
        ROI=secondimaginable

        pairs, multi_slice = makeAxisIndexPairs(axis, index, index_mode=index_mode)

        if multi_slice:
            figure_title = title
            panel_titles = titles
            if panel_titles is None and isinstance(title, (list, tuple, np.ndarray)):
                panel_titles = title
                figure_title = None
            if panel_titles is None:
                axes = [axis_value for axis_value, _ in pairs]
                if len(set(axes)) > 1:
                    panel_titles = [f"Axis {axis_value} Slice {index_value}" for axis_value, index_value in pairs]
                else:
                    panel_titles = [f"Slice {index_value}" for _, index_value in pairs]

            images = []
            labelmaps = []
            for axis_value, index_value in pairs:
                images.append(getImaginableSliceNumpy(IM, axis_value, index_value))
                labelmaps.append(getImaginableSliceNumpy(ROI, axis_value, index_value))

            return overlayNumpyImageAndNumpyLabelmapGrid(
                images,
                labelmaps,
                image_cmap=image_cmap,
                labelmap_cmap=labelmap_cmap,
                alpha_value=alpha_value,
                image_vmax=image_vmax,
                image_vmin=image_vmin,
                labelmap_vmax=labelmap_vmax,
                labelmap_vmin=labelmap_vmin,
                show=show,
                save=save,
                title=figure_title,
                titles=panel_titles,
                labelmap_name=labelmap_name,
                figsize=figsize,
                colorbar=False if colorbar is None else colorbar,
            )
        
        if titles is not None and title is None:
            if isinstance(titles, (list, tuple, np.ndarray)):
                flat_titles = np.asarray(titles, dtype=object).ravel()
                title = None if flat_titles.size == 0 else flat_titles[0]
            else:
                title = titles

        axis_value, index_value = pairs[0]
        im = getImaginableSliceNumpy(IM, axis_value, index_value)
        im2 = getImaginableSliceNumpy(ROI, axis_value, index_value)

        # getImaginableSliceNumpy() returns a 2D numpy slice in (Y,X) ordering
        # which is directly compatible with matplotlib.imshow (rows, cols).
        # Previously code used .T here (legacy from v2 conventions) — remove it.
        return overlayNumpyImageAndNumpyLabelmap(
            im, im2,
            image_cmap=image_cmap,
            labelmap_cmap=labelmap_cmap,
            alpha_value=alpha_value,
            image_vmax=image_vmax,
            image_vmin=image_vmin,
            labelmap_vmax=labelmap_vmax,
            labelmap_vmin=labelmap_vmin,
            show=show,
            save=save,
            title=title,
            labelmap_name=labelmap_name,
            colorbar=True if colorbar is None else colorbar
        )


    def overlayAbleImage(self,secondimaginable, axis,index,image_cmap='gray', labelmap_cmap='jet', alpha_value=0.5, image_vmin=None, image_vmax=None, labelmap_vmin=None, labelmap_vmax=None,as_base64=False,data_uri=False,save=None,origin='lower',title=None,titles=None,ncols=None,tile_gap=0,slice_offsets=None,title_font_size=12,title_padding=2,title_color=(255,255,255,255),background=(0,0,0,255),index_mode='auto'):
        """Return only the image+overlay raster for one or more slices.

        Returns an ``(H, W, 4)`` uint8 RGBA array by default. If
        ``as_base64=True``, returns a PNG base64 string instead. Pass a
        list/array/range of axes or indices, or pass ``slice_offsets`` with a
        scalar center index, to create a tight 2.5D montage. In
        ``index_mode='auto'``, multi-axis plus a 3D index point uses each axis
        coordinate from that point. Use ``index_mode='cartesian'`` to combine
        every requested axis with every requested index.
        """
        pairs, multi_slice = makeAxisIndexPairs(axis, index, slice_offsets=slice_offsets, index_mode=index_mode)

        if multi_slice:
            figure_title = title
            panel_titles = titles
            if panel_titles is None and isinstance(title, (list, tuple, np.ndarray)):
                panel_titles = title
                figure_title = None

            images = []
            labelmaps = []
            for axis_value, index_value in pairs:
                images.append(getImaginableSliceNumpy(self, axis_value, index_value))
                labelmaps.append(getImaginableSliceNumpy(secondimaginable, axis_value, index_value))

            return overlayNumpyImageAndNumpyLabelmapGridToImage(
                images,
                labelmaps,
                image_cmap=image_cmap,
                labelmap_cmap=labelmap_cmap,
                alpha_value=alpha_value,
                image_vmax=image_vmax,
                image_vmin=image_vmin,
                labelmap_vmax=labelmap_vmax,
                labelmap_vmin=labelmap_vmin,
                as_base64=as_base64,
                data_uri=data_uri,
                save=save,
                origin=origin,
                title=figure_title,
                titles=panel_titles,
                ncols=ncols,
                tile_gap=tile_gap,
                title_font_size=title_font_size,
                title_padding=title_padding,
                title_color=title_color,
                background=background
            )

        if titles is not None and title is None:
            if isinstance(titles, (list, tuple, np.ndarray)):
                flat_titles = np.asarray(titles, dtype=object).ravel()
                title = None if flat_titles.size == 0 else flat_titles[0]
            else:
                title = titles

        axis_value, index_value = pairs[0]
        im = getImaginableSliceNumpy(self, axis_value, index_value)
        im2 = getImaginableSliceNumpy(secondimaginable, axis_value, index_value)

        return overlayNumpyImageAndNumpyLabelmapToImage(
            im,
            im2,
            image_cmap=image_cmap,
            labelmap_cmap=labelmap_cmap,
            alpha_value=alpha_value,
            image_vmax=image_vmax,
            image_vmin=image_vmin,
            labelmap_vmax=labelmap_vmax,
            labelmap_vmin=labelmap_vmin,
            as_base64=as_base64,
            data_uri=data_uri,
            save=save,
            origin=origin,
            title=title,
            title_font_size=title_font_size,
            title_padding=title_padding,
            title_color=title_color,
            background=background
        )


    def filterValues(self,values):
        L=self.getDuplicate()
        C=sitk.Image(L.getImage())
        C=C*0
        for s in values:
            m=L.getImage()==s
            C+=m*s
        return self.setImage(C)
    
    def cropToBoundingBox(self):
        BL,BU=self.getBoundingBox()
        # getBoundingBox returns (z,y,x) numpy order, cropImage expects ITK (x,y,z)
        BL=[int(a) for a in reversed(BL)]
        BU=[int(a) for a in reversed(BU)]
        return self.cropImage(BL,BU)
    
    def resampleOnCanonicalSpace(self, interpolator=None, useNearestNeighborExtrapolator=None, bgvalue=0.0):
        """
        Resample image to canonical LPS orientation with axis-aligned grid.
        
        This comprehensive method handles both oblique acquisitions and 
        non-LPS orientations, ensuring the result has:
        - Direction matrix: identity (1,0,0, 0,1,0, 0,0,1)
        - Anatomical orientation: LPS (Left-Posterior-Superior)
        
        The method automatically detects if the image is oblique and applies
        the appropriate transformation:
        1. If oblique: resamples to axis-aligned grid first (uses interpolation)
        2. Then reorients to LPS anatomical convention (permute/flip if needed)
        
        Args:
            interpolator: Interpolation method for resampling (default: linear)
            useNearestNeighborExtrapolator: Extrapolator for out-of-bounds (default: False)
            bgvalue: Background value for regions outside original image (default: 0.0)
        
        Returns:
            self (for chaining)
            
        Example:
            >>> img = Imaginable(imagepath='oblique_scan.nii.gz')
            >>> img.resampleOnCanonicalSpace()
            >>> # Now in LPS with identity direction matrix
            >>> arr = img.getImageAsNumpy()  # (Z,Y,X) with predictable axes
            
        Note:
            - For oblique images: resampling with interpolation is applied
            - For axis-aligned images: only permutation/flipping (no interpolation)
            - Physical coordinates (mm) are always preserved
            - Deprecated alias: Use resampleToAxisAligned() + reorientToLPS() for explicit control
        """
        # Resolve interpolator/extrapolator defaults from class attributes
        if interpolator is None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator is None:
            useNearestNeighborExtrapolator = self.dfltuseNearestNeighborExtrapolator

        # Step 1: Handle oblique acquisitions if necessary
        if not self.isAxisAligned():
            # Image is oblique - need to resample to axis-aligned grid first
            self.resampleToAxisAligned(interpolator, useNearestNeighborExtrapolator, bgvalue)
        
        # Step 2: Ensure LPS anatomical orientation
        # If already LPS and axis-aligned, this is a no-op
        self.reorientToLPS()
        
        return self
    
    def setImageFromNumpy(self, nparray, refimage=None, vector=False, spacing=None, origin=None, direction=None):
        """
        Set the image from a numpy array in standard (Z, Y, X) ordering.
        
        Args:
            nparray: Numpy array in (Z, Y, X) order for 3D, (Y, X) for 2D
            refimage: Reference image to copy metadata from
            vector: Whether the array represents a vector image
            spacing: Physical spacing if not using refimage
            origin: Physical origin if not using refimage
            direction: Direction matrix if not using refimage
        
        Returns:
            self
        
        Note: Changed in v3! Now expects (Z,Y,X) ordering to match getImageAsNumpy().
              Use setImageFromNumpyXYZ() if you have (X,Y,Z) ordered arrays.
        
        Example:
            >>> arr = np.random.rand(100, 200, 300)  # (Z, Y, X)
            >>> img.setImageFromNumpy(arr)
        """
        # Array is already in (Z,Y,X), pass directly to sitk.GetImageFromArray
        nda = sitk.GetImageFromArray(nparray, isVector=vector)
        if refimage:
            REF = getmeTheSimpleITKImage(refimage)
            if np.array_equiv(REF.GetSize(), nparray.shape[::-1]):  # ITK size is (X,Y,Z)
                nda.CopyInformation(REF)
            else:
                nda = setSITKImageInforFromImage(nda, REF)
        elif ((spacing) and (origin) and (direction)):
            nda = setSITKImageInfo(nda, spacing=spacing, origin=origin, direction=direction)
        elif self.isImageSet():
            if(self.getImage()):
                r, o, d = getSITKImageInfo(getmeTheSimpleITKImage(self))
                nda = setSITKImageInfo(nda, spacing=r, origin=o, direction=d)            
        self.setImage(nda, 'image set from numpy array (Z,Y,X)!')
        return self

    def setImageFromNumpyZYX(self, nparray, refimage=None, vector=False, spacing=None, origin=None, direction=None):
        """
        Set the image from a numpy array in (Z, Y, X) ordering.
        Alias for setImageFromNumpy() for explicit clarity.
        
        Args:
            nparray: Numpy array in (Z, Y, X) order for 3D, (Y, X) for 2D
            refimage: Reference image to copy metadata from
            vector: Whether the array represents a vector image
            spacing: Physical spacing if not using refimage
            origin: Physical origin if not using refimage
            direction: Direction matrix if not using refimage
        
        Returns:
            self
        """
        return self.setImageFromNumpy(nparray, refimage, vector, spacing, origin, direction)
    
    def setImageFromNumpyXYZ(self, nparray, refimage=None, vector=False, spacing=None, origin=None, direction=None):
        """
        Set the image from a numpy array in (X, Y, Z) ordering.
        
        DEPRECATED: Provided for backward compatibility only.
        The old setImageFromNumpy() expected this ordering in v2.
        
        Args:
            nparray: Numpy array in (X, Y, Z) order - NON-STANDARD
        
        Returns:
            self
        """
        # Transpose from (X,Y,Z) to (Z,Y,X) before setting
        L = list(range(len(nparray.shape)))
        L.reverse()
        o = np.transpose(nparray, L)
        return self.setImageFromNumpy(o, refimage, vector, spacing, origin, direction)
    
    def getImageDirection(self):
        image=self.getImage()
        return image.GetDirection()
    
    def setImageDirection(self,direction):
        image=self.getImage()
        image.SetDirection(direction)
        self.setImage(image)
        return self

    def getImageSpacing(self):
        image=self.getImage()
        return image.GetSpacing()

    def changeImageSpacing(self,spacing,interpolator=None,useNearestNeighborExtrapolator=None,bgvalue=0.0):
        if interpolator == None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator ==None:
            useNearestNeighborExtrapolator=self.dfltuseNearestNeighborExtrapolator

        image=self.getImage()
        d=self.getImageSpacing()
        s=self.getImageSize()
        newSize = [round((sz-1)*spc/s) for sz,spc,s in zip(image.GetSize(), image.GetSpacing(), spacing)]
        t=sitk.Transform()
        t.SetIdentity()
        mess=f'spacing changed from {d} to {spacing}'
        
        self.setImage(sitk.Resample(image, newSize, t,  interpolator, self.getImageOrigin(), spacing, self.getImageDirection(), bgvalue,image.GetPixelIDValue(),useNearestNeighborExtrapolator),mess)
        return self
    
    def changeImageDirection(self, direction, interpolator=None, useNearestNeighborExtrapolator=None, bgvalue=0.0):
        """
        Resample image to a new direction matrix (e.g., axis-aligned grid).
        
        This method resamples the image data onto a new grid with the specified
        direction matrix. Useful for converting oblique acquisitions to axis-aligned
        grids with standard direction cosines like (1,0,0, 0,1,0, 0,0,1).
        
        IMPORTANT: This resamples the data, which may introduce interpolation artifacts.
        The image size, spacing, and origin are preserved, but the voxel values are
        interpolated onto the new grid orientation.
        
        Args:
            direction: Target direction matrix (9-tuple for 3D, 4-tuple for 2D)
                      Example: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) for axis-aligned
            interpolator: Interpolation method (default: linear)
            useNearestNeighborExtrapolator: Whether to use nearest neighbor for out-of-bounds
            bgvalue: Background value for regions outside original image
        
        Returns:
            self (for chaining)
            
        Example:
            >>> # Resample oblique image to axis-aligned grid
            >>> img.changeImageDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
            
        Note:
            - Use dicomOrient() if you want to reorient without resampling (permute/flip only)
            - This method resamples, so it may change voxel values slightly due to interpolation
        """
        if interpolator is None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator is None:
            useNearestNeighborExtrapolator = self.dfltuseNearestNeighborExtrapolator

        image = self.getImage()
        old_direction = self.getImageDirection()
        
        # Create identity transform (no spatial translation/rotation)
        t = sitk.Transform()
        t.SetIdentity()
        
        mess = f'direction changed from {old_direction} to {direction}'
        
        # Resample with new direction
        resampled = sitk.Resample(
            image, 
            self.getImageSize(),  # Keep same size
            t,  # Identity transform
            interpolator, 
            self.getImageOrigin(),  # Keep same origin
            self.getImageSpacing(),  # Keep same spacing
            direction,  # New direction
            bgvalue, 
            image.GetPixelIDValue(),
            useNearestNeighborExtrapolator
        )
        
        self.setImage(resampled, mess)
        return self
    
    def resampleToAxisAligned(self, interpolator=None, useNearestNeighborExtrapolator=None, bgvalue=0.0):
        """
        Resample oblique/rotated image to axis-aligned grid with identity direction matrix.
        
        This is particularly useful for oblique acquisitions (e.g., oblique MRI scans)
        where the direction cosines are not aligned with the standard axes. After this
        operation, the direction matrix will be identity: (1,0,0, 0,1,0, 0,0,1).
        
        The image is resampled so that:
        - Voxel i-axis aligns with physical X-axis (Left→Right in LPS)
        - Voxel j-axis aligns with physical Y-axis (Posterior→Anterior in LPS)
        - Voxel k-axis aligns with physical Z-axis (Inferior→Superior in LPS)
        
        Args:
            interpolator: Interpolation method (default: linear)
            useNearestNeighborExtrapolator: Whether to use nearest neighbor for out-of-bounds
            bgvalue: Background value for regions outside original image
        
        Returns:
            self (for chaining)
            
        Example:
            >>> # Load oblique MRI scan
            >>> img = Imaginable(imagepath='oblique_scan.nii.gz')
            >>> print(img.getDirectionCosines())
            (0.866, 0.5, 0.0, -0.5, 0.866, 0.0, 0.0, 0.0, 1.0)  # Oblique!
            >>> 
            >>> # Resample to axis-aligned
            >>> img.resampleToAxisAligned()
            >>> print(img.getDirectionCosines())
            (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)  # Axis-aligned!
            
        Note:
            - This resamples the data, introducing interpolation
            - Physical coordinates (mm) are preserved
            - Size, spacing, origin maintained
            - For axis-aligned images, this is a no-op (direction already identity)
        """
        # Standard identity direction matrix for 3D
        if self.getImageDimension() == 3:
            identity_direction = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        elif self.getImageDimension() == 2:
            identity_direction = (1.0, 0.0, 0.0, 1.0)
        else:
            raise ValueError(f"Unsupported dimension: {self.getImageDimension()}")
        
        return self.changeImageDirection(identity_direction, interpolator, useNearestNeighborExtrapolator, bgvalue)
    
    def isAxisAligned(self, tolerance=1e-6):
        """
        Check if image has axis-aligned direction cosines (near-identity matrix).
        
        Args:
            tolerance: Tolerance for considering values as 0 or 1
        
        Returns:
            bool: True if direction matrix is close to identity
            
        Example:
            >>> img.isAxisAligned()
            False
            >>> img.resampleToAxisAligned()
            >>> img.isAxisAligned()
            True
        """
        direction = self.getImageDirection()
        
        if self.getImageDimension() == 3:
            expected = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        elif self.getImageDimension() == 2:
            expected = (1.0, 0.0, 0.0, 1.0)
        else:
            return False
        
        # Check if all elements are within tolerance
        for actual, expect in zip(direction, expected):
            if abs(actual - expect) > tolerance:
                return False
        return True

    def dicomOrient(self,orientation='LPS'):
        """
        Reorient the image to a specified anatomical orientation.
        
        This method DOES modify the numpy array data! It physically reorients
        the image volume by permuting and potentially flipping axes so that
        the anatomical directions align with the specified orientation code.
        
        Args:
            orientation: Three-letter code specifying anatomical axes.
                        Common codes:
                        - 'LPS': Left-Posterior-Superior (DICOM standard)
                        - 'RAS': Right-Anterior-Superior (NIfTI/neuroimaging)
                        - 'LAS': Left-Anterior-Superior
                        - 'RPI': Right-Posterior-Inferior
        
        Returns:
            self (for chaining)
            
        Note:
            After calling this method:
            - The numpy array (from getImageAsNumpy()) will be rearranged
            - The direction matrix will be updated to nearly identity
            - Physical world coordinates are preserved (same anatomy in same locations)
            - The image size may change if axes are permuted
        """
        self.setImage(sitk.DICOMOrient(self.getImage(),orientation),'image oriented to '+orientation)
        return self
    
    def reorientToLPS(self):
        """
        Convenience method to reorient image to LPS (Left-Posterior-Superior).
        
        LPS is the DICOM standard orientation where:
        - First axis (X): Left to Right (L→R)
        - Second axis (Y): Posterior to Anterior (P→A)
        - Third axis (Z): Inferior to Superior (I→S)
        
        After calling this, array indices approximately map to:
        - array[k,j,i] where i increases L→R, j increases P→A, k increases I→S
        
        Returns:
            self (for chaining)
        """
        return self.dicomOrient('LPS')
    
    def reorientToRAS(self):
        """
        Convenience method to reorient image to RAS (Right-Anterior-Superior).
        
        RAS is common in neuroimaging (NIfTI) where:
        - First axis (X): Right to Left (R→L)
        - Second axis (Y): Anterior to Posterior (A→P)
        - Third axis (Z): Inferior to Superior (I→S)
        
        Returns:
            self (for chaining)
        """
        return self.dicomOrient('RAS')
    
    def reorientToRPI(self):
        """
        Convenience method to reorient image to RPI (Right-Posterior-Inferior).
        
        RPI orientation where:
        - First axis (X): Right to Left (R→L)
        - Second axis (Y): Posterior to Anterior (P→A)
        - Third axis (Z): Superior to Inferior (S→I)
        
        Returns:
            self (for chaining)
        """
        return self.dicomOrient('RPI')
    
    def getOrientationCode(self):
        """
        Get the current anatomical orientation code of the image.
        
        Returns a three-letter code (e.g., 'LPS', 'RAS') describing which
        anatomical direction each axis increases toward.
        
        Returns:
            str: Three-letter orientation code (e.g., 'LPS', 'RAS', 'RPI')
            
        Example:
            >>> img.getOrientationCode()
            'LPS'  # means X increases Left→Right, Y increases Posterior→Anterior, Z increases Inferior→Superior
        """
        return sitk.DICOMOrientImageFilter_GetOrientationFromDirectionCosines(self.getImageDirection())
    
    def getDirectionCosines(self):
        """
        Get the direction cosine matrix as a tuple.
        
        The direction matrix defines how voxel indices map to physical coordinates.
        For a 3D image, this is a 9-element tuple representing a 3x3 matrix:
        [dir_x0, dir_x1, dir_x2, dir_y0, dir_y1, dir_y2, dir_z0, dir_z1, dir_z2]
        
        Where:
        - (dir_x0, dir_x1, dir_x2) = direction cosines for first axis (columns)
        - (dir_y0, dir_y1, dir_y2) = direction cosines for second axis (rows)
        - (dir_z0, dir_z1, dir_z2) = direction cosines for third axis (slices)
        
        Returns:
            tuple: Direction cosine matrix (9 elements for 3D, 4 for 2D)
            
        Example:
            >>> img.getDirectionCosines()
            (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)  # Identity = LPS orientation
        """
        return self.getImageDirection()
    
    def setDirectionCosines(self, direction):
        """
        Set the direction cosine matrix.
        
        WARNING: This only changes the metadata, NOT the numpy array data!
        If you want to physically reorient the image data, use dicomOrient() instead.
        
        Args:
            direction: Tuple or list of direction cosines
                      (9 elements for 3D, 4 elements for 2D)
        
        Returns:
            self (for chaining)
            
        See Also:
            - dicomOrient(): To physically reorient the image data
            - getDirectionCosines(): To get current direction matrix
        """
        return self.setImageDirection(direction)

    def setImageSpacing(self,spacing):
        image=self.getImage()
        d=self.getImageSpacing()
        image.SetSpacing(spacing)
        self.setImage(image,f'spacing changed from {d} to {spacing}')
    
    def getImageOrigin(self):
        image=self.getImage()
        return image.GetOrigin()
    
    def setImageOrigin(self,o):
        image=self.getImage()
        d=self.getImageOrigin()
        image.SetOrigin(o)
        self.setImage(image,f'origin changed from {d} to {o}')
        
    def getImageDimension(self):
        image=self.getImage()
        return image.GetDimension()
    
    def getImageNumberOfComponentsPerPixel(self):
        image=self.getImage()
        return image.GetNumberOfComponentsPerPixel()

    def getImagePixelTypeAsString(self):
        image=self.getImage()
        return image.GetPixelIDValue(), image.GetPixelIDTypeAsString()

    def getImagePixelTypeAsID(self):
        image=self.getImage()
        return image.GetPixelIDValue()

    def getImageSize(self,index=None):
        image=self.getImage()
        S=list(image.GetSize())
        if not index==None:
            S=S[index]
        return S

         
    
    def getVoxelVolume(self):
        """get the volume of a voxel in the imaginable

        Returns:/discover
            float: voxel volume
        """        
        # image=self.getImage()
        return np.prod(self.getImageSpacing())

    def getNumberOfVoxels(self):
        """get the number of voxels in the imaginable

        Returns:
            int: N voxels
        """        
        return np.prod(self.getImageSize())

    # Backward compatibility alias
    getNumberofVoxels = getNumberOfVoxels

    def getVolume(self):
        """get the volume of the imaginable

        Returns:
            float: volume
        """        
        return self.getVoxelVolume() * self.getNumberOfVoxels()

    def __tellme__(self,m,t=None):
        if t:
            if t=='destroy xxx999':
                return
            
        else:
            self.log.append(m,t)

        if(self.getVerbose()):
            print(m)
    def whathappened(self):
        self.log.getWhatHappened()

    def _print_describe(self, info):
        """Pretty-print a describe() dictionary."""
        max_key = max(len(k) for k in info) if info else 0
        for k, v in info.items():
            print(f"  {k:<{max_key}} : {v}")

    def describe(self):
        """Print a concise summary of the image.

        Returns:
            dict: key image properties
        """
        info = {}
        try:
            info['class'] = type(self).__name__
            img = self.getImage()
            if img is not None:
                info['size (x,y,z)'] = list(img.GetSize())
                info['spacing'] = [round(s, 4) for s in img.GetSpacing()]
                info['origin'] = [round(o, 4) for o in img.GetOrigin()]
                info['dimension'] = img.GetDimension()
                info['pixel_type'] = img.GetPixelIDTypeAsString()
                info['num_voxels'] = int(self.getNumberOfVoxels())
                info['voxel_volume'] = round(float(self.getVoxelVolume()), 6)
                info['total_volume'] = round(float(self.getVolume()), 4)
            else:
                info['image'] = 'Not set'
        except Exception as e:
            info['error'] = str(e)
        self._print_describe(info)
        return info
    def __del__(self):
        self.__tellme__("I'm being automatically destroyed. Goodbye!",'destroy xxx999')

    def forkDuplicate(self):  
        return copy.deepcopy(self)

#    def getDuplicate(self):
#        PP=self.__class__()
#        PP.setImage(self.getImage())
#        return PP

    def getDuplicate(self):
        PP=self.__class__(image=self.getImage())
        # PP.setImage(self.getImage())
        return PP
            
    def resampleOnTargetImage(self,target,interpolator = None,default_value = 0,useNearestNeighborExtrapolator=None):
        """Resample the image on a target image

        Args:
            target (_type_): _description_
            interpolator (_type_, optional): _description_. Defaults to None.
            default_value (int, optional): _description_. Defaults to 0.
            useNearestNeighborExtrapolator (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        target=getmeTheSimpleITKImage(target)
        if interpolator == None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator ==None:
            useNearestNeighborExtrapolator=self.dfltuseNearestNeighborExtrapolator       
        if not default_value:
            default_value=0.0
        theT =sitk.Transform()
        theT.SetIdentity()
        #check if the pixeltype is complex
        if (self.getImagePixelTypeAsID() == sitk.sitkComplexFloat32) or (self.getImagePixelTypeAsID() == sitk.sitkComplexFloat64):
            real_image = sitk.ComplexToReal(self.getImage())
            imag_image = sitk.ComplexToImaginary(self.getImage())
            re=sitk.Resample(real_image,target,theT,interpolator,default_value,sitk.sitkUnknown,useNearestNeighborExtrapolator)
            im=sitk.Resample(imag_image,target,theT,interpolator,default_value,sitk.sitkUnknown,useNearestNeighborExtrapolator)
            self.setImage(sitk.RealAndImaginaryToComplex(re, im),'resampled on target image')
        else:            
            self.setImage(sitk.Resample(self.getImage(),target,theT,interpolator,default_value,sitk.sitkUnknown,useNearestNeighborExtrapolator),'resampled on target image')
        return self
    def getCoordinatesFromIndex(self,P):
        """
        DEPRECATED: Use getPhysicalPointFromITKIndex() for clarity.
        Converts ITK index (i,j,k) = (x,y,z) to physical point (x,y,z) in mm.
        """
        image=self.getImage()
        return image.TransformContinuousIndexToPhysicalPoint(P)
    
    def getIndexFromCoordinates(self,I):
        """
        DEPRECATED: Use getITKIndexFromPhysicalPoint() for clarity.
        Converts physical point (x,y,z) in mm to ITK index (i,j,k) = (x,y,z).
        """
        image=self.getImage()
        return image.TransformPhysicalPointToIndex(I)

        # return image.TransformPhysicalPointToContinuousIndex(I)
    
    def getPhysicalPointFromArrayIndex(self, kji_index):
        """
        Convert array index (k,j,i) = (z,y,x) to physical point (x,y,z) in mm.
        
        Bridges numpy/PyTorch array indexing with ITK physical coordinates.
        
        Args:
            kji_index: Tuple/list of array indices
                      For 3D: (k, j, i) = (z_idx, y_idx, x_idx)
                      For 2D: (j, i) = (y_idx, x_idx)
        
        Returns:
            tuple: Physical coordinates in mm
                  For 3D: (x_mm, y_mm, z_mm)
                  For 2D: (x_mm, y_mm)
        
        Example:
            >>> arr = img.getImageAsNumpy()  # Shape: (Z, Y, X)
            >>> k, j, i = arr.shape[0]//2, arr.shape[1]//2, arr.shape[2]//2
            >>> physical_point = img.getPhysicalPointFromArrayIndex((k, j, i))
            >>> print(f"Center voxel at {physical_point} mm")
        """
        image = self.getImage()
        if self.getImageDimension() == 3:
            i, j, k = float(kji_index[2]), float(kji_index[1]), float(kji_index[0])
            return image.TransformContinuousIndexToPhysicalPoint([i, j, k])
        elif self.getImageDimension() == 2:
            i, j = float(kji_index[1]), float(kji_index[0])
            return image.TransformContinuousIndexToPhysicalPoint([i, j])
        else:
            return image.TransformContinuousIndexToPhysicalPoint([float(kji_index[0])])
    
    def getArrayIndexFromPhysicalPoint(self, xyz_point):
        """
        Convert physical point (x,y,z) in mm to array index (k,j,i) = (z,y,x).
        
        Bridges ITK physical coordinates with numpy/PyTorch array indexing.
        
        Args:
            xyz_point: Tuple/list of physical coordinates in mm
                      For 3D: (x_mm, y_mm, z_mm)
                      For 2D: (x_mm, y_mm)
        
        Returns:
            tuple: Array indices
                  For 3D: (k, j, i) = (z_idx, y_idx, x_idx)
                  For 2D: (j, i) = (y_idx, x_idx)
        
        Example:
            >>> physical_point = (10.5, 20.3, 30.7)  # mm
            >>> k, j, i = img.getArrayIndexFromPhysicalPoint(physical_point)
            >>> arr = img.getImageAsNumpy()
            >>> value = arr[k, j, i]
        """
        image = self.getImage()
        itk_index = image.TransformPhysicalPointToIndex(xyz_point)
        if self.getImageDimension() == 3:
            return (itk_index[2], itk_index[1], itk_index[0])
        elif self.getImageDimension() == 2:
            return (itk_index[1], itk_index[0])
        else:
            return (itk_index[0],)
    
    def getPhysicalPointFromITKIndex(self, ijk_index):
        """
        Convert ITK index (i,j,k) = (x,y,z) to physical point (x,y,z) in mm.
        Same as getCoordinatesFromIndex() but with clearer name.
        
        Args:
            ijk_index: ITK index (i, j, k) = (x_idx, y_idx, z_idx)
        
        Returns:
            tuple: Physical coordinates (x_mm, y_mm, z_mm)
        """
        return self.getImage().TransformContinuousIndexToPhysicalPoint(ijk_index)
    
    def getITKIndexFromPhysicalPoint(self, xyz_point):
        """
        Convert physical point (x,y,z) in mm to ITK index (i,j,k) = (x,y,z).
        Same as getIndexFromCoordinates() but with clearer name.
        
        Args:
            xyz_point: Physical coordinates (x_mm, y_mm, z_mm)
        
        Returns:
            tuple: ITK index (i, j, k) = (x_idx, y_idx, z_idx)
        """
        return self.getImage().TransformPhysicalPointToIndex(xyz_point)
    def changeImageSize(self,newSize,interpolator= None,bgvalue=0.0,useNearestNeighborExtrapolator=None):
        if interpolator == None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator ==None:
            useNearestNeighborExtrapolator=self.dfltuseNearestNeighborExtrapolator
        
        image=self.getImage()
        dimension = image.GetDimension()
 
        reference_origin = image.GetOrigin()
        reference_direction = image.GetDirection()        

        reference_physical_size = np.zeros(image.GetDimension())
        reference_physical_size[:] = [round((sz-1)*spc) if sz*spc>mx  else mx for sz,spc,mx in zip(image.GetSize(), image.GetSpacing(), reference_physical_size)]
        reference_spacing = [ round(phys_sz/(sz-1), self.settings["spacingMinSize"]) for sz,phys_sz in zip(newSize, reference_physical_size) ]
        transform =sitk.Transform()
        transform.SetIdentity()
        mess=f"Image resized from {image.GetSize()} - {image.GetSpacing()} to {newSize} - {reference_spacing}"
        self.setImage(sitk.Resample(image, newSize, transform,  interpolator, reference_origin, reference_spacing, reference_direction, bgvalue,image.GetPixelIDValue(),useNearestNeighborExtrapolator),mess)
        return self
    def padImage(self,lower_padding,upper_padding,padding_value=0):
        """
        Pad an image with a constant value.

        Args:
            lower_padding (_type_): _description_
            upper_padding (_type_): _description_
            padding_value (_type_, optional): _description_. Defaults to 0.
        Returns:
            _type_: _description_
        """
        image=self.getImage()
        # Create a padding filter
        padding_filter = sitk.ConstantPadImageFilter()

        # Set the padding sizes for each dimension
        padding_filter.SetPadLowerBound(lower_padding)  # Lower padding
        padding_filter.SetPadUpperBound(upper_padding)  # Upper padding

        padding_filter.SetConstant(padding_value)

        # Apply the padding filter to the image
        self.setImage(padding_filter.Execute(image),f'image padded by {lower_padding} and {upper_padding} with constant {padding_value}')
        return  self        
    def getPaddedImage(self, padding, padding_value=0):
        """
        Return a padded copy of the current image.

        Parameters
        ----------
        padding : int or sequence of int
            Symmetric padding applied to the lower and upper bounds.
        padding_value : int or float, optional
            Constant used to fill the padded region.

        Returns
        -------
        Imaginable
            A duplicate image with the requested padding applied.
        """
        if isinstance(padding, numbers.Number):
            padding = [int(padding)] * self.getImageDimension()
        padded = self.getDuplicate()
        padded.padImage(list(padding), list(padding), padding_value=padding_value)
        return padded
    def cropImage(self,lowerB,upperB,coordinates=None):
        PP="voxels"
        image=self.getImage()
        U=self.getImageSize()
        if coordinates:
            PP="coordinates"
            lowerB=list(self.getIndexFromCoordinates(lowerB))
            upperB=list(self.getIndexFromCoordinates(upperB))
        

        for t in range(len(upperB)):
            if (((upperB[t]==0) and (isinstance(upperB[t],int))) or (np.isnan(upperB[t]) if isinstance(upperB[t], (int, float)) else False)):
                upperB[t]=U[t]

       
        S=[u-l for l,u in zip(lowerB,upperB)]
        crop = sitk.ExtractImageFilter()
        crop.SetSize(S)
        crop.SetIndex(lowerB)
        cropped_image = crop.Execute(image)
        m=f'image cropped to {PP} {lowerB} {upperB}'
        return self.setImage(cropped_image,m)
 
    def __transformImage__(self,transform,interpolator = None,reference_image=None ,default_value = 0,useNearestNeighborExtrapolator=None,fromregistration=False):
    # Output image Origin, Spacing, Size, Direction are taken from the reference
    # image in this call to Resample
        if interpolator == None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator ==None:
            useNearestNeighborExtrapolator=self.dfltuseNearestNeighborExtrapolator       
        if not default_value:
            default_value=0.0
        theT =sitk.Transform()
        if (isinstance(transform, tuple) or isinstance(transform,list)):
            for t in transform:
                theT = sitk.CompositeTransform([theT,t])
            transform = theT

        if not reference_image:
            reference_image=self.getImage()

        if fromregistration:
            return sitk.Resample(self.getImage(), reference_image, transform, interpolator, default_value,sitk.sitkUnknown,useNearestNeighborExtrapolator)
        else:
            return sitk.Resample(self.getImage(), reference_image, transform.GetInverse(), interpolator, default_value,sitk.sitkUnknown,useNearestNeighborExtrapolator)



    def translateImage(self,T,interpolator = sitk.sitkLinear,reference_image=None ,default_value = 0,useNearestNeighborExtrapolator=False):
        dimension=self.getImageDimension()
        translation = sitk.TranslationTransform(dimension, T)
        return self.setImage(self.__transformImage__(translation,interpolator,reference_image,default_value,useNearestNeighborExtrapolator),f"tranlated of {T}")

    def scaleImage(self,S,center=None,centerindex=False,interpolator = None,reference_image=None ,default_value = 0,useNearestNeighborExtrapolator=None):
        dimension=self.getImageDimension()
        if interpolator == None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator ==None:
            useNearestNeighborExtrapolator=self.dfltuseNearestNeighborExtrapolator

        if not center:
            center=self.getImageCenterCoordinate()
        else:
            if centerindex:
                center=self.getCoordinatesFromIndex(center)
        theS=[1]*dimension
        for i,s in enumerate(S):
            if (isinstance(s,int) and (s==0)):
                theS[i]=1
            else:
                theS[i]=s
        transform = sitk.ScaleTransform(dimension)
        transform.SetScale(theS)
        transform.SetCenter(center)
        return self.setImage(self.__transformImage__(transform,interpolator,reference_image,default_value,useNearestNeighborExtrapolator),f"scaled of {theS}")

    def transform(self,T,interpolator = None,reference_image=None ,default_value = 0,useNearestNeighborExtrapolator=None):
        if interpolator == None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator ==None:
            useNearestNeighborExtrapolator=self.dfltuseNearestNeighborExtrapolator
        return self.setImage(self.__transformImage__(T,interpolator,reference_image,default_value,useNearestNeighborExtrapolator),f"tranlated of {T}")

    def transformImageAffine(self,A,translation=None,center=[],centerindex=False,interpolator = None,reference_image=None ,default_value = 0.0,useNearestNeighborExtrapolator=None):
        if interpolator == None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator ==None:
            useNearestNeighborExtrapolator=self.dfltuseNearestNeighborExtrapolator
        dimension=self.getImageDimension()
        
        if not center:
            center=self.getImageCenterCoordinate()
        else:
            if centerindex:
                center=self.getCoordinatesFromIndex(center)

        transform = sitk.AffineTransform(dimension)
        transform.SetMatrix(A)
        transform.SetCenter(center)
        if translation is not None:
            transform.SetTranslation(translation)
        return self.setImage(self.__transformImage__(transform,interpolator,reference_image,default_value),f"scaled of {A}")

   
   
    def getImageCenterIndex(self):
        return [round(t/2) for t in self.getImageSize()]

    def getImageCenterCoordinate(self):
        return self.getCoordinatesFromIndex(self.getImageCenterIndex())

    def rotateImage(self,rotation=None,center=None, centerindex=False,translation=None,interpolator = None,reference_image=None ,default_value = 0.0,useNearestNeighborExtrapolator=None, angle=None):

        """
        This function rotates an image across each of the x, y, z axes by theta_x, theta_y, and theta_z degrees
        respectively and resamples it to be isotropic.
        :param image: An sitk 3D image
        :param angles: [theta_x in degrees, theta_y in degrees, param theta_z in degrees]
        :center: center of rotation (position of the center in mm)
        :centerindex: False if center is in cm True if it's in voxels 
        :translation: mm of translation 

        :return: The rotated image
        
        """

        dimension = self.getImageDimension()
        if rotation is None:
            if angle is None:
                raise ValueError("Either rotation or angle must be provided.")
            rotation = [0.0] * dimension
            rotation[min(2, dimension - 1)] = float(angle)
        elif isinstance(rotation, numbers.Number):
            scalar_rotation = float(rotation)
            rotation = [0.0] * dimension
            rotation[min(2, dimension - 1)] = scalar_rotation

        if interpolator == None:
            interpolator = self.dfltInterpolator
        if useNearestNeighborExtrapolator ==None:
            useNearestNeighborExtrapolator=self.dfltuseNearestNeighborExtrapolator
        if not center:
            center=[int(x) for x in self.getImageCenterCoordinate()]
        else:
            if centerindex:
                center=self.getCoordinatesFromIndex(center)

        R=[np.deg2rad(r) for r in rotation]
        if not translation:
            T=[0.0]*self.getImageDimension()
        else:
            T=translation
        if dimension==3:
            
            transform = sitk.Euler3DTransform (center,R[0], 
                                                R[1], 
                                                R[2],T)

        elif dimension==2:
            transform = sitk.Euler2DTransform ()
            transform.SetAngle(R[0])
            transform.SetCenter(center)
            transform.SetTranslation(T)
        
        

        return self.setImage(self.__transformImage__(transform,interpolator,reference_image,default_value,useNearestNeighborExtrapolator),f"rotated of {R} and {T}")

    # ========================================================================
    # DEFORMATION & REGISTRATION METHODS
    # ========================================================================

    def applyTransform(self, transform, target_image=None, interpolator=None, default_value=0):
        """
        Apply a registration transform to the image, with intelligent handling based on image type.

        This universal method works on all Imaginable subclasses:
        - For integer/label images (Roiable, LabelMapable): uses label-preserving transform (nearest-neighbor)
        - For continuous images: uses the standard transform with interpolation
        - For vector fields (Fieldable): handles vector-aware transforms
        - Optionally resamples the result to a target image geometry

        Parameters
        ----------
        transform : str or sitk.Transform
            Path to transform file (.tfm, .h5) or SimpleITK Transform object
        target_image : str or sitk.Image, optional
            Target geometry to resample the warped image into. If provided,
            the transformed image will be resampled to match this reference image's
            origin, spacing, size, and direction. Useful for registration output
            that needs to be in a specific space.
        interpolator : str, optional
            Interpolation method: 'linear', 'nearest', 'gaussian', 'bspline'.
            If None, uses default interpolator. For label images, this is ignored
            in favor of nearest-neighbor to preserve label integrity.
        default_value : float, default=0
            Pixel value for regions outside the image domain

        Returns
        -------
        self : Imaginable
            Self for method chaining

        Example
        -------
        >>> # Continuous image with interpolation
        >>> img = SITKImaginable('moving.nii.gz')
        >>> img.applyTransform('transform.tfm', interpolator='linear')
        >>> img.write('warped.nii.gz')
        
        >>> # Label/ROI with label-preserving transform
        >>> roi = Roiable('segmentation.nii.gz')
        >>> roi.applyTransform('transform.tfm')  # Automatically uses nearest-neighbor
        >>> roi.write('warped_roi.nii.gz')
        
        >>> # Resample to reference space
        >>> moving = SITKImaginable('moving.nii.gz')
        >>> fixed = SITKImaginable('fixed.nii.gz')
        >>> moving.applyTransform('transform.tfm', target_image=fixed.getImage())
        >>> moving.write('warped_to_fixed.nii.gz')
        """
        from . import deformations

        image = self.getImage()
        # Defensive: if no image is set, just return self
        if image is None:
            return self

        # Determine pixel type and prefer label-preserving transform for
        # common integer types.
        try:
            pixid = image.GetPixelID()
        except Exception:
            pixid = None

        integer_pixel_types = {
            sitk.sitkUInt8, sitk.sitkInt8,
            sitk.sitkUInt16, sitk.sitkInt16,
            sitk.sitkUInt32, sitk.sitkInt32
        }

        warped = None
        if pixid in integer_pixel_types:
            # Label/ROI image: use label-preserving transform
            try:
                warped = deformations.apply_transform_to_labels(
                    image,
                    transform,
                    target_image=target_image
                )
            except Exception:
                # Fallback to generic transform if label-specific call fails
                if interpolator is None:
                    interpolator_map = {
                        sitk.sitkLinear: 'linear',
                        sitk.sitkNearestNeighbor: 'nearest',
                        sitk.sitkGaussian: 'gaussian',
                        sitk.sitkBSpline: 'bspline',
                    }
                    interpolator = interpolator_map.get(self.dfltInterpolator, 'linear')
                
                warped = deformations.apply_transform(
                    image,
                    transform,
                    target_image=target_image,
                    interpolator=interpolator,
                    default_pixel_value=default_value
                )
        else:
            # Continuous image: use regular transform with interpolation
            if interpolator is None:
                interpolator_map = {
                    sitk.sitkLinear: 'linear',
                    sitk.sitkNearestNeighbor: 'nearest',
                    sitk.sitkGaussian: 'gaussian',
                    sitk.sitkBSpline: 'bspline',
                }
                interpolator = interpolator_map.get(self.dfltInterpolator, 'linear')
            
            warped = deformations.apply_transform(
                image,
                transform,
                target_image=target_image,
                interpolator=interpolator,
                default_pixel_value=default_value
            )

        return self.setImage(warped, f"applied transform from {transform if isinstance(transform, str) else 'transform object'}")

    def applyDisplacementField(self, displacement_field, target_image=None, interpolator=None, default_value=0):
        """
        Apply a displacement field to warp the image.

        Displacement fields can come from registration algorithms like ANTs, elastix, or custom tools.

        Parameters
        ----------
        displacement_field : str or sitk.Image
            Path to displacement field file (.mha, .nii.gz) or SimpleITK vector image
        target_image : str or sitk.Image, optional
            Target geometry reference. If None, uses displacement field geometry
        interpolator : str, optional
            Interpolation method: 'linear', 'nearest', 'gaussian', 'bspline'.
            If None, uses default interpolator
        default_value : float, default=0
            Pixel value for regions outside image domain

        Returns
        -------
        self : Imaginable
            Self for method chaining

        Example
        -------
        >>> img = SITKImaginable('moving.nii.gz')
        >>> img.applyDisplacementField('deformation.mha', target_image='fixed.nii.gz')
        >>> img.write('warped.nii.gz')
        """
        from . import deformations
        
        if interpolator is None:
            interpolator_map = {
                sitk.sitkLinear: 'linear',
                sitk.sitkNearestNeighbor: 'nearest',
                sitk.sitkGaussian: 'gaussian',
                sitk.sitkBSpline: 'bspline',
            }
            interpolator = interpolator_map.get(self.dfltInterpolator, 'linear')
        
        warped = deformations.apply_deformation_field(
            self.getImage(),
            displacement_field,
            target_image=target_image,
            interpolator=interpolator,
            default_pixel_value=default_value
        )
        
        return self.setImage(warped, f"applied displacement field from {displacement_field if isinstance(displacement_field, str) else 'field object'}")

    def warpImage(self, displacement_field, **kwargs):
        """
        Alias for applyDisplacementField. Warp image using a displacement field.

        Parameters
        ----------
        displacement_field : str or sitk.Image
            Displacement field
        **kwargs
            Additional arguments passed to applyDisplacementField

        Returns
        -------
        self : Imaginable
        """
        return self.applyDisplacementField(displacement_field, **kwargs)

    def alignGeometry(self, reference_image):
        """
        Align image geometry (origin, spacing, direction) to match a reference image.

        Useful for fixing displacement fields or images with incorrect metadata that
        was lost during processing.

        Parameters
        ----------
        reference_image : str or sitk.Image
            Reference image with correct geometry

        Returns
        -------
        self : Imaginable
            Self for method chaining

        Example
        -------
        >>> df = SITKImaginable('deform.mha')
        >>> fixed = SITKImaginable('fixed.nii.gz')
        >>> df.alignGeometry(fixed.getImage())
        >>> df.write('deform_aligned.mha')
        """
        from . import deformations
        
        aligned = deformations.align_geometry(self.getImage(), reference_image)
        return self.setImage(aligned, "geometry aligned to reference")

    def invertDisplacementField(self, max_iterations=100, mean_error_tolerance=1e-3):
        """
        Invert the displacement field for reverse warping.

        Useful for forward-backward consistency checks and inverse transformations.

        Parameters
        ----------
        max_iterations : int, default=100
            Maximum iterations for inversion algorithm
        mean_error_tolerance : float, default=1e-3
            Tolerance for convergence

        Returns
        -------
        self : Imaginable
            Self with inverted displacement field

        Example
        -------
        >>> df = SITKImaginable('forward_deform.mha')
        >>> df.invertDisplacementField()
        >>> df.write('backward_deform.mha')
        """
        from . import deformations
        
        inverted = deformations.invert_displacement_field(
            self.getImage(),
            max_iterations=max_iterations,
            mean_error_tolerance=mean_error_tolerance
        )
        
        return self.setImage(inverted, "displacement field inverted")

    def convertTransformToField(self, transform):
        """
        Convert any registration transform to a displacement field on this image's grid.

        Converts rigid, affine, B-spline, or composite transforms into a dense
        displacement field matching the geometry of this image.

        Parameters
        ----------
        transform : str or sitk.Transform
            Path to transform file (.tfm, .h5, .txt) or SimpleITK Transform object.

        Returns
        -------
        Fieldable
            Vector image (displacement field) on this image's grid.

        Example
        -------
        >>> img = SITKImaginable('image.nii.gz')
        >>> df = img.convertTransformToField('registration.tfm')
        >>> df.write('displacement_field.mha')
        """
        from . import deformations
        image = self.getImage()
        field = deformations.transform_to_displacement_field(
            transform,
            output_size=image.GetSize(),
            output_origin=image.GetOrigin(),
            output_spacing=image.GetSpacing(),
            output_direction=image.GetDirection(),
        )
        return Fieldable(image=field)

    def composeTransforms(self, transforms, inverse_flags=None):
        """
        Create a composite transform from multiple transforms and apply it.

        Parameters
        ----------
        transforms : list of str or sitk.Transform
            Transforms to chain (applied in order).
        inverse_flags : list of bool, optional
            If provided, each True entry inverts the corresponding transform.

        Returns
        -------
        self
            For method chaining.

        Example
        -------
        >>> img.composeTransforms(['rigid.tfm', 'bspline.tfm'])
        """
        from . import deformations
        composite = deformations.create_composite_transform(transforms, inverse_flags)
        return self.applyTransform(composite)

    def changePixelType(self,dtype):
        return self.setImage(sitk.Cast(self.getImage(),dtype),f'casted to {dtype}')
    def cast(self,dtype):
        if dtype=="float":
            dtype=sitk.sitkFloat32
        if dtype=="int32":
            dtype=sitk.sitkInt32
        if dtype=="uint32":
            dtype=sitk.sitkUInt32
        if dtype=="uint8":
            dtype=sitk.sitkUInt8
        if dtype=="int8":
            dtype=sitk.sitkInt8
        if dtype=="int16":
            dtype=sitk.sitkInt16
        if dtype=="uint16":
            dtype=sitk.sitkUInt16
        if dtype=="float64":
            dtype=sitk.sitkFloat64
        if dtype=="complex":
            dtype=sitk.sitkComplexFloat32
        if dtype=="complex64":
            dtype=sitk.sitkComplexFloat64
        if dtype=="complex32":
            dtype=sitk.sitkComplexFloat32            
        return self.changePixelType(dtype)
    def getPossiblePixelTypes(self):
        return "complex32,complex64,float,float64,int16,int32,int8,uint16,uint32,uint8"
    
    def getNumberOfNonZeroVoxels(self):
        return np.count_nonzero(self.getImageAsNumpyZYX())

    def applyAbs(self):
        m=self.getImage()
        ABS=sitk.AbsImageFilter()
        self.setImage(ABS.Execute(m))
        return self


    def applyModulus(self):
        m=self.getImage()
        try:
            ABS=sitk.ComplexToModulusImageFilter()
            self.setImage(ABS.Execute(m))
        except:
            try:
                ABS2=sitk.ModulusImageFilter()
                self.setImage(ABS2.Execute(m))
            except:
                raise Exception("can't apply abs")
        return self
    import numbers

    def add(self, toadd):
        w="Image"
        if not isinstance(toadd,numbers.Number):
            w=str(toadd)
        return self.__filterSelfAndImage__(sitk.AddImageFilter(),toadd,f'add {w}')

    def addImage(self, toadd):
        """Backward-compatible alias for :meth:`add`."""
        return self.add(toadd)

    def multiply(self, toadd):
        w="Image"
        if not isinstance(toadd,numbers.Number):
            w=str(toadd)
        return self.__filterSelfAndImageMat__(sitk.MultiplyImageFilter(),toadd,f'multiply {w}')

    def multiplyImage(self, toadd):
        """Backward-compatible alias for :meth:`multiply`."""
        return self.multiply(toadd)



    def subtract(self, toadd):
        w="Image"
        if not isinstance(toadd,numbers.Number):
            w=str(toadd)
        return self.__filterSelfAndImage__(sitk.SubtractImageFilter(),toadd,f'subtract {w}')

    def subtractImage(self, toadd):
        """Backward-compatible alias for :meth:`subtract`."""
        return self.subtract(toadd)

    def divide(self, toadd):
        w="Image"

        if not isinstance(toadd,numbers.Number):
            w=str(toadd)
        return self.__filterSelfAndImageMat__(sitk.DivideImageFilter(),toadd,f'divide {w}')

    def divideImage(self, toadd):
        """Backward-compatible alias for :meth:`divide`."""
        return self.divide(toadd)


    def __filterSelfAndImage__(self,filter,toadd,message):
        if not isinstance(toadd,numbers.Number):
            sitkimage=getmeTheSimpleITKImage(toadd)
            L=SITKImaginable(image=sitkimage)
            L.resampleOnTargetImage(self.getImage())
            toadd=L.getImage()
            
        else:
            toadd=float(toadd)
        try:
            self.setImage(filter.Execute(self.getImage(),toadd),message)
        except:
            raise Exception(f"Can't {message}")
        return self

    def __filterSelfAndImageMat__(self,filter,toadd,message):
        #get original pixel id
        S=self.getDuplicate()
        O=self.getImagePixelTypeAsID()
        # cast to float
        S.changePixelType(sitk.sitkFloat32)

        if not isinstance(toadd,numbers.Number):
            sitkimage=getmeTheSimpleITKImage(toadd)
            L=SITKImaginable(image=sitkimage)
            L.resampleOnTargetImage(self.getImage())
            L.changePixelType(sitk.sitkFloat32)
            toadd=L.getImage()
            
        else:
            toadd=float(toadd)
        try:
            OUT=filter.Execute(S.getImage(),toadd)
            S.setImage(OUT)
            S.changePixelType(O)
            self.setImage(S.getImage(),message)
        except:
            raise Exception(f"Can't {message}")
        return self
    
    def getCornersCoordinates(self):
        
        return get_image_corners_coordinates(self.getImage())
        
    def isInsidePoint(self,P):
        size=self.getImageSize()
        V=self.getIndexFromCoordinates(P)
        #pretend it's inside
        O=True 

        for a in range(len(P)):
            if (V[a]<0 ) | (V[a]>size[a]):
                O=False
                break

        return O

    def isInsideIndex(self,V):
        size=self.getImageSize()
        #pretend it's inside
        O=True 

        for a in range(len(V)):
            if (V[a]<0 ) | (V[a]>size[a]):
                O=False
                break

        return O
    
    def getImageUniqueValues(self,exclude=[]):
        O=set(np.unique(self.getImageAsNumpy().flatten()))
        for e in exclude:
            O.discard(e)
        return O

    def getMaximumValue(self):
        image=self.getImage()
        filter = sitk.MinimumMaximumImageFilter()
        filter.Execute(image)
        return filter.GetMaximum()
    
    def getMeanValue(self):
        image=self.getImage()
        filter = sitk.StatisticsImageFilter()
        filter.Execute(image)
        return filter.GetMean()

    def getVariance(self):
        image=self.getImage()
        filter = sitk.StatisticsImageFilter()
        filter.Execute(image)
        return filter.GetVariance()
    
    def getSum(self):
        image=self.getImage()
        filter = sitk.StatisticsImageFilter()
        filter.Execute(image)
        return filter.GetSum()
        
    
    def getRoiableValuesUpper(self,th):
        return Roiable(image=self.getImage()>th)



    def getBoundingBox(self, exclude=[0]):
        """
        Returns the bounding box of non-excluded values in array index space.
        
        Returns min and max indices in numpy/array ordering: (k, j, i) = (z, y, x)
        
        Args:
            exclude: List of values to exclude (e.g., background values). Default: [0]
        
        Returns:
            tuple: ((k_min, j_min, i_min), (k_max, j_max, i_max)) as numpy arrays
                  Returns None if no non-excluded values found
        
        Example:
            >>> bbox = img.getBoundingBox(exclude=[0])
            >>> if bbox is not None:
            ...     (k_min, j_min, i_min), (k_max, j_max, i_max) = bbox
            ...     arr = img.getImageAsNumpy()
            ...     cropped = arr[k_min:k_max+1, j_min:j_max+1, i_min:i_max+1]
        
        Note: Changed in v3! Now returns indices in (k,j,i) = (z,y,x) order to match
              numpy arrays. Previously returned (x,y,z) order.
        """
        N = self.getImageAsNumpy()  # Now returns (Z,Y,X) in v3
        mask = np.isin(N, exclude, invert=True)
        roi_indices = np.argwhere(mask)
        
        if len(roi_indices) == 0:
            return None
        
        min_coords = np.min(roi_indices, axis=0)
        max_coords = np.max(roi_indices, axis=0)
        bounding_box = (min_coords, max_coords)

        return bounding_box

    def getStdValue(self):
        image=self.getImage()
        filter = sitk.StatisticsImageFilter()
        filter.Execute(image)
        return filter.GetSigma()

    def getMinimumValue(self):
        image=self.getImage()
        filter = sitk.MinimumMaximumImageFilter()
        filter.Execute(image)
        return filter.GetMinimum()


    def __RegionExtractor__(self,size,index):
        Extractor = sitk.ExtractImageFilter()
        Extractor.SetSize(size)
        Extractor.SetIndex(index)
        return Extractor.Execute(self.getImage())

    def instantiateAnotherAble(self):
        able = self.__class__
        return able()

    def getSliceNormalKAsNumpy(self,slice):
        O=self.instantiateAnotherAble()
        out=self.getSliceNormalK(slice)
        O.setImage(out)
        return O.getImageAsNumpy()

    def getSliceNormalJAsNumpy(self,slice):
        O=self.instantiateAnotherAble()
        out=self.getSliceNormalJ(slice)
        O.setImage(out)
        return O.getImageAsNumpy()

    def getSliceNormalIAsNumpy(self,slice):
        O=self.instantiateAnotherAble()
        out=self.getSliceNormalI(slice)
        O.setImage(out)
        return O.getImageAsNumpy()

    def getSliceNormalK(self,slice):
        slice = int(slice)
        size = list(self.getImageSize())
        size[2] = 0
        index = [0, 0, slice]
        out = self.__RegionExtractor__(size,index)
        return out
    def getSliceNormalJ(self,slice):
        slice = int(slice)
        size = list(self.getImageSize())
        size[1] = 0
        index = [0, slice,0]
        out = self.__RegionExtractor__(size,index)
        return out

    def getSliceNormalI(self,slice):
        slice = int(slice)
        size = list(self.getImageSize())
        size[0] = 0
        index = [slice,0,0]
        out = self.__RegionExtractor__(size,index)
        return out

    def isImaginable(self):
        if isinstance(self,SITKImaginable):
            return True
        for b in self.__class__.__bases__:
            if isinstance(b,SITKImaginable):
                return True
        try:            
            return isinstance(self.getImage(),sitk.Image)
        except:
            return False
    def isSITKImaginable(self):
        return self.isImaginable()

    def viewK(self,km=[True,True]):
        fig, ax = plt.subplots(1, 1)
        tracker = IndexViewer(ax,self,km)
        fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
        fig.canvas.mpl_connect('button_press_event', tracker.onclick)
        plt.show()
    def viewJ(self,km=[True,True]):
        fig, ax = plt.subplots(1, 1)
        tracker = IndexViewer(ax,self,km,normal=1)
        fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
        fig.canvas.mpl_connect('button_press_event', tracker.onclick)
        plt.show()
    def viewI(self,km=[True,True]):
        fig, ax = plt.subplots(1, 1)
        tracker = IndexViewer(ax,self,km,normal=0)
        fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
        fig.canvas.mpl_connect('button_press_event', tracker.onclick)
        plt.show()
    def viewAxial(self):
        if self.getImageDimension()==3:
            k,km=self.__getOrientation__()
            if k==[0,1,2]:
                self.viewK([km[0],km[1]])
        elif self.getImageDimension()==2:
            self.view2D()
    def view2D(self):
            o=self.getImageAsNumpyZYX()
            # if self.UD:
            #     o=np.flipud(o)
            # if self.LR:
            #     o=np.fliplr(o)
            fig, ax = plt.subplots(1, 1)
            im = ax.imshow(o)
            ratio=self.getImageSpacing()
            ratio=ratio[0]/ratio[1]
            ax.set_aspect(ratio)
    def __getOrientation__(self):
        L=self.getImageDirection()

        k=[np.argmax(np.abs(L[0:3])),np.argmax(np.abs(L[3:6])),np.argmax(np.abs(L[6:]))]
        km=list(map(lambda x: x>0, [L[k[0]],L[k[1]+3],L[k[2]+6]]))
        return k, km

    def viewSagittal(self):
        k,km=self.__getOrientation__()
        if k==[0,1,2]:
            self.viewI([km[1],km[2]])

    def viewCoronal(self):
        k,km=self.__getOrientation__()
        if k==[0,1,2]:
            self.viewJ([km[0],km[1]])
    # registration parameters file are written already to be used with the resampling filter so not need to inverse the transform
    # this is good since there's no inverse at the moment for bsplinetransform

    def transformFromRegitration(self,T,interpolator = sitk.sitkLinear,reference_image=None ,default_value = 0,useNearestNeighborExtrapolator=False):
        return self.setImage(self.__transformImage__(T,interpolator,reference_image,default_value,useNearestNeighborExtrapolator,fromregistration=True),f"transformed of {T}")
    
    # def mask(self,mask):
    #     mask=getmeTheSimpleITKImage(mask)
        
    #     self.setImage(,'masked with image')
    def __applyImageToImageFilter__(self,f,cast=False):
        if cast:
            S=self.getDuplicate()
            O=self.getImagePixelTypeAsID()
            # cast to float
            S.changePixelType(sitk.sitkFloat32)
            S.setImage(f.Execute(S.getImage()))
            S.changePixelType(O)
            self.setImage(S.getImage())
        else:
            self.setImage(f.Execute(self.getImage()))

        return self

    def sharpen(self):
        f = sitk.LaplacianSharpeningImageFilter()
        return self.__applyImageToImageFilter__(f)
    
    def denoise(self,timestep=None,numberOfIterations=10,stencilRadius=None):
        # https://itk.org/Doxygen/html/Examples_2Filtering_2CurvatureFlowImageFilter_8cxx-example.html
        if stencilRadius==None:
            f = sitk.CurvatureFlowImageFilter()
            cast=False
        else:
            f=sitk.MinMaxCurvatureFlowImageFilter()
            cast=True
        f.SetNumberOfIterations(numberOfIterations)
        if timestep==None:
            if self.getImageDimension()==3:
                timestep=0.0625
            elif self.getImageDimension()==2:
                timestep=0.125
            else:
                raise Exception("please set a timestep")
        f.SetTimeStep(timestep)
        return self.__applyImageToImageFilter__(f,cast=cast)
    
    def plotOverlay(self, overlay=None, alpha=0.5, title=None, slice_idx=None, **kwargs):
        """
        Display image with optional overlay using interactive viewer.
        
        Parameters
        ----------
        overlay : sitk.Image or Imaginable, optional
            Overlay image (will be resampled to match this image)
        alpha : float, default=0.5
            Overlay opacity (0-1)
        title : str, optional
            Figure title
        slice_idx : int, optional
            Slice index for 3D images (middle slice if None)
        **kwargs
            Additional arguments passed to viewer
        
        Returns
        -------
        viewer : PlotViewer
            Viewer instance
        
        Examples
        --------
        >>> img = Imaginable('image.nii.gz')
        >>> overlay = Imaginable('segmentation.nii.gz')
        >>> img.plotOverlay(overlay, alpha=0.6)
        """
        try:
            from .plotable import plotOverlay
        except ImportError:
            from plotable import plotOverlay
        
        if title is None:
            title = "Image Viewer"
        
        return plotOverlay(self, overlay=overlay, alpha=alpha, 
                          title=title, slice_idx=slice_idx, **kwargs)

    def viewInteractive(self, overlays=None, orientation=2, slice_idx=None, 
                       title=None, figsize=(14, 10), cmap='gray'):
        """
        Open an interactive GUI viewer with advanced controls.
        
        Features:
        - Orientation selection (axial, sagittal, coronal)
        - Slice navigation with auto-center
        - Multiple overlay layers with individual opacity control
        - Real-time updates
        
        Parameters
        ----------
        overlays : Imaginable, array, or list, optional
            Single or multiple overlays to display on top of main image
            Can be:
            - Single Imaginable or array
            - List of [Imaginable/array, ...]
        orientation : int, default=2
            Initial viewing orientation:
            - 0: Axial (XY plane)
            - 1: Sagittal (YZ plane)
            - 2: Coronal (XZ plane)
        slice_idx : int, optional
            Initial slice index. If None, uses center slice.
        title : str, optional
            Window title. Auto-generated if None.
        figsize : tuple, default=(14, 10)
            Figure size in inches (width, height)
        cmap : str, default='gray'
            Colormap for primary image
            
        Returns
        -------
        InteractiveViewer
            Viewer instance (can be used to update/manipulate viewer)
            
        Example
        -------
        **Single overlay:**
        
        >>> img = Imaginable('image.nii.gz')
        >>> seg = Imaginable('segmentation.nii.gz')
        >>> img.viewInteractive(overlays=seg, orientation=2)
        
        **Multiple overlays:**
        
        >>> img = Imaginable('image.nii.gz')
        >>> seg1 = Imaginable('seg1.nii.gz')
        >>> seg2 = Imaginable('seg2.nii.gz')
        >>> img.viewInteractive(overlays=[seg1, seg2], orientation=0)
        
        **With numpy arrays:**
        
        >>> img = Imaginable('image.nii.gz')
        >>> mask = np.zeros((256, 256))
        >>> mask[50:200, 50:200] = 1
        >>> img.viewInteractive(overlays=mask, slice_idx=100)
        """
        try:
            from .interactive_viewer import InteractiveViewer
        except ImportError:
            from interactive_viewer import InteractiveViewer
        
        if title is None:
            title = "Image Viewer - Interactive"
        
        viewer = InteractiveViewer(self.getImage(), title=title, 
                                   figsize=figsize, cmap=cmap)
        viewer.current_orientation = orientation
        viewer.current_slice = slice_idx or viewer._get_center_slice(orientation)
        
        # Add overlays
        if overlays is not None:
            if isinstance(overlays, (list, tuple)):
                viewer.add_overlays(overlays)
            else:
                viewer.add_overlay(overlays)
        
        viewer.show()
        return viewer

    def extractRepresentativeSlices(self, planes='all', offsets=[-10, 0, 10], verbose=False):
        """
        Extract representative 2D slices from 3 orthogonal planes around center-of-gravity.
        
        Useful for quick preview, batch processing, or input to vision models.
        
        Args:
            planes (str or list): Which planes to extract. Options:
                - 'all': All 3 planes (sagittal, coronal, axial) [default]
                - list of plane indices [0, 1, 2] or names ['sagittal', 'coronal', 'axial']
            offsets (list): Slice offsets from center-of-gravity (in mm). Default: [-10, 0, 10]
                            Produces 3 slices per plane × N planes
            verbose (bool): Print debug information
            
        Returns:
            dict: Contains:
                - 'slices': List of numpy arrays (2D slices in order: plane0_offset0, plane0_offset1, etc.)
                - 'plane_names': List of plane names for each slice group
                - 'offsets': The offsets used
                - 'center_of_gravity': Physical coordinates of center-of-gravity
                - 'center_of_gravity_index': Index coordinates of center-of-gravity
                
        Example:
            >>> img = Imaginable('mri_scan.nii.gz')
            >>> result = img.extractRepresentativeSlices(planes='all', offsets=[-5, 0, 5])
            >>> slices = result['slices']  # List of 9 numpy arrays (3 planes × 3 offsets)
            >>> for i, s in enumerate(slices):
            ...     print(f"Slice {i}: shape {s.shape}")
        """
        # Create a working copy oriented to LPS with isotropic spacing
        working = copy.deepcopy(self)
        working.dicomOrient('LPS')
        working.changeImageSpacing((1.0, 1.0, 1.0))
        
        # Get center of gravity or geometric center
        img = working.getImage()
        img_size = img.GetSize()
        
        # Try to compute center of gravity from binary mask
        try:
            # Cast to int32 for LabelShapeStatistics (doesn't support float64 in 3D)
            caster = sitk.CastImageFilter()
            caster.SetOutputPixelType(sitk.sitkInt32)
            img_int = caster.Execute(img)
            
            # Threshold to get foreground
            threshold_filter = sitk.BinaryThresholdImageFilter()
            threshold_filter.SetLowerThreshold(1)
            threshold_filter.SetUpperThreshold(255)
            img_binary = threshold_filter.Execute(img_int)
            
            stats = sitk.LabelShapeStatisticsImageFilter()
            stats.Execute(img_binary)
            
            # If label 1 exists, use its centroid
            if 1 in stats.GetLabels():
                center_of_gravity = stats.GetCentroid(1)
            else:
                # Fallback to geometric center
                center_of_gravity = tuple(s / 2.0 for s in img_size)
        except:
            # Fallback to geometric center if label statistics fails
            center_of_gravity = tuple(s / 2.0 for s in img_size)
        cog_index = working.getIndexFromCoordinates(center_of_gravity)
        
        if verbose:
            print(f"Center of gravity (physical): {center_of_gravity}")
            print(f"Center of gravity (index): {cog_index}")
            print(f"Image size: {working.getImageSize()}")
        
        # Determine which planes to extract
        plane_indices = []
        plane_names = ['sagittal', 'coronal', 'axial']
        
        if isinstance(planes, str):
            if planes == 'all':
                plane_indices = [0, 1, 2]
            else:
                raise ValueError(f"planes must be 'all' or a list of indices/names, got: {planes}")
        else:
            plane_indices = planes
        
        # Extract slices
        slices_list = []
        plane_labels = []
        
        for plane_idx in plane_indices:
            for offset in offsets:
                try:
                    slice_index = int(cog_index[plane_idx] + offset)
                    img_size = working.getImageSize()[plane_idx]
                    
                    # Check bounds
                    if slice_index < 0 or slice_index >= img_size:
                        if verbose:
                            print(f"Skipping {plane_names[plane_idx]} offset {offset}: "
                                  f"index {slice_index} out of bounds [0, {img_size})")
                        continue
                    
                    # Extract slice as numpy
                    try:
                        from .utils import getImaginableSliceNumpy
                    except ImportError:
                        from utils import getImaginableSliceNumpy
                    
                    slice_array = getImaginableSliceNumpy(working, plane_idx, slice_index)
                    
                    if slice_array.shape[0] >= 10 and slice_array.shape[1] >= 10:  # Minimum size
                        slices_list.append(slice_array)
                        plane_labels.append((plane_names[plane_idx], offset))
                        if verbose:
                            print(f"Extracted {plane_names[plane_idx]} offset {offset}: shape {slice_array.shape}")
                    else:
                        if verbose:
                            print(f"Skipped {plane_names[plane_idx]} offset {offset}: too small ({slice_array.shape})")
                        
                except Exception as e:
                    if verbose:
                        print(f"Error extracting {plane_names[plane_idx]} offset {offset}: {e}")
                    continue
        
        if len(slices_list) == 0:
            if verbose:
                print("Warning: No valid slices extracted!")
            return {
                'slices': [],
                'plane_names': [],
                'offsets': offsets,
                'center_of_gravity': center_of_gravity,
                'center_of_gravity_index': cog_index
            }
        
        return {
            'slices': slices_list,
            'plane_names': plane_labels,
            'offsets': offsets,
            'center_of_gravity': center_of_gravity,
            'center_of_gravity_index': cog_index
        }

    def renderIsosurface(self, isosurface_value=None, component_index=0, time_index=0, 
                        color=(1.0, 0.0, 0.0), opacity=1.0, show=True, title=None):
        """
        Render an isosurface of the image using VTK.
        
        For continuous images (Imaginable, SITKImaginable, Fieldable):
            Creates a 3D isosurface at the specified value.
        For vector fields (Fieldable):
            Can render magnitude or extract specific component/time.
        For ROIs (Roiable):
            Renders the boundary of the ROI at value 0.5 (between 0 and 1).
        
        Parameters
        ----------
        isosurface_value : float, optional
            The isovalue at which to create the surface. If None:
            - For continuous images: uses mean intensity
            - For ROIs: uses 0.5 (boundary between foreground/background)
            - For vector fields: uses magnitude mean
        component_index : int, default=0
            For multi-component images, which component to render.
            For vector fields, 0=magnitude, 1+=individual components.
        time_index : int, default=0
            For 4D images (3D + time), which time frame to render.
        color : tuple, default=(1.0, 0.0, 0.0)
            RGB color for the surface (0-1 range). Default: red.
        opacity : float, default=1.0
            Surface opacity (0-1). Default: fully opaque.
        show : bool, default=True
            If True, displays the isosurface in an interactive VTK window.
            If False, returns the actor without showing.
        title : str, optional
            Window title. Auto-generated if None.
        
        Returns
        -------
        vtk.vtkActor or tuple
            If show=False: returns (actor, renderer, window)
            If show=True: returns the actor after display
        
        Example
        -------
        **Continuous image:**
        
        >>> img = SITKImaginable('mri_scan.nii.gz')
        >>> # Render at 50% of intensity range
        >>> img.renderIsosurface(isosurface_value=100)
        
        **ROI/segmentation:**
        
        >>> roi = Roiable('segmentation.nii.gz')
        >>> # Render boundary (automatically uses 0.5)
        >>> roi.renderIsosurface(color=(0.0, 1.0, 0.0))
        
        **Vector field magnitude:**
        
        >>> vec = Fieldable('displacement_field.nii.gz')
        >>> # Render isosurface of displacement magnitude
        >>> vec.renderIsosurface(component_index=0, isosurface_value=5.0)
        
        **Without displaying (for batch processing):**
        
        >>> img = SITKImaginable('image.nii.gz')
        >>> actor, renderer, window = img.renderIsosurface(show=False)
        >>> # Manipulate actor/renderer/camera as needed
        
        Note
        ----
        Requires VTK to be installed. ROI images must be binary or will be thresholded.
        """
        import vtk
        
        # Get the image, handling multi-component/4D cases
        image = self.getImage()
        
        # Extract component or time frame if needed
        if image.GetNumberOfComponentsPerPixel() > 1 and component_index > 0:
            # Extract specific component from multi-component image
            extractor = sitk.VectorIndexSelectionCastImageFilter()
            extractor.SetIndex(component_index)
            image = extractor.Execute(image)
        
        if image.GetDimension() == 4:
            # Extract time frame from 4D image
            slicer = sitk.ExtractImageFilter()
            size = list(image.GetSize())
            size[3] = 0  # Remove time dimension
            slicer.SetSize(size)
            
            index = [0, 0, 0, time_index]
            slicer.SetIndex(index)
            image = slicer.Execute(image)
        
        # Convert to VTK
        try:
            from .meshable import sitk2vtk
        except ImportError:
            from meshable import sitk2vtk
        
        vtk_image = sitk2vtk(image)
        
        # Determine isosurface value if not provided
        if isosurface_value is None:
            # Compute statistics to find good default isosurface value
            stats_filter = sitk.StatisticsImageFilter()
            stats_filter.Execute(image)
            mean_val = stats_filter.GetMean()
            
            # For ROIs (Roiable), use 0.5 if not specified
            try:
                from .imaginable import Roiable
                if isinstance(self, Roiable):
                    isosurface_value = 0.5
                else:
                    isosurface_value = mean_val
            except:
                isosurface_value = mean_val
        
        # Create marching cubes isosurface
        mc_filter = vtk.vtkMarchingCubes()
        mc_filter.SetInputData(vtk_image)
        mc_filter.SetValue(0, isosurface_value)  # Set isovalue
        mc_filter.Update()
        
        polydata = mc_filter.GetOutput()
        
        # Create mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        mapper.ScalarVisibilityOff()  # Don't color by scalars, use actor color
        
        # Create actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*color)
        actor.GetProperty().SetOpacity(opacity)
        
        if not show:
            return (actor, None, None)  # Return actor without displaying
        
        # Create renderer and window
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)
        renderer.SetBackground(0.1, 0.1, 0.1)  # Dark gray background
        renderer.ResetCamera()
        
        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)
        render_window.SetSize(800, 600)
        
        if title is None:
            title = f"Isosurface (value={isosurface_value:.2f})"
        render_window.SetWindowName(title)
        
        # Add interactor
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(render_window)
        
        # Use trackball camera style for better interaction
        style = vtk.vtkInteractorStyleTrackballCamera()
        interactor.SetInteractorStyle(style)
        
        # Start interactive rendering
        interactor.Initialize()
        render_window.Render()
        interactor.Start()
        
        return actor

    # ========================================================================
    # SEGMENTATION: THRESHOLDING (on Imaginable)
    # ========================================================================

    def segmentOtsu(self, n_bins=128):
        """
        Segment using Otsu's automatic threshold.

        Parameters
        ----------
        n_bins : int
            Number of histogram bins for threshold computation.

        Returns
        -------
        Roiable
            Binary segmentation as a new Roiable.

        Example
        -------
        >>> roi = img.segmentOtsu()
        """
        from . import segmentation as seg
        result = seg.otsu_threshold(self.getImage(), n_bins=n_bins)
        r = Roiable()
        r.setImage(result, 'segmented via Otsu threshold')
        return r

    def segmentMultiOtsu(self, n_thresholds=2, n_bins=256):
        """
        Segment using multi-level Otsu thresholding.

        Parameters
        ----------
        n_thresholds : int
            Number of thresholds (produces n_thresholds + 1 classes).
        n_bins : int
            Number of histogram bins.

        Returns
        -------
        LabelMapable
            Multi-label segmentation.
        """
        from . import segmentation as seg
        result = seg.multi_otsu_threshold(self.getImage(), n_thresholds, n_bins)
        lm = LabelMapable()
        lm.setImage(result, f'segmented via multi-Otsu ({n_thresholds} thresholds)')
        return lm

    def segmentLi(self, n_bins=128):
        """
        Segment using Li's minimum cross-entropy threshold.

        Parameters
        ----------
        n_bins : int
            Number of histogram bins for threshold computation.

        Returns
        -------
        Roiable
            Binary segmentation.
        """
        from . import segmentation as seg
        result = seg.li_threshold(self.getImage(), n_bins=n_bins)
        r = Roiable()
        r.setImage(result, 'segmented via Li threshold')
        return r

    def segmentYen(self, n_bins=128):
        """
        Segment using Yen's entropy-based threshold.

        Parameters
        ----------
        n_bins : int
            Number of histogram bins for threshold computation.

        Returns
        -------
        Roiable
            Binary segmentation.
        """
        from . import segmentation as seg
        result = seg.yen_threshold(self.getImage(), n_bins=n_bins)
        r = Roiable()
        r.setImage(result, 'segmented via Yen threshold')
        return r

    def segmentTriangle(self, n_bins=128):
        """
        Segment using the triangle (Zack) threshold.

        Parameters
        ----------
        n_bins : int
            Number of histogram bins for threshold computation.

        Returns
        -------
        Roiable
            Binary segmentation.
        """
        from . import segmentation as seg
        result = seg.triangle_threshold(self.getImage(), n_bins=n_bins)
        r = Roiable()
        r.setImage(result, 'segmented via triangle threshold')
        return r

    def segmentHuang(self, n_bins=128):
        """
        Segment using Huang's fuzzy threshold.

        Parameters
        ----------
        n_bins : int
            Number of histogram bins for threshold computation.

        Returns
        -------
        Roiable
            Binary segmentation.
        """
        from . import segmentation as seg
        result = seg.huang_threshold(self.getImage(), n_bins=n_bins)
        r = Roiable()
        r.setImage(result, 'segmented via Huang threshold')
        return r

    def segmentThreshold(self, lower=0.0, upper=1.0):
        """
        Segment by applying a manual intensity threshold.

        Parameters
        ----------
        lower : float
            Lower intensity bound (inclusive).
        upper : float
            Upper intensity bound (inclusive).

        Returns
        -------
        Roiable
            Binary segmentation.
        """
        from . import segmentation as seg
        result = seg.manual_threshold(self.getImage(), lower, upper)
        r = Roiable()
        r.setImage(result, f'segmented via threshold [{lower}, {upper}]')
        return r

    def segmentConnectedThreshold(self, seed_roi, lower=None, upper=None,
                                    n_seeds=200, face_connected=True):
        """
        Region growing with explicit intensity bounds from a seed ROI.

        Parameters
        ----------
        seed_roi : Roiable or sitk.Image
            Binary seed region.
        lower : float, optional
            Lower intensity bound (auto if None).
        upper : float, optional
            Upper intensity bound (auto if None).
        n_seeds : int
            Maximum seed points.
        face_connected : bool
            If True, use 6-connectivity. If False, use 26-connectivity.

        Returns
        -------
        Roiable
            Grown region.
        """
        from . import segmentation as seg
        seed = seg._to_sitk(seed_roi)
        result = seg.connected_threshold_grow(
            self.getImage(), seed, lower=lower, upper=upper,
            n_seeds=n_seeds, face_connected=face_connected,
        )
        r = Roiable()
        r.setImage(result, 'segmented via connected threshold')
        return r

    def segmentNeighbourhoodConnected(self, seed_roi, lower=None, upper=None,
                                       radius=1, n_seeds=200):
        """
        Region growing with neighbourhood connectivity.

        Parameters
        ----------
        seed_roi : Roiable or sitk.Image
            Binary seed region.
        lower : float, optional
            Lower intensity bound (auto if None).
        upper : float, optional
            Upper intensity bound (auto if None).
        radius : int
            Neighbourhood radius.
        n_seeds : int
            Maximum seed points.

        Returns
        -------
        Roiable
            Grown region.
        """
        from . import segmentation as seg
        seed = seg._to_sitk(seed_roi)
        result = seg.neighbourhood_connected_grow(
            self.getImage(), seed, lower=lower, upper=upper,
            radius=radius, n_seeds=n_seeds,
        )
        r = Roiable()
        r.setImage(result, 'segmented via neighbourhood connected')
        return r

    def segmentIsolatedConnected(self, seed1_roi, seed2_roi, n_seeds=50):
        """
        Find the threshold separating two seed regions.

        Parameters
        ----------
        seed1_roi : Roiable or sitk.Image
            Target region seeds.
        seed2_roi : Roiable or sitk.Image
            Excluded region seeds.
        n_seeds : int
            Maximum seeds per region.

        Returns
        -------
        Roiable
            Segmentation of seed1 region.
        """
        from . import segmentation as seg
        s1 = seg._to_sitk(seed1_roi)
        s2 = seg._to_sitk(seed2_roi)
        result = seg.isolated_connected_grow(self.getImage(), s1, s2, n_seeds)
        r = Roiable()
        r.setImage(result, 'segmented via isolated connected')
        return r

    def segmentMorphologicalWatershed(self, level=0.1, fully_connected=False):
        """
        Morphological watershed segmentation (marker-free).

        Parameters
        ----------
        level : float
            Flooding level — higher produces fewer basins.
        fully_connected : bool
            Use 26-connectivity vs 6-connectivity.

        Returns
        -------
        LabelMapable
            Label image of watershed basins.
        """
        from . import segmentation as seg
        result = seg.morphological_watershed(self.getImage(), level, fully_connected)
        lm = LabelMapable()
        lm.setImage(result, f'morphological watershed (level={level})')
        return lm

    def segmentWatershedFromMarkers(self, markers, fully_connected=False):
        """
        Watershed segmentation driven by marker labels.

        Parameters
        ----------
        markers : LabelMapable or sitk.Image
            Integer marker image (each label seeds a basin).
        fully_connected : bool
            Use 26-connectivity vs 6-connectivity.

        Returns
        -------
        LabelMapable
            Label image of watershed basins.
        """
        from . import segmentation as seg
        mk = seg._to_sitk(markers)
        result = seg.morphological_watershed_from_markers(
            self.getImage(), mk, fully_connected,
        )
        lm = LabelMapable()
        lm.setImage(result, 'watershed from markers')
        return lm

    # ========================================================================
    # PREPROCESSING
    # ========================================================================

    def correctBiasField(self, mask=None, shrink_factor=4, n_iterations=None,
                         convergence_threshold=0.001, spline_order=3):
        """
        Apply N4 bias field correction (in-place).

        Corrects low-frequency intensity inhomogeneity (e.g., MRI coil bias).

        Parameters
        ----------
        mask : Roiable or sitk.Image, optional
            Mask for bias estimation (Otsu if None).
        shrink_factor : int
            Downsample factor for speed.
        n_iterations : list, optional
            Iterations per fitting level. Length controls number of
            fitting levels (default [50,50,50,50] = 4 levels).
        convergence_threshold : float
            Convergence threshold.
        spline_order : int
            B-spline order for bias field estimation (default 3).

        Returns
        -------
        self
        """
        from . import segmentation as seg
        mk = seg._to_sitk(mask) if mask is not None else None
        result = seg.n4_bias_field_correction(
            self.getImage(), mask=mk, shrink_factor=shrink_factor,
            n_iterations=n_iterations, convergence_threshold=convergence_threshold,
            spline_order=spline_order,
        )
        return self.setImage(result, 'N4 bias field corrected')

    def smoothAnisotropic(self, iterations=5, time_step=0.0625, conductance=3.0):
        """
        Apply curvature anisotropic diffusion smoothing (in-place).

        Smooths while preserving edges.

        Parameters
        ----------
        iterations : int
            Diffusion iterations.
        time_step : float
            Time step per iteration.
        conductance : float
            Conductance parameter.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        result = seg.anisotropic_diffusion(
            self.getImage(), iterations=iterations,
            time_step=time_step, conductance=conductance,
        )
        return self.setImage(result, 'anisotropic diffusion smoothed')

    def getEdgeMap(self, sigma=1.0):
        """
        Compute gradient-magnitude edge potential normalised to [0, 1].

        Parameters
        ----------
        sigma : float
            Gaussian sigma in mm.

        Returns
        -------
        Imaginable
            Edge map as a new Imaginable.
        """
        from . import segmentation as seg
        result = seg.compute_edge_map(self.getImage(), sigma)
        edge = Imaginable()
        edge.setImage(result, f'edge map (sigma={sigma})')
        return edge


def maskSITKImage(r,maskingvalue=1,foreground=1,outsidevalue=0):
    return sitk.Mask(r, sitk.Cast(foreground,sitk.sitkInt16), maskingValue=maskingvalue, outsideValue=outsidevalue)

    
    
    
    
def getDirectiontransform(image):
    dimension=image.getImageDimension()
    cosines = sitk.AffineTransform(dimension)
    cosines.SetCenter(image.getImageOrigin())
    return cosines


class SITKImaginable(Imaginable):
    """Thin ``Imaginable`` subclass kept for explicit scalar-image usage."""
    pass

class Roiable(Imaginable):
    """
    A class For Region of Interest
    if you have more than a label in your ROI, use Labelmapable
    ROI can have a value but when getimage is called the output has 1 for region in ROI and 0 for outside
    
    """

    def __init__(self, filename=None, image=None, verbose=False,roivalue=None):
        """
        filename : str
            the filename of the ROI
        image : sitk.image
            the simple itk image that sotres the roi
        roivalue : int
            roi value
        """
        super().__init__(filename, image, verbose)
        self.dfltInterpolator=sitk.sitkNearestNeighbor
        self.dfltuseNearestNeighborExtrapolator=True
        if roivalue:
            self.setImage(self.getImage()==roivalue,'now the mask is equal to one')

    def getCenterOfGravityCoordinates(self):
        """
        Get the center of gravity of the ROI
        Returns:
            _type_: _description_
        """
        label_image = sitk.Cast(self.getImage(), sitk.sitkInt32)
        feature_image = self.getImage()
        label_statistic = sitk.LabelIntensityStatisticsImageFilter()
        label_statistic.Execute(label_image, feature_image)
        center_gravity = label_statistic.GetCenterOfGravity(1)
        return center_gravity

    def getCentroidCoordinates(self):
        """
        Get the Centroid of the ROI
        Returns:
            _type_: _description_
        """
        label_image = sitk.Cast(self.getImage(), sitk.sitkInt32)
        feature_image = self.getImage()
        label_statistic = sitk.LabelIntensityStatisticsImageFilter()
        label_statistic.Execute(label_image, feature_image)
        center_gravity = label_statistic.GetCentroid(1)
        return center_gravity
    
    def  getCenterOfGravityIndex(self):
        """
        Get the center of gravity of the ROI
        Returns:
            _type_: _description_
        """
        center = self.getIndexFromCoordinates(self.getCenterOfGravityCoordinates())
        return center

    def getCentroidIndex(self):
        """
        Get the Centroid of the ROI
        """
        Centroid_coordinate = self.getIndexFromCoordinates(self.getCentroidCoordinates())
        return Centroid_coordinate
    def dilateRadius(self,radius=2):
        return self.__derodeRadius__(radius,False)
    def erodeRadius(self,radius=2):
        return self.__derodeRadius__(radius)
    def __derodeRadius__(self,radius=2,erode=True):
        image=self.getImage()

        if erode:
            filter = sitk.BinaryErodeImageFilter()
        else:
            filter = sitk.BinaryDilateImageFilter()
        filter.SetKernelRadius ( radius )
        o=filter.Execute(image>0)
        self.setImage(o,'erode')
        return self
    
    def removeSmallObj(self,voxel_threshold=50,connectivity=26):
        mask =measure.label(self.getImageAsNumpy()==1)
        mask = morphology.remove_small_objects(mask, voxel_threshold,connectivity=connectivity)
        mask[np.where(mask>0)]=1
        self.setImageFromNumpy(mask,refimage=self.getImage())
        return self

    def removeHoles(self,voxel_threshold=50,connectivity=26):
        mask =measure.label(self.getImageAsNumpy()==1)
        mask = morphology.remove_small_holes(mask, voxel_threshold,connectivity=connectivity)
        mask[np.where(mask>0)]=1
        self.setImageFromNumpy(mask,refimage=self.getImage())
        return self
    def keepBiggestObj(self,connectivity=26):
        labelled = measure.label(self.getImageAsNumpy()==1)
        rp = measure.regionprops(labelled)
        # get size of largest cluster
        size = max([i.area for i in rp])
        # remove everything smaller than largest
        mask = morphology.remove_small_objects(labelled, min_size=size-1,connectivity=connectivity)
        mask[np.where(mask>0)]=1
        self.setImageFromNumpy(mask)
        return self

    # ========================================================================
    # ROI-SPECIFIC DEFORMATION METHODS
    # ========================================================================

    def warpROI(self, displacement_field, target_image=None):
        """
        Apply a displacement field to warp this ROI/mask, preserving label values.

        Uses nearest-neighbor interpolation to maintain ROI integrity across all labels.

        Parameters
        ----------
        displacement_field : str or sitk.Image
            Path to displacement field file (.mha, .nii.gz) or SimpleITK vector image
        target_image : str or sitk.Image, optional
            Target geometry reference

        Returns
        -------
        self : Roiable
            Self for method chaining

        Example
        -------
        >>> roi = Roiable('segmentation.nii.gz')
        >>> roi.warpROI('deformation.mha', target_image='fixed.nii.gz')
        >>> roi.write('warped_roi.nii.gz')
        """
        from . import deformations
        
        warped = deformations.apply_deformation_field_to_labels(
            self.getImage(),
            displacement_field,
            target_image=target_image,
            default_label=0
        )
        
        return self.setImage(warped, f"applied displacement field to ROI from {displacement_field if isinstance(displacement_field, str) else 'field object'}")

    def applyTransformToROI(self, transform, target_image=None):
        """
        Backward-compatible alias for ``applyTransform`` on ROI data.

        Parameters
        ----------
        transform : str or sitk.Transform
            Transform file or transform object.
        target_image : str or sitk.Image, optional
            Optional reference geometry for the warped ROI.

        Returns
        -------
        Roiable
            Self for method chaining.
        """
        return self.applyTransform(transform, target_image=target_image)

    # ========================================================================
    # SEGMENTATION REFINEMENT METHODS
    # ========================================================================

    def refineWatershed(self, image, height_map=None, erosion_iters=3,
                        dilation_iters=5, min_voxels=50):
        """
        Refine ROI boundaries using marker-based watershed segmentation.

        Markers are created from an eroded interior (foreground) and a
        dilated exterior (background).

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Reference intensity image (e.g. MRI scan).
        height_map : Imaginable or sitk.Image, optional
            Custom height / cost map.  If *None*, the gradient magnitude
            of *image* is used.
        erosion_iters : int
            Erosion iterations for foreground markers (default 3).
        dilation_iters : int
            Dilation iterations for background markers (default 5).
        min_voxels : int
            Remove components smaller than this (default 50).

        Returns
        -------
        self
            For method chaining.

        Example
        -------
        >>> roi = Roiable('mask.nii.gz')
        >>> roi.refineWatershed(Imaginable('scan.nii.gz'))
        """
        from . import segmentation as seg
        img = seg._to_sitk(image)
        hm = seg._to_sitk(height_map) if height_map is not None else None
        result = seg.watershed_refine(
            self.getImage(), img, height_map=hm,
            erosion_iters=erosion_iters, dilation_iters=dilation_iters,
            min_voxels=min_voxels,
        )
        return self.setImage(result, 'refined via watershed')

    
    
    
    def getParaViewSurface(
        self,
        space="lps",
        smooth_iterations=0,
        relaxation_factor=0.1,
    ):
        """
        Return the ROI as a coordinate-aware vtkPolyData surface.
        """
        from .meshable import binary_mask_to_polydata

        return binary_mask_to_polydata(
            self.getImage(),
            space=space,
            smooth_iterations=smooth_iterations,
            relaxation_factor=relaxation_factor,
        )


    def writeParaView(
        self,
        output_path,
        mode="surface",
        space="lps",
        stride=1,
        array_name="mask",
        smooth_iterations=0,
        relaxation_factor=0.1,
    ):
        """
        Export the ROI for ParaView.

        mode="surface"
            Write a .vtp surface.

        mode="volume"
            Write a .vts voxel volume.
        """
        from .meshable import write_vtk_dataset

        mode = mode.lower()

        if mode == "volume":
            return super().writeParaView(
                output_path,
                space=space,
                stride=stride,
                array_name=array_name,
            )

        if mode != "surface":
            raise ValueError(
                "Roiable mode must be 'surface' or 'volume'."
            )

        surface = self.getParaViewSurface(
            space=space,
            smooth_iterations=smooth_iterations,
            relaxation_factor=relaxation_factor,
        )

        return write_vtk_dataset(
            surface,
            output_path,
        )


    def refineRegionGrowing(self, image, multiplier=2.5, neighborhood_radius=1,
                            n_iterations=3, max_distance_mm=10.0, n_seeds=100,
                            prob_map=None, prob_threshold=0.1, min_voxels=50):
        """
        Refine ROI via confidence-connected region growing.

        Seeds are sampled from the eroded ROI interior and the region
        grows into intensity-similar voxels.

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Reference intensity image.
        multiplier : float
            Standard-deviation multiplier for intensity acceptance range.
        neighborhood_radius : int
            Radius for local statistics.
        n_iterations : int
            Region-growing iterations.
        max_distance_mm : float
            Maximum growth distance from original ROI surface (mm).
        n_seeds : int
            Maximum number of seed points.
        prob_map : Imaginable or sitk.Image, optional
            Probability map — grown region is intersected with values
            above *prob_threshold*.
        prob_threshold : float
            Threshold on *prob_map*.
        min_voxels : int
            Remove components smaller than this.

        Returns
        -------
        self

        Example
        -------
        >>> roi.refineRegionGrowing(scan, multiplier=2.0, max_distance_mm=5.0)
        """
        from . import segmentation as seg
        img = seg._to_sitk(image)
        pm = seg._to_sitk(prob_map) if prob_map is not None else None
        result = seg.region_growing_refine(
            self.getImage(), img, multiplier=multiplier,
            neighborhood_radius=neighborhood_radius, n_iterations=n_iterations,
            max_distance_mm=max_distance_mm, n_seeds=n_seeds,
            prob_map=pm, prob_threshold=prob_threshold, min_voxels=min_voxels,
        )
        return self.setImage(result, 'refined via region growing')

    def refineGeodesicActiveContour(self, image, propagation=0.5, curvature=0.3,
                                    advection=1.0, iterations=50,
                                    rms_tolerance=0.001, allow_shrink=False,
                                    speed_image=None):
        """
        Refine ROI boundaries using a geodesic active contour level-set.

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Reference intensity image for edge computation.
        propagation : float
            Balloon force (positive → expand, negative → shrink).
        curvature : float
            Smoothing force (higher = smoother boundaries).
        advection : float
            Edge-attraction force.
        iterations : int
            Maximum GAC iterations.
        rms_tolerance : float
            Convergence threshold.
        allow_shrink : bool
            If *False*, the result is the union with the seed
            (expansion-only).
        speed_image : Imaginable or sitk.Image, optional
            Pre-computed speed / feature image.

        Returns
        -------
        self

        Example
        -------
        >>> roi.refineGeodesicActiveContour(scan, propagation=0.3, curvature=0.8)
        """
        from . import segmentation as seg
        img = seg._to_sitk(image)
        si = seg._to_sitk(speed_image) if speed_image is not None else None
        result = seg.geodesic_active_contour_refine(
            self.getImage(), img, propagation=propagation,
            curvature=curvature, advection=advection, iterations=iterations,
            rms_tolerance=rms_tolerance, allow_shrink=allow_shrink,
            speed_image=si,
        )
        return self.setImage(result, 'refined via geodesic active contour')

    def expandByProbability(self, prob_map, threshold=0.25, max_layers=5,
                            max_distance_mm=None, fill_holes=True,
                            min_voxels=50):
        """
        Expand ROI layer-by-layer, accepting only voxels above a
        probability threshold.

        Parameters
        ----------
        prob_map : Imaginable or sitk.Image
            Probability / confidence map (float [0, 1]).
        threshold : float
            Minimum probability for a voxel to be included.
        max_layers : int
            Maximum expansion layers.
        max_distance_mm : float, optional
            Hard distance cap from original ROI surface.
        fill_holes : bool
            Fill holes after expansion.
        min_voxels : int
            Remove components smaller than this.

        Returns
        -------
        self

        Example
        -------
        >>> roi.expandByProbability(tissue_prob, threshold=0.3, max_layers=3)
        """
        from . import segmentation as seg
        pm = seg._to_sitk(prob_map)
        result = seg.probability_expansion(
            self.getImage(), pm, threshold=threshold, max_layers=max_layers,
            max_distance_mm=max_distance_mm, fill_holes=fill_holes,
            min_voxels=min_voxels,
        )
        return self.setImage(result, 'expanded by probability')

    def shrinkByProbability(self, prob_map, threshold=0.15, max_layers=5,
                            min_preserve_fraction=0.5, fill_holes=True):
        """
        Shrink ROI by removing low-probability boundary voxels.

        Parameters
        ----------
        prob_map : Imaginable or sitk.Image
            Probability / confidence map (float [0, 1]).
        threshold : float
            Boundary voxels below this probability are removed.
        max_layers : int
            Maximum shrinkage layers.
        min_preserve_fraction : float
            Stop if volume falls below this fraction of the original.
        fill_holes : bool
            Fill holes after shrinkage.

        Returns
        -------
        self

        Example
        -------
        >>> roi.shrinkByProbability(tissue_prob, threshold=0.2, max_layers=3)
        """
        from . import segmentation as seg
        pm = seg._to_sitk(prob_map)
        result = seg.probability_shrinkage(
            self.getImage(), pm, threshold=threshold, max_layers=max_layers,
            min_preserve_fraction=min_preserve_fraction, fill_holes=fill_holes,
        )
        return self.setImage(result, 'shrunk by probability')

    # ========================================================================
    # LEVEL-SET REFINEMENTS
    # ========================================================================

    def refineThresholdLevelSet(self, image, lower_threshold=0.1,
                                upper_threshold=0.9, propagation=1.0,
                                curvature=1.0, iterations=100,
                                rms_tolerance=0.02, allow_shrink=True):
        """
        Refine ROI using a threshold-based level-set.

        The contour expands into voxels whose normalised intensity
        falls within [lower_threshold, upper_threshold].

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Reference intensity image.
        lower_threshold : float
            Lower intensity bound (normalised [0, 1]).
        upper_threshold : float
            Upper intensity bound (normalised [0, 1]).
        propagation : float
            Balloon force.
        curvature : float
            Smoothing force.
        iterations : int
            Maximum iterations.
        rms_tolerance : float
            Convergence threshold.
        allow_shrink : bool
            If False, result is unioned with seed.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        img = seg._to_sitk(image)
        result = seg.threshold_level_set_refine(
            self.getImage(), img,
            lower_threshold=lower_threshold, upper_threshold=upper_threshold,
            propagation=propagation, curvature=curvature,
            iterations=iterations, rms_tolerance=rms_tolerance,
            allow_shrink=allow_shrink,
        )
        return self.setImage(result, 'refined via threshold level-set')

    def refineLaplacianLevelSet(self, image, propagation=1.0, curvature=1.0,
                                iterations=100, rms_tolerance=0.02,
                                allow_shrink=True):
        """
        Refine ROI using a Laplacian-based level-set.

        Drives the contour towards zero-crossings of the Laplacian.

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Reference intensity image.
        propagation : float
            Balloon force.
        curvature : float
            Smoothing force.
        iterations : int
            Maximum iterations.
        rms_tolerance : float
            Convergence threshold.
        allow_shrink : bool
            If False, result is unioned with seed.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        img = seg._to_sitk(image)
        result = seg.laplacian_level_set_refine(
            self.getImage(), img,
            propagation=propagation, curvature=curvature,
            iterations=iterations, rms_tolerance=rms_tolerance,
            allow_shrink=allow_shrink,
        )
        return self.setImage(result, 'refined via Laplacian level-set')

    def refineShapeDetectionLevelSet(self, image, propagation=1.0, curvature=0.5,
                                     iterations=100, rms_tolerance=0.02,
                                     sigma_mm=1.0, allow_shrink=True):
        """
        Refine ROI using a shape detection level-set.

        Similar to geodesic active contour but without advection.

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Reference intensity image.
        propagation : float
            Balloon force.
        curvature : float
            Smoothing force.
        iterations : int
            Maximum iterations.
        rms_tolerance : float
            Convergence threshold.
        sigma_mm : float
            Gaussian sigma for edge computation.
        allow_shrink : bool
            If False, result is unioned with seed.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        img = seg._to_sitk(image)
        result = seg.shape_detection_level_set_refine(
            self.getImage(), img,
            propagation=propagation, curvature=curvature,
            iterations=iterations, rms_tolerance=rms_tolerance,
            sigma_mm=sigma_mm, allow_shrink=allow_shrink,
        )
        return self.setImage(result, 'refined via shape detection level-set')

    def refineChanVese(self, image, lambda1=1.0, lambda2=1.0,
                       curvature_weight=0.0, area_weight=0.0,
                       volume_weight=0.0, iterations=100,
                       rms_tolerance=0.02, allow_shrink=True):
        """
        Refine ROI using Chan-Vese (region-based) level-set.

        Works well for images with weak or absent edges.

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Reference intensity image.
        lambda1 : float
            Inside-region variance weight.
        lambda2 : float
            Outside-region variance weight.
        curvature_weight : float
            Curvature regularisation.
        area_weight : float
            Area penalty.
        volume_weight : float
            Volume penalty.
        iterations : int
            Maximum iterations.
        rms_tolerance : float
            Convergence threshold.
        allow_shrink : bool
            If False, result is unioned with seed.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        img = seg._to_sitk(image)
        result = seg.chan_vese_refine(
            self.getImage(), img,
            lambda1=lambda1, lambda2=lambda2,
            curvature_weight=curvature_weight, area_weight=area_weight,
            volume_weight=volume_weight, iterations=iterations,
            rms_tolerance=rms_tolerance, allow_shrink=allow_shrink,
        )
        return self.setImage(result, 'refined via Chan-Vese')

    # ========================================================================
    # DISTANCE AND MASK OPERATIONS
    # ========================================================================

    def constrainByDistance(self, reference_roi, max_distance_mm=5.0,
                           exclude_interior=False):
        """
        Clip ROI to stay within a distance from a reference ROI.

        Parameters
        ----------
        reference_roi : Roiable or sitk.Image
            Reference ROI.
        max_distance_mm : float
            Maximum allowed distance from reference surface.
        exclude_interior : bool
            If True, also remove voxels inside the reference.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        ref = seg._to_sitk(reference_roi)
        result = seg.constrain_by_distance(
            self.getImage(), ref,
            max_distance_mm=max_distance_mm,
            exclude_interior=exclude_interior,
        )
        return self.setImage(result, f'constrained by distance ({max_distance_mm}mm)')

    def subtractMask(self, mask_to_remove):
        """
        Remove voxels that overlap with another mask.

        Parameters
        ----------
        mask_to_remove : Roiable or sitk.Image
            Binary mask of voxels to remove.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        mk = seg._to_sitk(mask_to_remove)
        result = seg.subtract_mask(self.getImage(), mk)
        return self.setImage(result, 'subtracted mask')

    def intersectWith(self, other):
        """
        Keep only voxels present in both this ROI and other.

        Parameters
        ----------
        other : Roiable or sitk.Image
            Second binary ROI.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        o = seg._to_sitk(other)
        result = seg.intersect_masks(self.getImage(), o)
        return self.setImage(result, 'intersected with mask')

    def unionWith(self, other):
        """
        Combine this ROI with another (logical OR).

        Parameters
        ----------
        other : Roiable or sitk.Image
            Second binary ROI.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        o = seg._to_sitk(other)
        result = seg.union_masks(self.getImage(), o)
        return self.setImage(result, 'unioned with mask')

    # ========================================================================
    # MORPHOLOGICAL OPERATIONS (mm-based)
    # ========================================================================

    def erodeMM(self, radius_mm=1.0):
        """
        Erode the ROI by a radius in mm.

        Parameters
        ----------
        radius_mm : float
            Erosion radius in mm.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        result = seg.binary_erode(self.getImage(), radius_mm)
        return self.setImage(result, f'eroded {radius_mm}mm')

    def dilateMM(self, radius_mm=1.0):
        """
        Dilate the ROI by a radius in mm.

        Parameters
        ----------
        radius_mm : float
            Dilation radius in mm.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        result = seg.binary_dilate(self.getImage(), radius_mm)
        return self.setImage(result, f'dilated {radius_mm}mm')

    def openMM(self, radius_mm=1.0):
        """
        Morphological opening (erosion then dilation) in mm.

        Removes small protrusions and disconnected fragments.

        Parameters
        ----------
        radius_mm : float
            Structuring element radius in mm.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        result = seg.binary_open(self.getImage(), radius_mm)
        return self.setImage(result, f'opened {radius_mm}mm')

    def closeMM(self, radius_mm=1.0):
        """
        Morphological closing (dilation then erosion) in mm.

        Fills small holes and gaps.

        Parameters
        ----------
        radius_mm : float
            Structuring element radius in mm.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        result = seg.binary_close(self.getImage(), radius_mm)
        return self.setImage(result, f'closed {radius_mm}mm')

    def fillBinaryHoles(self):
        """
        Fill all enclosed holes inside the ROI.

        Unlike ``removeHoles`` (which targets small holes), this fills
        *every* internal cavity regardless of size.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        result = seg.fill_binary_holes(self.getImage())
        return self.setImage(result, 'binary holes filled')

    def getDistanceMap(self):
        """
        Compute the Euclidean distance transform from the ROI surface.

        Returns
        -------
        Imaginable
            Float image where each voxel contains the distance (mm) to the
            nearest ROI surface voxel.  Zero inside the ROI.
        """
        from . import segmentation as seg
        dist = seg.compute_distance_map(self.getImage())
        result = Imaginable(image=dist)
        return result

    def getShell(self, width_mm=5.0):
        """
        Compute a shell (annular band) around the ROI surface.

        Parameters
        ----------
        width_mm : float
            Shell width in mm (extends outward from ROI surface).

        Returns
        -------
        Roiable
            Binary mask of the shell region.
        """
        from . import segmentation as seg
        shell = seg.compute_shell(self.getImage(), width_mm=width_mm)
        return Roiable(image=shell)

    @staticmethod
    def fuseSTAPLE(segmentations, confidence_threshold=0.5, min_voxels=50):
        """
        Fuse multiple binary segmentations using the STAPLE algorithm.

        Parameters
        ----------
        segmentations : list of Roiable or sitk.Image
            Binary segmentations to fuse.
        confidence_threshold : float
            Threshold on STAPLE probability for the final output.
        min_voxels : int
            Remove components smaller than this.

        Returns
        -------
        Roiable
            Fused binary ROI.

        Example
        -------
        >>> fused = Roiable.fuseSTAPLE([roi1, roi2, roi3])
        """
        from . import segmentation as seg
        result = seg.staple_fusion(
            segmentations, confidence_threshold=confidence_threshold,
            min_voxels=min_voxels,
        )
        return Roiable(image=result)

    # ========================================================================
    # COMPARISON & SURFACE METRICS
    # ========================================================================

    def compareTo(self, other):
        """
        Full comparison against another ROI: overlap + surface distance metrics.

        Parameters
        ----------
        other : Roiable or sitk.Image
            Reference ROI to compare against.

        Returns
        -------
        dict
            Dice, Jaccard, VolumeSimilarity, FalseDiscoveryRate,
            FalseNegativeError, FalsePositiveError, MeanSurfaceDist_mm,
            RMSDist_mm, HD95_mm, SurfaceDice_1mm.

        Example
        -------
        >>> metrics = predicted_roi.compareTo(ground_truth_roi)
        >>> print(f"Dice = {metrics['Dice']:.3f}")
        """
        from . import metrics as met
        ref = getmeTheSimpleITKImage(other)
        return met.compare_segmentations(self.getImage(), ref)

    def getSurfaceDistances(self, other):
        """
        Symmetric surface distance metrics against another ROI.

        Parameters
        ----------
        other : Roiable or sitk.Image
            Reference ROI.

        Returns
        -------
        dict
            MeanSurfaceDist_mm, RMSDist_mm, HD95_mm, SurfaceDice_1mm.
        """
        from . import metrics as met
        ref = getmeTheSimpleITKImage(other)
        return met.compute_surface_distances(self.getImage(), ref)

    def exportSurfaceWithError(self, other, out_dir='/tmp/surface_error',
                               highlight_percentile=95):
        """
        Export this ROI and a reference as PLY/VTP meshes with per-vertex error.

        Requires ``pyvista`` and ``scikit-image``.

        Parameters
        ----------
        other : Roiable or sitk.Image
            Reference (ground-truth) ROI.
        out_dir : str
            Output directory for mesh files.
        highlight_percentile : int
            Percentile threshold for the high-error sub-mesh.

        Example
        -------
        >>> pred.exportSurfaceWithError(gt, out_dir='/g/paraview')
        """
        from . import metrics as met
        ref_img = getmeTheSimpleITKImage(other)
        met.export_surface_with_error(
            reference_arr=ref_img,
            test_arr=self.getImage(),
            spacing=self.getImageSpacing(),
            out_dir=out_dir,
            highlight_percentile=highlight_percentile,
        )

    # ========================================================================
    # QUALITY SCORES (no ground truth needed)
    # ========================================================================

    def getEdgeAlignmentScore(self, image, band_width=2, sigma_mm=1.0):
        """
        How well the ROI boundary aligns with image edges.

        Higher score means better alignment with actual anatomical boundaries.

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Grayscale image (e.g. MRI).
        band_width : int
            Surface band width in voxels.
        sigma_mm : float
            Gaussian sigma for gradient computation.

        Returns
        -------
        float
            Edge alignment score in [0, 1].
        """
        from . import metrics as met
        img = getmeTheSimpleITKImage(image)
        return met.compute_edge_alignment_score(
            self.getImage(), img, band_width=band_width, sigma_mm=sigma_mm)

    def getCompactnessScore(self):
        """
        Volume / surface-area ratio.  Higher = more compact / smoother.

        Returns
        -------
        float
        """
        from . import metrics as met
        return met.compute_compactness_score(self.getImage())

    def getConnectedComponentCount(self):
        """
        Number of connected components.  Ideally 1.

        Returns
        -------
        int
        """
        from . import metrics as met
        return met.compute_connectivity_score(self.getImage())

    # ========================================================================
    # MORPHOMETRIC ANALYSIS
    # ========================================================================

    def getMorphometrics(self):
        """
        Comprehensive morphometric analysis of the ROI.

        Returns a dict containing PCA-based extents (length, width, thickness),
        max Feret diameter, max inscribed-sphere thickness, and local
        thickness statistics.

        Returns
        -------
        dict
            Keys: principal (length/width/thickness), max_feret_mm,
            max_inscribed_thickness_mm, thickness_stats (mean/median/std/max/p25/p75).

        Example
        -------
        >>> roi = Roiable('cartilage.nii.gz')
        >>> m = roi.getMorphometrics()
        >>> print(f"Length={m['principal']['length']:.1f} mm")
        """
        from . import metrics as met
        mask = self.getImageAsNumpy() > 0
        spacing = self.getImageSpacing()  # (x, y, z)
        # metrics expect ZYX spacing
        sp_zyx = tuple(reversed(spacing))

        principal = met.roi_principal_extents(mask, sp_zyx)
        result = {
            'principal': {
                'length': principal['length'],
                'width': principal['width'],
                'thickness': principal['thickness'],
            },
            'max_feret_mm': met.roi_max_feret(mask, sp_zyx),
            'max_inscribed_thickness_mm': met.roi_max_thickness_inscribed(mask, sp_zyx),
            'thickness_stats': met.roi_mean_thickness(mask, sp_zyx),
        }
        return result

    def getPrincipalExtents(self):
        """
        PCA-based oriented bounding-box extents (length, width, thickness).

        Returns
        -------
        dict with length, width, thickness (all in mm).
        """
        from . import metrics as met
        mask = self.getImageAsNumpy() > 0
        sp_zyx = tuple(reversed(self.getImageSpacing()))
        return met.roi_principal_extents(mask, sp_zyx)

    def getMaxFeret(self):
        """
        Maximum Feret diameter in mm (largest pairwise distance on
        the convex hull).

        Returns
        -------
        float
        """
        from . import metrics as met
        mask = self.getImageAsNumpy() > 0
        sp_zyx = tuple(reversed(self.getImageSpacing()))
        return met.roi_max_feret(mask, sp_zyx)

    def getMaxInscribedThickness(self):
        """
        Maximum inscribed-sphere thickness: ``2 * max(EDT)`` in mm.

        Returns
        -------
        float
        """
        from . import metrics as met
        mask = self.getImageAsNumpy() > 0
        sp_zyx = tuple(reversed(self.getImageSpacing()))
        return met.roi_max_thickness_inscribed(mask, sp_zyx)

    def getThicknessStats(self):
        """
        Local thickness statistics via EDT.

        Returns
        -------
        dict with mean, median, std, max, p25, p75 (all in mm).
        """
        from . import metrics as met
        mask = self.getImageAsNumpy() > 0
        sp_zyx = tuple(reversed(self.getImageSpacing()))
        return met.roi_mean_thickness(mask, sp_zyx)

    # ========================================================================
    # BINARY SMOOTHING
    # ========================================================================

    def smoothBinary(self, radius=1, mode='closing'):
        """
        Morphological smoothing of the binary ROI.

        Parameters
        ----------
        radius : int
            Kernel radius in voxels.
        mode : str
            ``'closing'`` (fill small gaps), ``'opening'`` (remove small bumps),
            or ``'both'`` (closing then opening).

        Returns
        -------
        self
        """
        img = self.getImage()
        if mode in ('closing', 'both'):
            img = sitk.BinaryMorphologicalClosing(img > 0, [radius] * img.GetDimension())
        if mode in ('opening', 'both'):
            img = sitk.BinaryMorphologicalOpening(img > 0, [radius] * img.GetDimension())
        return self.setImage(sitk.Cast(img, sitk.sitkUInt8), f'smoothBinary mode={mode} r={radius}')

    # ========================================================================
    # DISTANCE MAP & SURFACE METHODS
    # ========================================================================

    def getSurfaceMask(self):
        """
        Extract 1-voxel-thick surface contour of the ROI.

        Returns
        -------
        Roiable
            Binary surface mask.

        Example
        -------
        >>> surface = roi.getSurfaceMask()
        >>> surface.write('surface.nii.gz')
        """
        from . import metrics as met
        surf = met.surface_mask(self.getImage())
        return Roiable(image=surf)

    def getSignedDistanceMap(self):
        """
        Signed Maurer distance map from this ROI.

        Negative inside, positive outside, in mm (physical spacing).

        Returns
        -------
        Imaginable
            Float32 signed distance map.

        Example
        -------
        >>> sdm = roi.getSignedDistanceMap()
        >>> sdm.write('signed_dist.nii.gz')
        """
        from . import metrics as met
        sdm = met.signed_distance_map(self.getImage())
        return Imaginable(image=sdm)

    def getSurfaceDistanceMap(self, other):
        """
        Image where voxels on THIS surface contain distance (mm) to OTHER mask.

        All non-surface voxels are 0.

        Parameters
        ----------
        other : Roiable or sitk.Image
            Target mask to measure distance to.

        Returns
        -------
        Imaginable
            Float32 distance map (non-zero only on this ROI's surface).

        Example
        -------
        >>> dist_img = pre_roi.getSurfaceDistanceMap(post_roi)
        >>> dist_img.write('surface_dist.nii.gz')
        """
        from . import metrics as met
        other_img = getmeTheSimpleITKImage(other)
        dm = met.surface_distance_map_image(self.getImage(), other_img)
        return Imaginable(image=dm)

    def getChangeMaps(self, other):
        """
        Boolean change maps between this ROI (pre) and another (post).

        Parameters
        ----------
        other : Roiable or sitk.Image
            Post-registration mask (should be on the same grid).

        Returns
        -------
        dict
            ``{'removed': Roiable, 'added': Roiable, 'changed': Roiable}``

        Example
        -------
        >>> changes = pre_roi.getChangeMaps(post_roi_registered)
        >>> changes['removed'].write('removed.nii.gz')
        """
        from . import metrics as met
        other_img = getmeTheSimpleITKImage(other)
        maps = met.compute_change_maps(self.getImage(), other_img)
        return {k: Roiable(image=v) for k, v in maps.items()}

    def splitByConnectedComponents(self, sort_by='size'):
        """
        Split this ROI into its connected components.

        Parameters
        ----------
        sort_by : str
            ``'size'`` (descending voxel count) or ``'none'``.

        Returns
        -------
        list[Roiable]
            One Roiable per connected component.

        Example
        -------
        >>> parts = roi.splitByConnectedComponents()
        >>> largest = parts[0]
        """
        from . import metrics as met
        masks = met.split_connected_components(self.getImage(), sort_by=sort_by)
        return [Roiable(image=m) for m in masks]

    def splitLeftRight(self, axis=0):
        """
        Split into left and right parts using connected-component centroids.

        The two largest components are identified.  "Left" is the one with the
        smaller centroid coordinate along *axis* (default: X = left-right in LPS).

        Parameters
        ----------
        axis : int
            Physical axis index (0=X, 1=Y, 2=Z).

        Returns
        -------
        (left, right) : tuple[Roiable, Roiable]

        Example
        -------
        >>> left, right = femur_roi.splitLeftRight()
        """
        from . import metrics as met
        left, right = met.split_left_right(self.getImage(), axis=axis)
        return Roiable(image=left), Roiable(image=right)

    def exportSurfaceVTP(self, out_path, distance_to=None, scalar_name='dist_mm'):
        """
        Export marching-cubes surface mesh as VTP, optionally coloured by
        distance to another mask.

        Requires VTK (``pip install vtk``).

        Parameters
        ----------
        out_path : str
            Output ``.vtp`` file path.
        distance_to : Roiable or sitk.Image, optional
            If provided, each surface vertex is coloured by distance (mm)
            to this mask.
        scalar_name : str
            Name of the scalar array in the VTP file.

        Example
        -------
        >>> roi.exportSurfaceVTP('surface.vtp')
        >>> roi.exportSurfaceVTP('error.vtp', distance_to=other_roi)
        """
        from . import metrics as met
        dist_img = None
        if distance_to is not None:
            other_img = getmeTheSimpleITKImage(distance_to)
            dist_img = sitk.Abs(sitk.SignedMaurerDistanceMap(
                sitk.Cast(other_img > 0, sitk.sitkUInt8),
                insideIsPositive=False,
                squaredDistance=False,
                useImageSpacing=True,
            ))
        met.export_surface_vtp(
            self.getImage(), out_path,
            distance_image=dist_img,
            scalar_name=scalar_name,
        )

    def describe(self):
        """Print a concise summary of the ROI.

        Returns:
            dict: key ROI properties
        """
        info = super().describe() if hasattr(super(), 'describe') else {}
        try:
            info['roi_value'] = getattr(self, 'roiValue', None)
            arr = self.getImageAsNumpy()
            if arr is not None:
                info['unique_values'] = sorted(int(v) for v in np.unique(arr))
                info['non_zero_voxels'] = int(np.count_nonzero(arr))
        except Exception as e:
            info['roi_error'] = str(e)
        self._print_describe(info)
        return info



class LabelMapable(Imaginable):
    """
    Multi-label segmentation wrapper with label-preserving transforms.

    This class extends ``Imaginable`` with per-label extraction, priors,
    refinement, longitudinal comparison, and registration helpers.
    """
    def __init__(self, filename=None, image=None, verbose=False):
        super().__init__(filename, image, verbose)
        self.dfltInterpolator=sitk.sitkNearestNeighbor
        self.dfltuseNearestNeighborExtrapolator=True
    def noNearestNeighborExtrapolator(self):
        self.dfltuseNearestNeighborExtrapolator=False
        return self

    @classmethod
    def fromMasks(cls, masks, values, reference=None, remove_small_obj=None, overlap="overwrite"):
        """
        Create a LabelMapable from N binary/ROI masks and N label values.

        Args:
            masks: list of mask sources (file paths, Roiable, or Imaginable objects)
            values: list of integer label values, one per mask
            reference: reference image for resampling (sitk.Image or Imaginable);
                       defaults to the first mask's image
            remove_small_obj: if set, remove connected components smaller than this
                              voxel count from each mask before assigning
            overlap: how to handle overlapping masks
                - "overwrite": later masks replace earlier labels (default)
                - "keep_first": earlier labels are preserved
                - "error": raise ValueError if any masks overlap

        Returns:
            LabelMapable with labels assigned from the masks
        """
        import numpy as np

        if len(masks) != len(values):
            raise ValueError("masks and values must have the same length")

        if reference is None:
            first = Roiable(masks[0]) if isinstance(masks[0], str) else masks[0]
            reference = first.getImage()

        combined = None
        occupied = None

        for mask_src, value in zip(masks, values):
            roi = Roiable(mask_src) if isinstance(mask_src, str) else mask_src
            roi.resampleOnTargetImage(reference)

            mask = roi.getImageAsNumpy() > 0

            if remove_small_obj:
                tmp = Roiable(image=roi.getImage())
                tmp.setImageFromNumpy(mask.astype(np.uint8), refimage=reference)
                tmp.removeSmallObj(remove_small_obj)
                mask = tmp.getImageAsNumpy() > 0

            if combined is None:
                combined = np.zeros(mask.shape, dtype=np.uint16)
                occupied = np.zeros(mask.shape, dtype=bool)

            overlap_mask = mask & occupied
            if overlap == "error" and np.any(overlap_mask):
                raise ValueError(f"Mask overlap detected for label {value}")

            if overlap == "keep_first":
                mask = mask & ~occupied

            combined[mask] = int(value)
            occupied |= mask

        out = cls()
        out.setImageFromNumpy(combined, refimage=reference)
        out.cast("uint16")
        return out

    def getCenterOfGravityCoordinatesPerLabel(self):
        """
        Get the center of gravity of the labelmap per label
        Returns:
            _type_: _description_
        """
        label_image = sitk.Cast(self.getImage(), sitk.sitkInt32)
        feature_image = self.getImage()
        label_statistic = sitk.LabelIntensityStatisticsImageFilter()
        label_statistic.Execute(label_image, feature_image)
        
        centers_gravity = {}
        for label in label_statistic.GetLabels():
            centers_gravity[label] = label_statistic.GetCenterOfGravity(label)
        
        return centers_gravity
    def getCenterOfGravityCoordinates(self):
        """
        Get the center of gravity of the labelmap
        Returns:
            _type_: _description_
        """
        centers_gravity = self.getCenterOfGravityCoordinatesPerLabel()
        
        center_of_all = np.mean(list(centers_gravity.values()), axis=0)

        return center_of_all
 
    def getCentroidCoordinatesPerLabel(self):
        """
        Get the Centroid of the labelmap per label
        Returns:
            _type_: _description_
        """
        label_image = sitk.Cast(self.getImage(), sitk.sitkInt32)
        feature_image = self.getImage()
        label_statistic = sitk.LabelIntensityStatisticsImageFilter()
        label_statistic.Execute(label_image, feature_image)
        
        centers = {}
        for label in label_statistic.GetLabels():
            centers[label] = label_statistic.GetCentroid(label)
        
        return centers
    
    def getCenterOfGravityIndex(self):
        center = self.getIndexFromCoordinates(self.getCenterOfGravityCoordinates())
        return center
    
    def getCentroidCoordinates(self):
        """ 
        Get the Centroid of the labelmap

        Returns:
            _type_: _description_
        """
        centers = self.getCentroidCoordinatesPerLabel()
        center_of_all = np.mean(list(centers.values()), axis=0)
        return center_of_all
    
    def getCentroidIndex(self):
        Centroid = self.getIndexFromCoordinates(self.getCentroidCoordinates())
        return Centroid

    def describe(self):
        """Print a concise summary of the label map.

        Returns:
            dict: key label map properties
        """
        info = super().describe() if hasattr(super(), 'describe') else {}
        try:
            arr = self.getImageAsNumpy()
            if arr is not None:
                labels = sorted(int(v) for v in np.unique(arr))
                info['labels'] = labels
                info['num_labels'] = len(labels)
                info['non_zero_voxels'] = int(np.count_nonzero(arr))
        except Exception as e:
            info['labelmap_error'] = str(e)
        self._print_describe(info)
        return info

    # ========================================================================
    # PER-LABEL EXTRACTION AND MANIPULATION
    # ========================================================================

    def getLabels(self, exclude_background=True):
        """
        Return the list of unique label values in the label map.

        Parameters
        ----------
        exclude_background : bool
            If *True*, label 0 is excluded.

        Returns
        -------
        list of int
        """
        arr = self.getImageAsNumpy()
        labels = sorted(int(v) for v in np.unique(arr))
        if exclude_background and 0 in labels:
            labels.remove(0)
        return labels

    def extractLabel(self, label_value):
        """
        Extract a single label as a binary Roiable.

        Parameters
        ----------
        label_value : int
            The label value to extract.

        Returns
        -------
        Roiable
            Binary mask where the label equals *label_value*.

        Example
        -------
        >>> bone = labelmap.extractLabel(1)
        >>> bone.describe()
        """
        arr = (self.getImageAsNumpy() == label_value).astype(np.uint8)
        roi = Roiable()
        roi.setImageFromNumpy(arr, refimage=self.getImage())
        return roi

    def setLabel(self, label_value, roi):
        """
        Set / overwrite a single label from a binary Roiable.

        Existing voxels with *label_value* are first cleared, then the
        non-zero voxels of *roi* are written with *label_value*.

        Parameters
        ----------
        label_value : int
            Label value to write.
        roi : Roiable or sitk.Image
            Binary mask indicating where *label_value* should be placed.

        Returns
        -------
        self
        """
        from . import segmentation as seg
        roi_sitk = seg._to_sitk(roi)
        roi_sitk = seg._ensure_same_grid(roi_sitk, self.getImage(),
                                          sitk.sitkNearestNeighbor)
        arr = self.getImageAsNumpy().copy()
        roi_arr = sitk.GetArrayFromImage(roi_sitk) > 0
        # Clear existing label
        arr[arr == label_value] = 0
        # Write new label
        arr[roi_arr] = label_value
        self.setImageFromNumpy(arr, refimage=self.getImage())
        return self

    # ========================================================================
    # COMPARISON METHODS
    # ========================================================================

    def writeParaView(
        self,
        output_path,
        space="lps",
        labels=None,
        prefix="label",
        label_names=None,
        smooth_iterations=0,
        relaxation_factor=0.1,
    ):
        """
        Export a label map for ParaView.

        If output_path ends with ".vtm", all labels are written into one
        multiblock VTK file.

        Otherwise, output_path is treated as a directory and one ".vtp"
        surface file is written for each nonzero label.

        Parameters
        ----------
        output_path : str or pathlib.Path
            Output directory for separate VTP files, or a VTM filename.
        space : str
            Coordinate space: "lps", "ras", "fsl", or "index".
        labels : sequence of int, optional
            Labels to export. By default, all nonzero labels are exported.
        prefix : str
            Filename prefix for separate label files.
        label_names : dict, optional
            Mapping from label values to readable names.
        smooth_iterations : int
            Optional number of surface smoothing iterations.
        relaxation_factor : float
            Smoothing relaxation factor.

        Returns
        -------
        dict or str
            Dictionary of label-to-file mappings for separate files,
            or the VTM filename when writing a multiblock file.
        """
        from pathlib import Path
        import re

        import numpy as np
        import SimpleITK as sitk

        from .meshable import (
            label_map_to_multiblock,
            write_vtk_dataset,
        )

        output_path = Path(output_path)
        image = self.getImage()

        array = sitk.GetArrayViewFromImage(image)

        if labels is None:
            labels = [
                int(value)
                for value in np.unique(array)
                if int(value) != 0
            ]
        else:
            labels = [
                int(value)
                for value in labels
                if int(value) != 0
            ]

        if not labels:
            raise ValueError(
                "The LabelMapable contains no nonzero labels to export."
            )

        # --------------------------------------------------------------
        # Option 1: one VTM multiblock file
        # --------------------------------------------------------------
        if output_path.suffix.lower() == ".vtm":
            surfaces = label_map_to_multiblock(
                image,
                labels=labels,
                space=space,
                smooth_iterations=smooth_iterations,
                relaxation_factor=relaxation_factor,
            )

            return write_vtk_dataset(
                surfaces,
                output_path,
            )

        # Reject unsupported file extensions.
        if output_path.suffix:
            raise ValueError(
                "For a LabelMapable, output_path must either be a "
                "directory or a filename ending in '.vtm'."
            )

        # --------------------------------------------------------------
        # Option 2: one VTP file per label
        # --------------------------------------------------------------
        output_path.mkdir(
            parents=True,
            exist_ok=True,
        )

        label_names = label_names or {}
        written_files = {}

        for label_value in labels:
            mask = sitk.Cast(
                image == label_value,
                sitk.sitkUInt8,
            )

            if not np.any(
                sitk.GetArrayViewFromImage(mask)
            ):
                continue

            # Reuse the Roiable implementation.
            roi = Roiable(image=mask)

            label_name = label_names.get(
                label_value,
                f"label_{label_value}",
            )

            safe_name = re.sub(
                r"[^A-Za-z0-9_.-]+",
                "_",
                str(label_name),
            ).strip("_")

            output_filename = (
                output_path
                / f"{prefix}_{safe_name}.vtp"
            )

            written_files[label_value] = roi.writeParaView(
                output_filename,
                space=space,
                smooth_iterations=smooth_iterations,
                relaxation_factor=relaxation_factor,
            )

        return written_files

    def compareToByLabel(self, other, labels=None):
        """
        Per-label comparison against another label map.

        Returns overlap and surface distance metrics for each shared label.

        Parameters
        ----------
        other : LabelMapable or sitk.Image
            Reference label map.
        labels : list of int, optional
            Labels to compare.  If *None*, uses all non-zero labels from
            both label maps.

        Returns
        -------
        dict[int, dict]
            Mapping ``label_value → metrics_dict``.

        Example
        -------
        >>> results = pred_lm.compareToByLabel(gt_lm)
        >>> for label, m in results.items():
        ...     print(f"Label {label}: Dice={m['Dice']:.3f}")
        """
        from . import metrics as met

        other_img = getmeTheSimpleITKImage(other)
        other_img = met._ensure_same_geometry(self.getImage(), other_img)

        if labels is None:
            arr_self = self.getImageAsNumpy()
            arr_other = sitk.GetArrayFromImage(other_img)
            labels = sorted(set(
                int(v) for v in np.unique(arr_self) if v != 0
            ) | set(
                int(v) for v in np.unique(arr_other) if v != 0
            ))

        results = {}
        for lbl in labels:
            self_bin = sitk.Cast(sitk.Equal(self.getImage(), int(lbl)), sitk.sitkUInt8)
            other_bin = sitk.Cast(sitk.Equal(other_img, int(lbl)), sitk.sitkUInt8)
            # Skip if both empty
            s_sum = sitk.GetArrayFromImage(self_bin).sum()
            o_sum = sitk.GetArrayFromImage(other_bin).sum()
            if s_sum == 0 and o_sum == 0:
                continue
            results[lbl] = met.compare_segmentations(self_bin, other_bin)
        return results

    # ========================================================================
    # LABEL-MAP PRIORS
    # ========================================================================

    def buildPriors(self, tau=0.8, blur_sigma_mm=0.6, classes=None):
        """
        Build soft probability priors from this label map.

        For each foreground class a sigmoid of the signed distance map
        is computed and optionally blurred.  Returns a vector Imaginable
        whose channels are ``[background, class_1, class_2, …]``.

        Parameters
        ----------
        tau : float
            Sigmoid steepness (smaller → sharper boundaries).
        blur_sigma_mm : float
            Gaussian anti-aliasing blur in mm (0 = no blur).
        classes : list of int, optional
            Foreground labels to include.  If *None*, auto-detected.

        Returns
        -------
        prior_image : Imaginable
            Vector (multi-channel) image with probability priors.
        class_list : list[int]
            Ordered label list matching channels.

        Example
        -------
        >>> priors, cls = labelmap.buildPriors(tau=0.8)
        >>> priors.write('priors.nii.gz')
        """
        from . import metrics as met
        prior_vec, class_list = met.build_label_priors(
            self.getImage(), classes=classes, tau=tau,
            blur_sigma_mm=blur_sigma_mm,
        )
        result = Imaginable(image=prior_vec)
        return result, class_list

    @staticmethod
    def combineBinaryMasks(mask_list, priority_order=None):
        """
        Combine multiple binary masks into a single multi-label image.

        Parameters
        ----------
        mask_list : list of Roiable or sitk.Image
            Binary masks (same geometry).  >0 = foreground.
        priority_order : list of int, optional
            Overwrite order (indices into *mask_list*).  Last wins.

        Returns
        -------
        LabelMapable
            Multi-label image (labels 1..K).

        Example
        -------
        >>> lm = LabelMapable.combineBinaryMasks([bone_roi, cartilage_roi])
        """
        from . import metrics as met
        imgs = [getmeTheSimpleITKImage(m) for m in mask_list]
        lab_img = met.combine_binary_masks_to_label(imgs, priority_order)
        lm = LabelMapable(image=lab_img)
        return lm

    # ========================================================================
    # SEGMENTATION REFINEMENT METHODS
    # ========================================================================

    def refineLabel(self, label_value, image, method='watershed', **kwargs):
        """
        Refine a single label using the specified segmentation method.

        The label is extracted as a Roiable, refined, and written back.

        Parameters
        ----------
        label_value : int
            Label to refine.
        image : Imaginable or sitk.Image
            Reference intensity image.
        method : str
            One of ``'watershed'``, ``'region_growing'``, ``'gac'``,
            ``'expand'``, ``'shrink'``.
        **kwargs
            Extra parameters forwarded to the refinement method.

        Returns
        -------
        self

        Example
        -------
        >>> labelmap.refineLabel(1, scan, method='gac', propagation=0.3)
        """
        roi = self.extractLabel(label_value)
        method_map = {
            'watershed': roi.refineWatershed,
            'region_growing': roi.refineRegionGrowing,
            'gac': roi.refineGeodesicActiveContour,
            'expand': roi.expandByProbability,
            'shrink': roi.shrinkByProbability,
        }
        fn = method_map.get(method)
        if fn is None:
            raise ValueError(
                f"Unknown method '{method}'. "
                f"Choose from: {list(method_map.keys())}"
            )
        fn(image, **kwargs)
        self.setLabel(label_value, roi)
        return self

    def refineAllLabels(self, image, method='watershed',
                        resolve_overlaps=True, **kwargs):
        """
        Refine every non-zero label and optionally resolve overlaps.

        Parameters
        ----------
        image : Imaginable or sitk.Image
            Reference intensity image.
        method : str
            Refinement method (see ``refineLabel``).
        resolve_overlaps : bool
            If *True*, overlapping voxels between labels are assigned
            to the nearest label centroid after refinement.
        **kwargs
            Extra parameters forwarded to the refinement method.

        Returns
        -------
        self

        Example
        -------
        >>> labelmap.refineAllLabels(scan, method='gac', propagation=0.3)
        """
        labels = self.getLabels(exclude_background=True)
        for lbl in labels:
            self.refineLabel(lbl, image, method=method, **kwargs)
        if resolve_overlaps and len(labels) > 1:
            self.resolveOverlaps()
        return self

    def resolveOverlaps(self):
        """
        Resolve overlapping voxels between labels.

        Overlapping voxels are assigned to the label whose centroid
        is closest (Euclidean distance in mm).

        Returns
        -------
        self

        Example
        -------
        >>> labelmap.resolveOverlaps()
        """
        from . import segmentation as seg
        labels = self.getLabels(exclude_background=True)
        if len(labels) < 2:
            return self

        # Extract per-label binary ROIs
        label_rois = {}
        for lbl in labels:
            label_rois[lbl] = self.extractLabel(lbl).getImage()

        # Resolve
        resolved = seg.resolve_label_overlaps(label_rois, reference=self.getImage())

        # Rebuild label map
        arr = np.zeros_like(self.getImageAsNumpy(), dtype=self.getImageAsNumpy().dtype)
        for lbl in labels:
            roi_arr = sitk.GetArrayFromImage(resolved[lbl]) > 0
            arr[roi_arr] = lbl

        self.setImageFromNumpy(arr, refimage=self.getImage())
        return self

    # ========================================================================
    # MULTI-LABEL DEFORMATION METHODS
    # ========================================================================

    def warpLabelMap(self, displacement_field, target_image=None):
        """
        Apply a displacement field to warp this label map, preserving all label values.

        Uses nearest-neighbor interpolation to maintain label integrity across all labels.

        Parameters
        ----------
        displacement_field : str or sitk.Image
            Path to displacement field file (.mha, .nii.gz) or SimpleITK vector image
        target_image : str or sitk.Image, optional
            Target geometry reference

        Returns
        -------
        self : LabelMapable
            Self for method chaining

        Example
        -------
        >>> labels = LabelMapable('segmentation.nii.gz')
        >>> labels.warpLabelMap('deformation.mha', target_image='fixed.nii.gz')
        >>> labels.write('warped_labels.nii.gz')
        """
        from . import deformations
        
        warped = deformations.apply_deformation_field_to_labels(
            self.getImage(),
            displacement_field,
            target_image=target_image,
            default_label=0
        )
        
        return self.setImage(warped, f"applied displacement field to label map from {displacement_field if isinstance(displacement_field, str) else 'field object'}")

    def applyTransformToLabelMap(self, transform, target_image=None):
        """
        Backward-compatible alias for ``applyTransform`` on label maps.

        Parameters
        ----------
        transform : str or sitk.Transform
            Transform file or transform object.
        target_image : str or sitk.Image, optional
            Optional reference geometry for the warped label map.

        Returns
        -------
        LabelMapable
            Self for method chaining.
        """
        return self.applyTransform(transform, target_image=target_image)

    # ========================================================================
    # REGISTRATION & LONGITUDINAL ANALYSIS
    # ========================================================================

    def registerTo(self, other, method='rigid', roi_values=None, iterations=200):
        """
        Register this label map onto *other* using signed-distance maps.

        A union mask of all (or specified) labels is used for registration.
        The result is this label map resampled into *other*'s space.

        Parameters
        ----------
        other : LabelMapable or sitk.Image
            Fixed (target) label map.
        method : str
            ``'rigid'`` (Euler3D) or ``'affine'``.
        roi_values : list[int], optional
            Labels to include in the union mask for registration.
            If *None*, all non-zero labels are used.
        iterations : int
            Max optimiser iterations.

        Returns
        -------
        (registered, transform) : tuple[LabelMapable, sitk.Transform]
            ``registered`` is this label map warped into *other*'s space.

        Example
        -------
        >>> post_reg, T = post_labels.registerTo(pre_labels, method='rigid')
        """
        from . import metrics as met
        other_img = getmeTheSimpleITKImage(other)
        self_img = self.getImage()

        # Build union masks
        if roi_values is None:
            self_arr = sitk.GetArrayFromImage(self_img)
            roi_values = sorted(int(v) for v in np.unique(self_arr) if v != 0)

        def _union(img, vals):
            m = sitk.Equal(img, vals[0])
            for v in vals[1:]:
                m = m | sitk.Equal(img, v)
            return sitk.Cast(m, sitk.sitkUInt8)

        fixed_mask = _union(other_img, roi_values)
        moving_mask = _union(self_img, roi_values)

        # Register on signed distance maps
        _, transform = met.rigid_register_masks(
            fixed_mask, moving_mask,
            method=method, iterations=iterations,
        )

        # Resample entire label map with nearest-neighbor
        resampled = sitk.Resample(
            sitk.Cast(self_img, sitk.sitkUInt16),
            other_img,
            transform,
            sitk.sitkNearestNeighbor,
            0,
            sitk.sitkUInt16,
        )
        result = LabelMapable(image=resampled)
        return result, transform

    def getChangeMapsByLabel(self, other, labels=None):
        """
        Per-label change maps between this label map (pre) and *other* (post).

        Parameters
        ----------
        other : LabelMapable or sitk.Image
            Post-registration label map (should be on the same grid).
        labels : list[int], optional
            Labels to analyse.  If *None*, all non-zero labels from both
            images are used.

        Returns
        -------
        dict[int, dict[str, Roiable]]
            ``{label: {'removed': Roiable, 'added': Roiable, 'changed': Roiable}}``

        Example
        -------
        >>> changes = pre.getChangeMapsByLabel(post_registered)
        >>> changes[1]['removed'].write('roi1_removed.nii.gz')
        """
        from . import metrics as met
        other_img = getmeTheSimpleITKImage(other)
        other_img = met._ensure_same_geometry(self.getImage(), other_img)

        if labels is None:
            pre_arr = sitk.GetArrayFromImage(self.getImage())
            post_arr = sitk.GetArrayFromImage(other_img)
            labels = sorted(
                set(int(v) for v in np.unique(pre_arr) if v != 0)
                | set(int(v) for v in np.unique(post_arr) if v != 0)
            )

        result = {}
        for lbl in labels:
            pre_bin = sitk.Cast(sitk.Equal(self.getImage(), lbl), sitk.sitkUInt8)
            post_bin = sitk.Cast(sitk.Equal(other_img, lbl), sitk.sitkUInt8)
            maps = met.compute_change_maps(pre_bin, post_bin)
            result[lbl] = {k: Roiable(image=v) for k, v in maps.items()}
        return result

    def toRoiable(self, exclude_background=True):
        """
        Convert this label map into a list of Roiable objects, one per label.

        Returns
        -------
        list[Roiable]
            List of binary Roiable objects for each label.

        Example
        -------
        >>> rois = labelmap.toRoiable()
        >>> for roi in rois:
        ...     roi.describe()
        """
        rois = []
        for lbl in self.getLabels(exclude_background=exclude_background):
            rois.append(self.extractLabel(lbl))
        return rois

class LabelMapableROI(LabelMapable):
    """Old class for Labelmapable
    will be removed in the future
    use LabelMapable instead
    Args:
        LabelMapable (_type_): _description_
    """
    def __init__(self, filename=None, image=None, verbose=False,labelsvalues=None):
        super().__init__(filename, image, verbose)
        self.dfltInterpolator=sitk.sitkNearestNeighbor
        self.dfltuseNearestNeighborExtrapolator=True
        self.labelsvalues=labelsvalues
        self.ROIS=[]
        if self.labelsvalues is None:
            self.labelsvalues=self.getImageUniqueValues(exclude=[0])

        for v in self.labelsvalues:
            self.ROIS.append(Roiable(filename,image,verbose,roivalue=v))
        

    def removeSmallObj(self,voxel_threshold=50,connectivity=26):
        for r in self.ROIS:
            r.removeSmallObj(voxel_threshold,connectivity)
        self.mergeLabels()
        return self

    def removeHoles(self,voxel_threshold=50,connectivity=26):
        for r in self.ROIS:
            r.removeHoles(voxel_threshold,connectivity)
        # self.mergeLabels()
        return self
    def keepBiggestObj(self,connectivity=26):
        for r in self.ROIS:
            r.keepBiggestObj(connectivity)
        # self.mergeLabels()
        return self

    def mergeLabels(self):
        LABELMAP = None
        for rd,v in zip(self.ROIS,self.labelsvalues):
            O=rd.getImageAsNumpy()
            try:
                if LABELMAP is None:
                    LABELMAP = np.zeros_like(O, dtype=np.float32)
                LABELMAP[np.where(O==1)]=v    
            except Exception as e:
                print(f"Error in mergeLabels: {e}")
                raise
        self.setImageFromNumpy(LABELMAP,refimage=super().getImage())
        return self

class Fieldable(Imaginable):
    """Vector/displacement-field flavored ``Imaginable``."""
    def __init__(self, filename=None, image=None, verbose=False):
        super().__init__(filename, image, verbose)
        self.dfltInterpolator=sitk.sitkNearestNeighbor
        self.dfltuseNearestNeighborExtrapolator=False

    def setImageFromNumpy(self, nparray, refimage=None, spacing=None, origin=None, direction=None):
        """
        Set vector field from numpy array in (Z, Y, X, components) order.
        
        For vector fields, the last dimension should be the vector components.
        For 3D: (Z, Y, X, 3) for a 3-component vector field
        
        Args:
            nparray: Numpy array in (Z, Y, X, components) order
            refimage: Reference image to copy metadata from
            spacing: Physical spacing if not using refimage
            origin: Physical origin if not using refimage
            direction: Direction matrix if not using refimage
        
        Returns:
            self
        
        Note: Changed in v3! Now expects (Z,Y,X,components) to match getImageAsNumpy().
        """
        vector = True
        # Array is already in (Z,Y,X,components), pass directly to sitk.GetImageFromArray
        nda = sitk.GetImageFromArray(nparray, isVector=vector)
        if refimage:
            REF = getmeTheSimpleITKImage(refimage)
            # For vector images, compare shape without components dimension
            if len(nparray.shape) > 3:
                img_shape = nparray.shape[:3]  # (Z,Y,X)
            else:
                img_shape = nparray.shape
            if np.array_equiv(REF.GetSize(), img_shape[::-1]):  # ITK size is (X,Y,Z)
                nda.CopyInformation(REF)
            else:
                nda = setSITKImageInforFromImage(nda, REF)
        elif ((spacing) and (origin) and (direction)):
            nda = setSITKImageInfo(nda, spacing=spacing, origin=origin, direction=direction)
        elif self.isImageSet():
            if(self.getImage()):
                r, o, d = getSITKImageInfo(getmeTheSimpleITKImage(self))
                nda = setSITKImageInfo(nda, spacing=r, origin=o, direction=d)            
        self.setImage(nda, 'vector field set from numpy array (Z,Y,X,components)!')
        return self

    def setImageFromNumpyZYX(self, nparray, refimage=None, spacing=None, origin=None, direction=None):
        """
        Alias for setImageFromNumpy() for explicit clarity.
        """
        return self.setImageFromNumpy(nparray, refimage, spacing, origin, direction)
    
    def setImageFromNumpyXYZ(self, nparray, refimage=None, spacing=None, origin=None, direction=None):
        """
        Set vector field from numpy array in (X, Y, Z, components) order.
        
        DEPRECATED: Provided for backward compatibility only.
        """
        # Transpose from (X,Y,Z,components) to (Z,Y,X,components)
        if len(nparray.shape) == 4:  # 3D vector field
            o = np.transpose(nparray, (2, 1, 0, 3))
        elif len(nparray.shape) == 3:  # Could be 2D vector or 3D scalar
            # Assume it's (X,Y,components) for 2D vector
            o = np.transpose(nparray, (1, 0, 2))
        else:
            L = list(range(len(nparray.shape)))
            L.reverse()
            o = np.transpose(nparray, L)
        return self.setImageFromNumpy(o, refimage, spacing, origin, direction)

    def describe(self):
        """Print a concise summary of the vector field.

        Returns:
            dict: key field properties
        """
        info = super().describe() if hasattr(super(), 'describe') else {}
        try:
            img = self.getImage()
            if img is not None:
                info['num_components'] = img.GetNumberOfComponentsPerPixel()
        except Exception as e:
            info['field_error'] = str(e)
        self._print_describe(info)
        return info


#     def toVtk(self):
#         return sitk2vtk(self.getImage(), debugOn=False)
if __name__=="__main__":
    A=Imaginable('/data/MYDATA/CARTILAGE_HIP/images/wo.nii')
    A.padImage([10,20,50],[10,10,50])
    A.writeImageAs('/g/a.nii')
    
    
    


    # A=Imaginable('/data/MYDATA/fulldixon-images/C-4/data/IN.nii')
    # B=Imaginable('/data/MYDATA/fulldixon-images/C-4/data/OUT.nii')
    # A.add(B)
    # A.divide(0.5)
    # A.writeImageAs('/g/_wo.nii')
    # A.reset()
    # A.subtract(B)
    # A.divide(0.5)
    # A.applyAbs()
    # A.writeImageAs('/g/_fo.nii')
