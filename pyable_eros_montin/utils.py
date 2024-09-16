try:
    import imaginable as ima
except:
    import pyable_eros_montin.imaginable as ima

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


from PIL import Image

import matplotlib
def saveSliceToImage(I,axis,index,fn,spacing=None):
    if spacing:
        I.changeImageSpacing(spacing)
    f=getImaginableSlice(I,axis=axis,index=index)
    matplotlib.image.imsave(fn, f,vmin=0.2,vmax=1,cmap='jet'  )


import matplotlib.pyplot as plt
import matplotlib.cm as cm

def overlayNumpyImageAndNumpyLabelmap(image, labelmap, image_cmap='gray', labelmap_cmap='jet', alpha_value=0.5, image_vmin=None, image_vmax=None, labelmap_vmin=None, labelmap_vmax=None,show=False,save=None,title=None):
    # Display the image as it is
    plt.imshow(image, cmap=image_cmap, vmin=image_vmin, vmax=image_vmax,origin='lower')

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
    labelmap_colored = cm.get_cmap(labelmap_cmap)(labelmap_norm(labelmap))


    # Create an alpha channel based on the labelmap
    alpha_channel = np.where(labelmap == 0, 0, alpha_value)

    # Replace the alpha channel in the colored labelmap
    labelmap_colored[..., 3] = alpha_channel

    # Overlay the labelmap on top of the image
    lbl=plt.imshow(labelmap_colored,origin='lower')
    # Create a ScalarMappable object for the colorbar
    sm = plt.cm.ScalarMappable(cmap=labelmap_cmap, norm=labelmap_norm)
    sm.set_array([])

    # Add the colorbar
    plt.colorbar(sm, label='Labelmap')

        # Invert x and y axes
    if title:
        plt.title(title)
    
    if save:
        plt.savefig(save,dpi=300)
    if show:
        plt.show()
    
    



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
    IM.overlayAble(R,2,40,show=True)
    
    
    

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
