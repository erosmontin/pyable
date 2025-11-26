from imaginable import Imaginable
import numpy as np
import matplotlib.pyplot as plt
A=Imaginable('/media/erosm/DATA/aging20/DATA/HCP/HCA6002236_V1_MR/aparc+aseg.nii.gz')

N=A.getImageAsNumpy()

plt.imshow(N[90,:,:])
plt.show()

A.reorientToRAS()
N_ras=A.getImageAsNumpy()
plt.imshow(N_ras[90,:,:])
plt.title('RAS')
plt.show()


A.undo()
A.dicomOrient('RAS')
N=A.getImageAsNumpy()
plt.imshow(N[90,:,:])
plt.title('RAS')
plt.show()


A.undo()
A.reorientToRPI()
N_ras=A.getImageAsNumpy()
plt.imshow(N_ras[90,:,:])
plt.title('RAS')
plt.show()


A.undo()
A.dicomOrient('RPI')
N=A.getImageAsNumpy()
plt.imshow(N[90,:,:])
plt.title('RPI')
plt.show()

A.undo()
A.dicomOrient('LPS')
N=A.getImageAsNumpy()
plt.imshow(N[90,:,:])
plt.title('LPS')
plt.show()

A.undo()
A.dicomOrient('LPI')
N=A.getImageAsNumpy()
plt.imshow(N[90,:,:])
plt.title('LPI')
plt.show()