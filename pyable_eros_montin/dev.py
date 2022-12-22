from imaginable import Imaginable,numpyToImaginable
from utils import wlt
import numpy as np

class Imaginable(Imaginable):




from pynico_eros_montin import pynico as pn

if __name__=="__main__":
    # FN='/g/LeftRight111/input/Leftp03-1.nii.gz'
    FN='/g/LeftRight111/input/Leftp02.nii.gz'
    K=Imaginable(filename=FN)
    WT=K.getWavelet('db1')
    for t in WT:
        t["map"].writeImageAs(pn.Pathable(FN).addPrefix(f'new{t["name"]}').changePath('/g').getPosition())
        

