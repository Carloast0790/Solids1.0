#!/home/carlos/miniconda3/bin/python3.8 -u
#w.array.py name.xyz 3 1
import numpy as np
from utils.libmoleculas import writexyzs, translate, translate_to_cm, uniontab, align_all_inertia_axis_x, molecular_radius_max
#------------------------------------------------------------------------------------------
def array_mol(moleculein,xmax,ymax):
    mol01=moleculein[0:xmax*ymax]
    #mol01=align_all_inertia_axis_x(mol01)
    rmax=max([molecular_radius_max(imol) for imol in mol01])
    fixdistance=rmax+8.0
    fiydistance=rmax+8.0
    imol=0
    yd=float(0.0)
    for yi in range(ymax):
        yd=yd-fiydistance
        xd=float(0.0)
        for xi in range(xmax):
            xd=xd+fixdistance   
            vector=np.array([xd, yd, 0.0])
            translate_to_cm(mol01[imol])        
            translate(mol01[imol],vector)
            imol=imol+1
    mol02=uniontab(mol01)
    writexyzs([mol02],'array.xyz')
#------------------------------------------------------------------------------------------
import sys
from utils.libmoleculas import readxyzs
name=sys.argv[1]
xmax=int(sys.argv[2])
ymax=int(sys.argv[3])
molx=readxyzs(name)
array_mol(molx,xmax,ymax)
exit()
