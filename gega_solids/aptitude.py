import numpy as np
from utils_solids.libmoleculas import sort_by_energy
#------------------------------------------------------------------------------------------
def get_aptitude(moleculein):
    if len(moleculein)==1:
       aptitude=[float(1.0)]
    else:
       aptitude=[]
       moleculetmp=sort_by_energy(moleculein,1)
       Emin=moleculetmp[0].e
       Emax=moleculetmp[-1].e
       for imol in moleculetmp:
           Ei=imol.e
           EFEi=(Ei-Emin)/(Emax-Emin)
           fi=0.5*(1.0-np.tanh(((2.0*EFEi)-1.0)))
           aptitude.append(float(fi))
    return aptitude
#------------------------------------------------------------------------------------------