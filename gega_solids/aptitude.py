import numpy as np
from utils_solids.libmoleculas import sort_by_energy
#------------------------------------------------------------------------------------------
def get_aptitude(xtalist_in):
   '''It obtains the aptitude of the structures received.

   in: xtalist_in (List); The structures of interest
   out: aptitude (List); Each entry will concatenate with the fitness of each structure
   '''
   if len(xtalist_in)==1:
      aptitude=[float(1.0)]
   else:
      aptitude=[]
      sort_xtal = sort_by_energy(xtalist_in,1)
      Emin = sort_xtal[0].e
      Emax = sort_xtal[-1].e
      for ixtal in sort_xtal:
         Ei = ixtal.e
         EFEi = (Ei - Emin) / (Emax - Emin)
         fi = 0.5 * (1.0 - np.tanh(((2.0 * EFEi) - 1.0)))
         aptitude.append(float(fi))
   return aptitude
#------------------------------------------------------------------------------------------