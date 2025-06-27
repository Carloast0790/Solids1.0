import numpy as np
from dscribe.descriptors import CoulombMatrix,  ACSF, SOAP, MBTR, LMBTR
#---descriptors----------------------------------------------------------------

def descriptors_MBTR (lista):
    des_list = []
    for imol in lista:
        atoms = imol.get_chemical_symbols()
        especies = set(atoms)
        geometria={"function": "inverse_distance"}
        grid={"min": 0, "max": 1, "sigma": 0.1, "n" : 100}
        pesaje={"function": "exp", "scale": 0.5, "threshold" :1e-3}
        mbtr = MBTR( species = especies, geometry = geometria, weighting = pesaje,  grid=grid, periodic=False, normalization="l2")
        des = mbtr.create(imol)
        des_list.append(des)
    return des_list

def descriptors_MBTR_C (lista):
    des_list = []
    for imol in lista:
         atoms = imol.get_chemical_symbols()
         especies = set(atoms)
         geometria_d ={"function": "inverse_distance"}
         geometria_c ={"function": "cosine"}
         grid={"min": 0, "max": 1, "sigma": 0.1, "n" : 100}
         #pesaje_d={"function": "inverse_square", "r_cut": 10, "threshold": 1e-3}
         pesaje={"function": "exp", "scale": 0.5, "threshold": 1e-3}
         mbtr_d = MBTR( species = especies, geometry = geometria_d , weighting = pesaje,  grid=grid, periodic=False, normalization="l2")
         mbtr_c = MBTR( species = especies, geometry = geometria_c , weighting = pesaje,  grid=grid, periodic=False, normalization="l2")
         des_d = mbtr_d.create(imol)
         des_c = mbtr_c.create(imol)
         des = np.concatenate((des_d, des_c))
         des_list.append(des)
    return des_list

def descriptors_LMBTR (lista):
    des_list = []
    for imol in lista:
         geometria={"function": "distance"}
         grid={"min": 0, "max": 5, "n": 100, "sigma": 0.1}
         pesaje={"function": "exp", "scale": 0.5, "threshold": 1e-3}
         atoms = imol.get_chemical_symbols()
         especies = set(atoms)
         lmbtr = LMBTR(species = especies, geometry=geometria, grid=grid, weighting=pesaje, periodic=False, normalization="l2")
         lmbtr_des = lmbtr.create(imol)
         des = lmbtr_des.flatten()
         des_list.append(des)
    return des_list

def descriptors_SOAP (lista):
    des_list = []
    for imol in lista:
        atoms = imol.get_chemical_symbols()
        especies = set(atoms)
        r_cut = 6.0
        n_max = 8
        l_max = 6
        soap = SOAP(species = especies, periodic=False, r_cut = r_cut, n_max=n_max, l_max=l_max)
        soap_des = soap.create(imol)
        des = soap_des.flatten()
        des_list.append(des)
    return des_list

def descriptors_CoulombMatrix (lista):
    des_list = []
    for imol in lista:
        filas = len(imol)
        cm = CoulombMatrix(n_atoms_max = filas)
        des = cm.create(imol)
        des_list.append(des)
    return des_list

def descriptors_ACSF (lista):
    des_list = []
    for imol in lista:
         g2params=[[1, 1], [1, 2], [1, 3]]
         g4params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, -1]]
         atoms = imol.get_chemical_symbols()
         especies = set(atoms)
         r_cut = 6.0
         acsf = ACSF(species = especies, r_cut = r_cut, g2_params=g2params, g4_params=g4params)
         acsf_des = acsf.create(imol)
         des = acsf_des.flatten()
         des_list.append(des)
    return des_list

#def descriptors_USR (imol):
    #des_list = []
    #for imol in lista:
        #des = USR(imol)
        #des_list.append(des)
    #return des_list