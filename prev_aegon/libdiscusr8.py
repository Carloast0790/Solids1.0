import os
import time
import numpy as np
from multiprocessing import Process
from aegon.libutils import readxyzs, writexyzs, sort_by_energy, four_points, prepare_folders, split_poscarlist
tolene=0.25
#------------------------------------------------------------------------------------------
def lastwo_central_moment(listxxx):
    xavg=np.mean(listxxx)
    res = listxxx - xavg
    ssq2 = np.mean(res**2) 
    ssq3 = np.mean(res**3) 
    return [ssq2,ssq3]
#------------------------------------------------------------------------------------------
#Ultrafast Shape Recognition (USR) Algorithm
def usrcreate(moleculein):
    ctd, cst, fct, ftf=four_points(moleculein)
    lctd,  lcst,  lfct,  lftf = [], [], [], []
    for iatom in moleculein:
        xyz=np.array(iatom.position)
        ictd1=np.linalg.norm(xyz-ctd)
        icst1=np.linalg.norm(xyz-cst)
        ifct1=np.linalg.norm(xyz-fct)
        iftf1=np.linalg.norm(xyz-ftf)
        lctd.append(ictd1)
        lcst.append(icst1)
        lfct.append(ifct1)
        lftf.append(iftf1)
    a1=lastwo_central_moment(lctd)
    a2=lastwo_central_moment(lcst)
    a3=lastwo_central_moment(lfct)
    a4=lastwo_central_moment(lftf)
    standard_USR=a1+a2+a3+a4
    return standard_USR
#------------------------------------------------------------------------------------------
def find_similar_elements(similarity_matrix, threshold):
    similar_elements_indices = []
    num_elements = similarity_matrix.shape[0]
    for i in range(num_elements):
        for j in range(i+1,num_elements):
            if similarity_matrix[i, j] >= threshold:
                similar_elements_indices.append(j)
    return similar_elements_indices
#------------------------------------------------------------------------------------------
def disc_USR(atoms, threshold):
    num_molecules = len(atoms)
    descriptors=[usrcreate(imol) for imol in atoms]
    n_features=8
    similar_elements_indices = []
    similarity_matrix = np.zeros((num_molecules, num_molecules))
    for i in range(num_molecules):
        vi=descriptors[i]
        for j in range(i, num_molecules):
            vj=descriptors[j]
            manhattan_distance=sum([np.absolute(a-b) for a,b in zip(vi,vj)])
            similarity=1.0/(1.0+manhattan_distance/float(n_features))
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
    similar_elements_indices=find_similar_elements(similarity_matrix, threshold)
    similar_elements_indices.sort()
    disimilars_atoms=[atoms[i] for i in range(num_molecules) if i not in similar_elements_indices]
    return disimilars_atoms
#------------------------------------------------------------------------------------------
def comparator_usr_conv(atoms0, threshold):
    start = time.time()
    ni=len(atoms0)
    atoms1=disc_USR(atoms0, threshold)
    nf=len(atoms1)
    end = time.time()
    print('USR comparison (-serial-) at %5.2f s [%d -> %d]' %(end - start, ni, nf))
    return atoms1
#------------------------------------------------------------------------------------------
def make_comparator_usr(ifolder, threshold):
    atoms0=readxyzs(ifolder+'/'+ifolder+'.xyz')
    atoms1=disc_USR(atoms0, threshold)
    writexyzs(atoms1,ifolder+'/'+ifolder+'_disc.xyz',1)
#------------------------------------------------------------------------------------------
def comparator_usr_parallel(poscarlist, threshold, nproc, base_name):
    poscarlist=sort_by_energy(poscarlist, 1)
    ni=len(poscarlist)
    if ni < 4*nproc:
        moleculeout=comparator_usr_conv(poscarlist, threshold)
        return moleculeout
    start = time.time()
    folderlist=prepare_folders(poscarlist, nproc, base_name)
    poscar_split_list=split_poscarlist(poscarlist, nproc)
    procs = []
   #dicc_term = {iposcar.info['i']: iposcar.info['c'] for iposcar in poscarlist}
    for ifolder, iposcars in zip(folderlist, poscar_split_list):
        writexyzs(iposcars, ifolder+'/'+ifolder+'.xyz',1)
        proc = Process(target=make_comparator_usr, args=(ifolder,threshold,))
        procs.append(proc)
        proc.start()
    for proc in procs:
        proc.join()
    center_mols, border_mols = [], []
    for ifolder in folderlist:
        molx=readxyzs(ifolder+'/'+ifolder+'_disc.xyz')
       #for imol in molx: imol.info['c']=dicc_term[imol.info['i']]
        molx=sort_by_energy(molx, 1)
        emin=molx[0].info['e']
        emax=molx[-1].info['e']
        for imol in molx:
            delw=imol.info['e'] - emin
            deup=emax - imol.info['e']
            if ( delw < tolene ) or ( deup < tolene ):
                border_mols.extend([imol])
            else:
                center_mols.extend([imol])
    border_mols=disc_USR(border_mols, threshold)
    moleculeout=center_mols+border_mols
    moleculeout=sort_by_energy(moleculeout, 1)
    os.system('rm -rf %sproc[0-9][0-9]' %(base_name))
    nf=len(moleculeout)
    end = time.time()
    print('USR comparison (parallel) at %5.2f s [%d -> %d]' %(end - start, ni, nf))
    return moleculeout
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#from dask.distributed import Client
#def usr_descriptors_two(atoms):
#    information=[usrcreate(imol) for imol in atoms]
#    return information
#------------------------------------------------------------------------------------------
#def disc_USR_two(atoms, threshold, nproc):
#    n_features=8
#    num_molecules = len(atoms)
#    client = Client(n_workers=nproc)
#    lista=split_poscarlist(atoms, nproc)
#    futures = client.map(usr_descriptors_two, lista)
#    listofthinks = client.gather(futures)
#    client.close()
#    descriptors=[]
#    for ithink in listofthinks: descriptors=descriptors+ithink
#    similar_elements_indices = []
#    similarity_matrix = np.zeros((num_molecules, num_molecules))
#    for i in range(num_molecules):
#        vi=descriptors[i]
#        for j in range(i, num_molecules):
#            vj=descriptors[j]
#            manhattan_distance=sum([np.absolute(a-b) for a,b in zip(vi,vj)])
#            similarity=1.0/(1.0+manhattan_distance/float(n_features))
#            similarity_matrix[i, j] = similarity
#            similarity_matrix[j, i] = similarity
#    similar_elements_indices=find_similar_elements(similarity_matrix, threshold)
#    similar_elements_indices.sort()
#    disimilars_atoms=[atoms[i] for i in range(num_molecules) if i not in similar_elements_indices]
#    return disimilars_atoms
#------------------------------------------------------------------------------------------
#def comparator_usr_parallel_two(atoms0, threshold, nproc):
#    start = time.time()
#    ni=len(atoms0)
#    atoms1=disc_USR_two(atoms0, threshold, nproc)
#    nf=len(atoms1)
#    end = time.time()
#    print('USR comparison (parallel descriptors) at %5.2f s [%d -> %d]' %(end - start, ni, nf))
#    return atoms1
#------------------------------------------------------------------------------------------