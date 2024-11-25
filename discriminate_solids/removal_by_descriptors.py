import numpy as np
from ase import Atoms
from dscribe.descriptors import ValleOganov
from inout_solids.getbilparam import get_a_str, get_a_float

log_file = get_a_str('output_file','solids_out.txt')
# simil_tol = get_a_float('similarity_tolerance',0.8)
#------------------------------------------------------------------------------------------
def solids2ase(solidslist):
    aselist = []
    for sxtal in solidslist:
        symbols=[]
        positions=[]
        for satm in sxtal.atoms:
            dx,dy,dz = satm.xc, satm.yc, satm.zc
            dx,dy,dz = round(dx,6),round(dy,6),round(dz,6)
            symbols.append(satm.s)
            positions.append([dx,dy,dz])
        axtal = Atoms(symbols=symbols, positions=positions, cell=sxtal.m, pbc=True)
        aselist.append(axtal)
    return aselist

#------------------------------------------------------------------------------------------
def descriptor_comparison_calculated(xtalist_in, tolerance):
    fopen = open(log_file,'a')
    print('\n-------------------------------------------------------------------',file=fopen)
    print('--------------- Duplicates Removal in Generation ------------------',file=fopen)
    print('\nTolerance Given: ' + str(tolerance),file=fopen)
    size_in = len(xtalist_in)
    xtalist_out = []
    aselist_in = solids2ase(xtalist_in)
    species = list(aselist_in[0].get_chemical_symbols())
    species = set(species)
    vo = ValleOganov(species=species, function='distance', n=100, sigma=1E-5, r_cut=10)
    descriptors = [vo.create(structure) for structure in aselist_in]
    disc_count = 0
    for i in range(len(descriptors)):
        stop_flag = False
        for j in range(i+1, len(descriptors)):
            norm_i = np.linalg.norm(descriptors[i])
            norm_j = np.linalg.norm(descriptors[j])
            dot_product = np.dot(descriptors[i], descriptors[j])
            similarity = dot_product / (norm_i * norm_j)
            if similarity >= tolerance:
                print('%s removed, too similar to %s, similarity = %.5f' %(xtalist_in[i].i,xtalist_in[j].i,similarity),file=fopen)
                disc_count = disc_count + 1 
                stop_flag = True
                break
        if stop_flag == False:
            xtalist_out.append(xtalist_in[i])
    print('\n'+str(disc_count)+' structures removed by similarity in generation comparison',file=fopen)
    fopen.close()
    return xtalist_out

#------------------------------------------------------------------------------------------
def descriptor_comparison_calculated_vs_pool(xtal_calculated, xtal_pool, tolerance):
    fopen = open(log_file,'a')
    print('\n-------------------------------------------------------------------',file=fopen)
    print('---------------- Duplicates Removal Gen vs Pool -------------------\n',file=fopen)
    different_calc = []
    xtalist_out = []
    calculated_ase = solids2ase(xtal_calculated)
    pool_ase = solids2ase(xtal_pool)
    species = list(calculated_ase[0].get_chemical_symbols())
    species = set(species)
    vo = ValleOganov(species=species, function='distance', n=100, sigma=1E-5, r_cut=10)
    descr_calc = [vo.create(structure) for structure in calculated_ase]
    descr_pool = [vo.create(structure) for structure in pool_ase]
    disc_count = 0
    for i in range(len(descr_calc)):
        stop_flag = False
        for j in range(len(descr_pool)):
            norm_i = np.linalg.norm(descr_calc[i])
            norm_j = np.linalg.norm(descr_pool[j])
            dot_product = np.dot(descr_calc[i], descr_pool[j])
            similarity = dot_product / (norm_i * norm_j)
            if similarity >= tolerance:
                print('%s removed, too similar to %s, similarity = %.5f' %(xtal_calculated[i].i,xtal_pool[j].i,similarity),file=fopen) 
                stop_flag = True
                disc_count = disc_count + 1
                break
        if stop_flag == False:
            different_calc.append(xtal_calculated[i])
    if different_calc:
        print('\n'+str(disc_count)+' structures removed by similarity in Gen vs Pool comparison'+'\n',file=fopen)
        xtalist_out.extend(different_calc)
    else:
        print('\nZero structures removed by similarity in Gen vs Pool comparison'+'\n',file=fopen)
    fopen.close()
    return xtalist_out
