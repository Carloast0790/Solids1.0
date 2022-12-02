import numpy as np
from utils.libmoleculas import sort_by_stoichiometry
from inout.getbilparam  import get_a_str
log_file=get_a_str('output_file','glomos_out.txt')
#------------------------------------------------------------------------------------------
from ase import Atoms
from ase.ga.ofp_comparator import OFPComparator
#------------------------------------------------------------------------------------------
def xyz2ase(moleculein):
    symbols, positions = [], []
    for iatom in moleculein.atoms:
        symbols.append(iatom.s)
        positions.append([iatom.xc, iatom.yc, iatom.zc])
    if np.any(moleculein.m):
        celda=moleculein.m
        moleculeout=Atoms(symbols=symbols, positions=positions, cell=celda, pbc=[2, 2, 2])
    else:
        moleculeout=Atoms(symbols=symbols, positions=positions)
    return moleculeout
#------------------------------------------------------------------------------------------
def compare_oganov_fp(moleculeinx, moleculeiny):
    if np.any(moleculeinx.m) and np.any(moleculeiny.m):
        mypbc=[True, True, True]
    else:
        mypbc=[False, False, False]
    mol=sort_by_stoichiometry([moleculeinx,moleculeiny])
    comp=OFPComparator(
    dE=1.0,
    cos_dist_max=5e-3,
    rcut=20.0,
    binwidth=0.05,
    pbc=mypbc,
    sigma=0.02,
    nsigma=4,
    recalculate=True)
    mola=xyz2ase(mol[0])
    molb=xyz2ase(mol[1])
    ans=comp.looks_like(mola, molb)
    return ans
#------------------------------------------------------------------------------------------
def discriminate_all_OFP(xtal_list):
    xtal_out = xtal_list.copy()
    org_len = len(xtal_list)
    for i,xtal1 in enumerate(xtal_list):
        flag = True
        print('----str ',xtal1.i,'----')
        for j in range(i+1,org_len):
            xtal2 = xtal_list[j]
            ans = compare_oganov_fp(xtal1,xtal2)
            print('comparison with',xtal2.i,'result',ans)
            if ans:
                print('removing',xtal1.i)
                xtal_out.remove(xtal1)
                break
    return xtal_out

#------------------------------------------------------------------------------------------
def run_sample():
    from vasp.libperiodicos import readposcars, writeposcars, expand_poscar
    x = readposcars('reported_rutilo.vasp')[0]
    y = readposcars('found_rutilo.vasp')[0]
    ans = compare_oganov_fp(x,x)
    print(ans)
run_sample()
