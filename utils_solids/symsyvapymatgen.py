import os
from utils_solids.atomic import get_atomic_number, get_chemical_symbol
from utils_solids.libmoleculas import Atom, Molecule, align_all_inertia_axis_x
from inout_solids.getbilparam import get_a_str
from pymatgen.core import Molecule as MoleculePMG
from pymatgen.symmetry.analyzer import PointGroupAnalyzer, iterative_symmetrize
log_file=get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
#https://pymatgen.org/pymatgen.symmetry.html
#------------------------------------------------------------------------------------------
def moleculeglomos2pymatgen(moleculeglomos):
    coords, atoms=[],[]
    for iatom in moleculeglomos.atoms:
        atoms.append(iatom.s)
        coords.append([iatom.xc, iatom.yc, iatom.zc])
    moleculepymatgen=MoleculePMG(atoms,coords)
    return moleculepymatgen
#------------------------------------------------------------------------------------------
def moleculepymatgen2glomos(moleculepymatgen, name, energy):
    atomic_numbers=moleculepymatgen.atomic_numbers
    cart_coords=moleculepymatgen.cart_coords
    moleculeglomos=Molecule(name,energy)
    for ii,iatom in enumerate(atomic_numbers):
        s=get_chemical_symbol(iatom)
        xc=cart_coords[ii,0]
        yc=cart_coords[ii,1]
        zc=cart_coords[ii,2]
        ai=Atom(s,xc,yc,zc)
        moleculeglomos.add_atom(ai)
    return moleculeglomos
#------------------------------------------------------------------------------------------
def point_group(moleculeglomos, tolerance=0.25, eigen_tolerance=1E-3, matrix_tolerance=0.1):
    molpmg=moleculeglomos2pymatgen(moleculeglomos)
    molx=PointGroupAnalyzer(molpmg, tolerance, eigen_tolerance, matrix_tolerance)
    pointgroup=molx.get_pointgroup()
    return pointgroup
#------------------------------------------------------------------------------------------
def symmetrize_molecule(moleculeglomos, tolerance=1e-05):
    molpmg=moleculeglomos2pymatgen(moleculeglomos)
    dictionary=iterative_symmetrize(molpmg, max_n=10000, tolerance=tolerance, epsilon=1e-06)
    molf=dictionary['sym_mol']
    name, energy=moleculeglomos.i, moleculeglomos.e
    molgms=moleculepymatgen2glomos(molf, name, energy)
    return molgms
#------------------------------------------------------------------------------------------
def round_molecule(moleculeglomos, snumbers):
    #tol=float(pow(10,-snumbers))
    for iatom in moleculeglomos.atoms:
        iatom.xc=round(iatom.xc,snumbers)
        iatom.yc=round(iatom.yc,snumbers)
        iatom.zc=round(iatom.zc,snumbers)
    return moleculeglomos
#------------------------------------------------------------------------------------------
def syva(moleculeglomos):
    natoms=moleculeglomos.n
    namein=moleculeglomos.i
    energy=moleculeglomos.e
    fh=open('syvamolinput',"w")
    print(moleculeglomos.i, file=fh)
    print(moleculeglomos.n, file=fh)
    for iatom in moleculeglomos.atoms:
        ani=get_atomic_number(iatom.s)
        print("%-2d %16.9f %16.9f %16.9f" %(ani, iatom.xc, iatom.yc, iatom.zc), file=fh)
    fh.close()
    os.system("syva all syvamolinput > syvamolout")
    syvafile=open('syvamolout','r')
    moleculeout=[]
    for line in syvafile:
        if "Optimized" in line:
            hls=line.split()
            pg=str(hls[1])
            moltmp=Molecule(namein,energy)
            line=syvafile.readline()
            for ii in range(natoms):
                line=syvafile.readline()
                ls=line.split()
                ss=get_chemical_symbol(int(ls[0]))
                xc,yc,zc = float(ls[1]), float(ls[2]), float(ls[3])
                ai=Atom(ss,xc,yc,zc)
                moltmp.add_atom(ai)
            moltmp.c=pg
            moleculeout.extend([moltmp])
    syvafile.close()
    os.system("rm -f syvamolinput syvamolout")
    return moleculeout[-1]
#------------------------------------------------------------------------------------------
def sym_syva(moleculelist):
    moleculeout2=[]
    for imol in moleculelist:
        mol1=syva(imol)
        comentario=mol1.c
        for x in [1e-03, 1e-06]:
            mol1=align_all_inertia_axis_x([mol1])[0]
            mol1=round_molecule(mol1,7)
            mol1=symmetrize_molecule(mol1, x)
        pg=point_group(mol1, x, 0.01, 0.1)
        #mol1.i=mol1.i+'_'+str(pg)
        fopen = open(log_file,'a')
        print('%-11s run syva:%-3s -> read pymatgen:%-3s' %(mol1.i, comentario, pg), file=fopen)
        fopen.close()
        moleculeout2.extend([mol1])
    return moleculeout2
#------------------------------------------------------------------------------------------
#from utils.libmoleculas  import readxyzs, writexyzs
#mol=readxyzs('summary.xyz')
#mol=sym_syva(mol)
#writexyzs(mol,'hola.xyz')
