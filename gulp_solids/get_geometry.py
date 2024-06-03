import os.path
import numpy as np
from utils_solids.atomic import get_chemical_symbol
from utils_solids.libmoleculas import Atom, Molecule, writexyzs
from vasp_solids.libperiodicos import direct2cartesian, writeposcars
from inout_solids.getbilparam import get_a_str
log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
def get_normaltermination_gulp(path,filename):
    if os.path.isfile(path+filename):
        file=open(path+filename,'r')
        for line in file:
            if '**** Optimisation achieved ****' in line:
                normal = 1
                break
            else:
                normal = 0
        file.close()
        return normal
    else:
        return False
#------------------------------------------------------------------------------------------
def get_energy_gulp(path,filename):
    if not os.path.isfile(path+filename):
        fopen = open(log_file,'a')
        print("The file %s does not exist." %(path+filename),file=fopen)
        fopen.close()
        exit()
    file = open(path+filename,'r')
    for line in file:
        if 'Final energy =' in line:
            ls = line.split()
            eneineV = float(ls[3])
        if 'Final enthalpy =' in line:
            ls = line.split()
            eneineV = float(ls[3])
    file.close()
    #eVtokcalpermol=float(23.0609)
    #eneinkcalpermol=eneineV*eVtokcalpermol
    #return eneinkcalpermol
    return eneineV
#------------------------------------------------------------------------------------------
def get_geometry_gulp(path, filename):
    if not os.path.isfile(path+filename):
        fopen = open(log_file,'a')
        print("The file %s does not exist." %(path+filename),file=fopen)
        fopen.close()
        exit()
    file = open(path+filename,'r')
    ans1 = 'Final Cartesian lattice vectors (Angstroms) :'
    ans2 = 'Cartesian lattice vectors (Angstroms) :'
    for line in file:
        if (ans1 in line) or (ans2 in line):
            line=file.readline()
            line=file.readline()
            ls = line.split()
            a1x, a1y, a1z=float(ls[0]), float(ls[1]), float(ls[2])
            line=file.readline()
            ls = line.split()
            a2x, a2y, a2z=float(ls[0]), float(ls[1]), float(ls[2])
            line=file.readline()
            ls = line.split()
            a3x, a3y, a3z=float(ls[0]), float(ls[1]), float(ls[2])
    matrix=np.array([[a1x, a1y, a1z],[a2x, a2y, a2z],[a3x, a3y, a3z]])
    file.close()
    name=filename.split('.')[0]
    energy=get_energy_gulp(path,filename)
    poscarx=Molecule(name,energy,matrix)
    file=open(path+filename,'r')
    for line in file:
        if 'Final fractional coordinates of atoms :' in line:
            for ii in range(5): line=file.readline()
            line=file.readline()
            ls = line.split()
            lenls=len(ls)
            while ( lenls == 7 ):
                s=str(ls[1])
                xd,yd,zd=float(ls[3]),float(ls[4]),float(ls[5])
                xc,yc,zc=direct2cartesian(xd, yd, zd, matrix)
                ai=Atom(s,xc,yc,zc)
                poscarx.add_atom(ai)
                line=file.readline()
                ls = line.split()
                lenls=len(ls)
    file.close()
    return poscarx
#------------------------------------------------------------------------------------------
def get_xt_geometry_gulp(path,filename):
    r = get_normaltermination_gulp(path,filename)
    if r is False or r == 0:
        return False
    else:
        poscarout = get_geometry_gulp(path,filename)
        poscarout.c = [r]
        return poscarout
#------------------------------------------------------------------------------------------
def get_all_xt_geometry_gulp(poscarlist, path):
    poscarout,count=[],0
    for ipos in poscarlist:
        iposname = ipos.i
        filename = iposname + '.got'
        pos01 = get_xt_geometry_gulp(path,filename)
        if pos01:
            poscarout.extend([pos01])
            count = count + 1
    if count==0: return False
    else: return poscarout
#------------------------------------------------------------------------------------------
def run_sample():
    from vasp.libperiodicos import writeposcars
    mol=get_xt_geometry_gulp('ejemplos/','anatasa.got')
    print(mol.c[0])
    writeposcars([mol],'out.vasp','C')
#run_sample()

