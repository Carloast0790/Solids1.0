import os
import os.path
import numpy as np
from utils_solids.libmoleculas import Atom, Molecule
from inout_solids.getbilparam import get_a_str
log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
def get_normaltermination_vasp(path):
    filename = 'OUTCAR'
    if os.path.isfile(path+filename):
        file = open(path+filename,'r')
        normal = 0
        for line in file:
            if "General timing and accounting informations for this job" in line: 
                normal = 1
                break
        file.close()
        return normal
    else:
        return False
#------------------------------------------------------------------------------------------
#eV = 1.602176620898*10^-19 Joule*(calorie/(4.184 Joule))*(kilo/1000)*((6.02214085774*10^23 permol)/NA)=23.0605 kilocalpermol/NA
##The relevant energy for molecules and atoms is energy without entropy.
##if "free  energy   TOTEN" in line:
##     eneline=line.split()
##     eneineV=float(eneline[4])
##The total energy will only be correct in the limit sigma->0
#------------------------------------------------------------------------------------------
def get_energy_vasp(path):
    outcarfile = open(path+'OUTCAR','r')
    eneineV = 0.0
    for line in outcarfile:
        if "energy  without entropy" in line:
             eneline = line.split()
             if len(eneline) > 5:
                eneineV = float(eneline[6])
    outcarfile.close()
    # eVtokcalpermol=float(23.0605)
    # eneinkcalpermol=eneineV*eVtokcalpermol
    return eneineV
#------------------------------------------------------------------------------------------
def get_geometry_vasp(path):
    filename='CONTCAR'
    print(path+filename)
    contcarfile=open(path+filename,'r')

    line=contcarfile.readline()
    name=line.strip()
    name=name.replace(" ", "")

    line=contcarfile.readline()
    if line=="":
        return False
    else:
        a0=float(line.split()[0])

    line=contcarfile.readline()
    a1x, a1y, a1z=map(float,line.split())
    line=contcarfile.readline()
    a2x, a2y, a2z=map(float,line.split())
    line=contcarfile.readline()
    a3x, a3y, a3z=map(float,line.split())
    matrix=np.array([[a0*a1x, a0*a1y, a0*a1z],[a0*a2x, a0*a2y, a0*a2z],[a0*a3x, a0*a3y, a0*a3z]])
    line=contcarfile.readline()
    elements=line.split()
    line=contcarfile.readline()
    ocupnumchar=line.split()
    ocupnuminte=list(map(int,ocupnumchar))
    natom=sum(ocupnuminte)
    liste, kk=[], 0
    for ii in ocupnuminte:
        for jj in range(ii):
            liste.append(elements[kk])
        kk=kk+1
    energy=get_energy_vasp(path)
    poscarx=Molecule(name, energy, matrix)
    ##VERIFICAR ESTA LINEA .....
    #namein=foldername.split('/')[-1]
    for line in contcarfile:
        sd=0
        if 'Selective dynamics' in line:
            line=contcarfile.readline()
            sd=1
        if 'Direct' in line:
            for iatom in range(natom):
                line=contcarfile.readline()
                vecxyz=line.split()
                s=liste[iatom]
                xd=float(vecxyz[0])
                yd=float(vecxyz[1])
                zd=float(vecxyz[2])
                xc=a0*(a1x*xd+a2x*yd+a3x*zd)
                yc=a0*(a1y*xd+a2y*yd+a3y*zd)
                zc=a0*(a1z*xd+a2z*yd+a3z*zd)
                if sd==1:
                    xf=str(vecxyz[3])
                    yf=str(vecxyz[4])
                    zf=str(vecxyz[5])
                    xatom=Atom(s,xc,yc,zc,xf,yf,zf)
                elif sd==0:
                    xatom=Atom(s,xc,yc,zc)
                poscarx.add_atom(xatom)
        if 'Cartesian' in line:
            for iatom in range(natom):
                line=contcarfile.readline()
                vecxyz=line.split()
                s=liste[iatom]
                xc=float(vecxyz[0])
                yc=float(vecxyz[1])
                zc=float(vecxyz[2])
                if sd==1:
                    xf=str(vecxyz[3])
                    yf=str(vecxyz[4])
                    zf=str(vecxyz[5])
                    xatom=Atom(s,xc,yc,zc,xf,yf,zf)
                elif sd==0:
                    xatom=Atom(s,xc,yc,zc)
                poscarx.add_atom(xatom)
    contcarfile.close()
    return poscarx
#------------------------------------------------------------------------------------------
def get_nt_geometry_vasp(path):
    filename = 'OUTCAR'
    r = get_normaltermination_vasp(path)
    if r is False:
        fopen = open(log_file,'a')
        print("%s does NOT EXIST!!" %(path+filename), file=fopen)
        fopen.close()
        return False
    elif r == 0:
        fopen = open(log_file,'a')
        print("ABNORMAL Termination in %s" %(path+filename), file=fopen)
        fopen.close()
        return False
    else:
        fopen = open(log_file,'a')
        print("Normal Termination in %s" %(path+filename), file=fopen)
        fopen.close()
        moleculeout=get_geometry_vasp(path)
        moleculeout.c=[r]
        return moleculeout
#------------------------------------------------------------------------------------------
def get_all_nt_geometry_vasp(generationfolder, moleculelist):
    moleculeout=[]
    for imol in moleculelist:
        foldername=imol.i
        path=generationfolder+foldername+'/'
        mol01=get_nt_geometry_vasp(path)
        if mol01 is not False:
            mol01.i=foldername
            moleculeout.extend([mol01])
    return moleculeout
#------------------------------------------------------------------------------------------
def get_at_geometry_vasp(path):
    r = get_normaltermination_vasp(path)
    if r is False: return False
    elif r == 0:
        if get_energy_vasp(path) == float(0.0): return False
        else:
            moleculeout = get_geometry_vasp(path)
            moleculeout.c = [r]
            return moleculeout
    else: return False
#------------------------------------------------------------------------------------------
def get_all_at_geometry_vasp(generationfolder, moleculelist, erase=0):
    moleculeout = []
    for imol in moleculelist:
        foldername=imol.i
        path=generationfolder+foldername+'/'
        mol01=get_at_geometry_vasp(path)
        if mol01 is not False:
            mol01.i=foldername
            moleculeout.extend([mol01])
            if (erase==1):
                os.remove(path+'OUTCAR')
                os.remove(path+foldername+'.out')
                fopen = open(log_file,'a')
                print("Move: %sCONTCAR -> %sPOSCAR" %(path, path), file=fopen)
                fopen.close()
                os.system('mv '+path+'CONTCAR '+path+'POSCAR')
    return moleculeout
#------------------------------------------------------------------------------------------
def get_xt_geometry_vasp(path):
    r = get_normaltermination_vasp(path)
    if r is False: return False
    elif r == 0:
        if get_energy_vasp(path)==float(0.0): return False
        else:
            moleculeout=get_geometry_vasp(path)
            moleculeout.c=[r]
            return moleculeout
    else:
        moleculeout=get_geometry_vasp(path)
        moleculeout.c=[r]
        return moleculeout
#------------------------------------------------------------------------------------------
def get_all_xt_geometry_vasp(generationfolder, moleculelist):
    moleculeout, count, fall=[], 0, 0
    for imol in moleculelist:
        foldername=imol.i
        path=generationfolder+foldername+'/'
        mol01=get_xt_geometry_vasp(path)
        if mol01 is not False:
            mol01.i=foldername
            moleculeout.extend([mol01])
            count=count+1
        else:
            fall=fall+1
    fopen = open(log_file,'a')
    print("Number of NON-recoverable calculations = %d" %(fall), file=fopen)
    print("Recoverable + SUCCESSFUL calculations  = %d" %(count), file=fopen)
    fopen.close()
    return moleculeout
#------------------------------------------------------------------------------------------
def get_mmt_from_outcar(path):
    filename='OUTCAR'
    mmt=float(0.0)
    if os.path.isfile(path+filename):
        outcarfile=open(path+filename,'r')
        for line in outcarfile:
            if "magnetization (x)" in line:
                count=0
                while (count < 1):
                    line=outcarfile.readline()
                    linetest=line.split()
                    if ('tot' in linetest) and ('#' not in linetest):
                        mmt=float(linetest[4])
                        count=1
    return mmt
#------------------------------------------------------------------------------------------
