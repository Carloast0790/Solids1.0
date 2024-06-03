import os.path
import numpy as np
import random
from utils_solids.atomic  import get_chemical_symbol, get_atomic_mass, get_covalent_radius
from utils_solids.general import euler_matrix, rodrigues_rotation_matrix
#------------------------------------------------------------------------------------------
class Atom:
    def __init__(self, symbol, xcartesian, ycartesian, zcartesian, xfix='T', yfix='T', zfix='T'):
        self.s=symbol
        self.xc=xcartesian
        self.yc=ycartesian
        self.zc=zcartesian
        self.xf=xfix
        self.yf=yfix
        self.zf=zfix
    def vec(self):
        vector=np.array([self.xc,self.yc,self.zc])
        return vector
#------------------------------------------------------------------------------------------
class Molecule:
    def __init__(self, name, energy, matrix=[], comments=[]):
        self.atoms = []
        self.i = name
        self.e = energy
        self.m = matrix
        self.c = comments
        self.n = 0
    def add_atom(self, atom):
        self.atoms.append(atom)
        natoms=len(self.atoms)
        self.n = natoms
    def clear(self):
        delattr(self, 'i')
        delattr(self, 'e')
        delattr(self, 'm')
        delattr(self, 'c')
        delattr(self, 'n')
        delattr(self, 'atoms')
#------------------------------------------------------------------------------------------
def readxyzs(filename):
    """ SAMPLE OF USE IN PYTHON:
    >>> from utils.libmoleculas import readxyzs
    >>> mol0=readxyzs('CoenzymeA.xyz')
    >>> print mol0[0].i, mol0[0].e, mol0[0].n
    """
    if not os.path.isfile(filename):
        print("The file",filename,"does not exist.")
        exit()
    file=open(filename,'r')
    imol=-1
    moleculeout=[]
    for line in file:
        ls=line.split()
        if len(ls)==1:
            natoms=int(ls[0])
            count=0
            imol=imol+1
            line=file.readline()
            ls=line.split()
            if len(ls)==0: name,energy='unknown', float(0.0)
            if len(ls)==1: name,energy=str(ls[0]),float(0.0)
            if len(ls)>=2: name,energy=str(ls[1]),float(ls[0])
            mol=Molecule(name,energy)
        if len(ls)==4:
            s=get_chemical_symbol(ls[0])
            xc,yc,zc=float(ls[1]),float(ls[2]),float(ls[3])
            ai=Atom(s,xc,yc,zc)
            mol.add_atom(ai)
            count=count+1
            if count==natoms: moleculeout.extend([mol])
    file.close()
    return moleculeout
#------------------------------------------------------------------------------------------
def readcc1s(filename):
    file=open(filename,'r')
    imol=-1
    moleculeout=[]
    for line in file:
        ls=line.split()
        if len(ls)==1:
            natoms=int(ls[0])
            count=0
            imol=imol+1
            name='unknown'+str(imol).zfill(3)
            energy=float(0.0)
            mol=Molecule(name,energy)
        if len(ls)>4:
            s=get_chemical_symbol(ls[0])
            xc,yc,zc=float(ls[2]),float(ls[3]),float(ls[4])
            ai=Atom(s,xc,yc,zc)
            mol.add_atom(ai)
            count=count+1
            if count==natoms: moleculeout.extend([mol])
    file.close()
    return moleculeout
#------------------------------------------------------------------------------------------
def molecular_radius_max(moleculein):
    moleculeout=copymol(moleculein)
    translate_to_cm(moleculeout)
    r=[]
    for iatom in moleculeout.atoms:
        xyz=np.array([iatom.xc, iatom.yc, iatom.zc])
        ri=get_covalent_radius(iatom.s)
        dist=np.linalg.norm(xyz)+ri
        r.append(dist)
    r.sort()
    rmax=r[-1]
    moleculeout.clear()
    return rmax
#------------------------------------------------------------------------------------------
def readmatrix(matrixfiletxt):
    mfile=open(matrixfiletxt,'r')
    ii=0
    matrix=np.zeros(shape=(3,3),dtype=np.float64)
    for line in mfile:
        ls=line.split()
        matrix[ii,0]=float(ls[0])
        matrix[ii,1]=float(ls[1])
        matrix[ii,2]=float(ls[2])
        ii=ii+1
    mfile.close()
    rmatrix=np.array(matrix)
    return(rmatrix)
#------------------------------------------------------------------------------------------
def writexyzs(moleculein,filename, in_log=0):
    """ SAMPLE OF USE IN PYTHON:
    >>> from utils.libmoleculas import readxyzs, writexyzs
    >>> mol0=readxyzs('CoenzymeA.xyz')
    >>> writexyzs(mol0,'w.xyz')
    >>> mol1=[mol0[0]]
    >>> writexyzs(mol1,'w.xyz')
    """
    fh=open(filename,"w")
    for xmol in moleculein:
        print(xmol.n, file=fh)
        print("%12.8f     %s" %(xmol.e, xmol.i), file=fh)
        for iatom in xmol.atoms:
            print("%-2s %16.9f %16.9f %16.9f" %(iatom.s, iatom.xc, iatom.yc, iatom.zc), file=fh)
    fh.close()
    if in_log==0:
        print("Writing %s" %(filename))
    if in_log==1:
        from inout.getbilparam import get_a_str
        log_file=get_a_str('output_file','solids_out.txt')
        fopen = open(log_file,'a')
        print("Writing %s" %(filename), file=fopen)
        fopen.close()
#------------------------------------------------------------------------------------------
def copymol(moleculein):
    moleculeout=Molecule(moleculein.i, moleculein.e, moleculein.m, moleculein.c)
    for iatom in moleculein.atoms:
        ss = iatom.s
        xx, yy, zz = iatom.xc, iatom.yc, iatom.zc
        xf, yf, zf = iatom.xf, iatom.yf, iatom.zf
        ai=Atom(ss, xx, yy, zz, xf, yf, zf)
        moleculeout.add_atom(ai)
    return moleculeout
#------------------------------------------------------------------------------------------
def rename_molecule(moleculein, basename, ndigist):
    for imol, xmol in enumerate(moleculein):
        xmol.i=basename+str(imol+1).zfill(ndigist)
    return moleculein
#------------------------------------------------------------------------------------------
def uniontab(moleculein):
    moleculeout=Molecule('uniontab', 0.0)
    for xmol in moleculein:
        for iatom in xmol.atoms:
            moleculeout.add_atom(iatom)
    return moleculeout
#------------------------------------------------------------------------------------------
def sort_by_energy(moleculelist, opt=0):
    moleculeout=[]
    if len(moleculelist) == 0: return moleculeout
    s=[[imol,xmol.e] for imol,xmol in enumerate(moleculelist)]
    t = sorted(s, key=lambda x: float(x[1]))
    energy_ref = t[0][1] if (opt==0) else float(0.0)
    for ii in t:
        moltmp=copymol(moleculelist[ii[0]])
        moltmp.e=ii[1] - energy_ref
        moleculeout.extend([moltmp])
    return moleculeout
#------------------------------------------------------------------------------------------
def centroid(moleculein):
    natoms=moleculein.n
    xc, yc, zc = 0.0, 0.0, 0.0
    for iatom in moleculein.atoms:
        xc=xc+(iatom.xc)
        yc=yc+(iatom.yc)
        zc=zc+(iatom.zc)
    vectorc=np.array([xc, yc, zc])/float(natoms)
    return vectorc
#------------------------------------------------------------------------------------------
def center_off_mass(moleculein, black_list=[]):
    xcm, ycm, zcm, mt=0.0, 0.0, 0.0, 0.0
    for iatom in moleculein.atoms:
        si, xi, yi, zi=iatom.s, iatom.xc, iatom.yc, iatom.zc
        if si not in black_list:
            mi=get_atomic_mass(si)
            mt=mt+mi
            xcm=xcm+(xi*mi)
            ycm=ycm+(yi*mi)
            zcm=zcm+(zi*mi)
    vectorcm=np.array([xcm, ycm, zcm])/mt
    return vectorcm
#------------------------------------------------------------------------------------------
def translate(moleculein, vector, black_list=[]):
    for iatom in moleculein.atoms:
        if iatom.s not in black_list:
            iatom.xc=iatom.xc+vector[0]
            iatom.yc=iatom.yc+vector[1]
            iatom.zc=iatom.zc+vector[2]
    return moleculein
#------------------------------------------------------------------------------------------
def translate_to_cm(moleculein, black_list=[]):
    vcm=center_off_mass(moleculein, black_list)
    translate(moleculein, -vcm)
    return moleculein
#------------------------------------------------------------------------------------------
def translate_to_centroid(moleculein):
    vc=centroid(moleculein)
    translate(moleculein, -vc)
    return moleculein
#------------------------------------------------------------------------------------------
def rotate_matrix(moleculein, matrixr):
    matrix=np.copy(moleculein.m)
    if len(matrix) != 0:
        matrixtmp=np.transpose(matrix)
        matrixnew=np.matmul(matrixr,matrixtmp) 
        moleculein.m=np.transpose(matrixnew)   
    for iatom in moleculein.atoms:
        xc, yc, zc=iatom.xc, iatom.yc, iatom.zc
        vector=np.array([xc, yc, zc])
        vri = np.matmul(matrixr, vector)
        iatom.xc = vri[0]  
        iatom.yc = vri[1]
        iatom.zc = vri[2]
    return moleculein
#------------------------------------------------------------------------------------------
def rotate_deg(moleculein, qdegx, qdegy, qdegz):
    eulerm=euler_matrix(qdegx, qdegy, qdegz)
    rotate_matrix(moleculein, eulerm)
    return moleculein
#------------------------------------------------------------------------------------------
#vorig=center_off_mass(moleculein)
#vorig=centroid(moleculein)
def rotate_vector_angle_deg(moleculein, kvector, qdeg):
    rodriguesrm=rodrigues_rotation_matrix(kvector, qdeg)
    rotate_matrix(moleculein, rodriguesrm)
    return moleculein
#------------------------------------------------------------------------------------------
def random_deg_angles():
    qmin,qmax=0,360
    qdegx=random.randint(qmin,qmax)
    qdegy=random.randint(qmin,qmax)
    qdegz=random.randint(qmin,qmax)
    return qdegx, qdegy, qdegz
#------------------------------------------------------------------------------------------
def rotate_random(moleculein):
    qdegx, qdegy, qdegz=random_deg_angles()
    rotate_deg(moleculein, qdegx, qdegy, qdegz)
    return moleculein
#------------------------------------------------------------------------------------------
def get_inertia_tensor(moleculein, black_list=[]):
    moltmp=copymol(moleculein)
    translate_to_cm(moltmp, black_list)
    ixx,iyy,izz=float(0.0),float(0.0),float(0.0)  
    ixy,ixz,iyz=float(0.0),float(0.0),float(0.0)
    for iatom in moltmp.atoms:
        si, xi, yi, zi=iatom.s, iatom.xc, iatom.yc, iatom.zc
        if si not in black_list:
            mi=get_atomic_mass(si)
            ixx=mi*(yi*yi+zi*zi)+ixx
            iyy=mi*(xi*xi+zi*zi)+iyy
            izz=mi*(xi*xi+yi*yi)+izz
            ixy=mi*xi*yi+ixy
            ixz=mi*xi*zi+ixz
            iyz=mi*yi*zi+iyz
    inertia_tensor=np.array([[ixx,-ixy,-ixz], [-ixy,iyy,-iyz], [-ixz,-iyz,izz]])
    return inertia_tensor
#------------------------------------------------------------------------------------------
def get_principal_rotation_inertias(moleculein, opt=0, black_list=[]):
    """ SAMPLE OF USE IN PYTHON:
    opt=0 Three principal rotation inertias, Ia, Ib, and Ic, (Ia <= Ib <= Ic)
    opt=1 Three principal unitary vector axes, ua, ub, and uc.
    """
    InertialMatrix=get_inertia_tensor(moleculein, black_list)
    evals, evecs = np.linalg.eigh(InertialMatrix)
    vx, vy, vz=evecs[:,0], evecs[:,1], evecs[:,2]
    evecs=[vx,vy,vz]
    s = zip(evals, evecs)
    t = sorted(s, key=lambda x: float(x[0]))
    if opt==0:
        ia, ib, ic = t[0][0], t[1][0], t[2][0]
        return ia, ib, ic
    if opt==1:
        ua, ub, uc = t[0][1], t[1][1], t[2][1]
        vec= np.cross(ua, ub) - uc
        det= np.dot(vec,vec)
        if det > float(2.0): uc=-uc
        return ua, ub, uc
    if opt==2:
        ixy=InertialMatrix[0,1]
        ixz=InertialMatrix[0,2]
        iyz=InertialMatrix[1,2]
        norma_extradiagonal=np.sqrt(2.0*ixy*ixy+2.0*ixz*ixz+2.0*iyz*iyz)
        return norma_extradiagonal
#------------------------------------------------------------------------------------------
def align_all_inertia_axis_x(moleculelist, black_list=[]):
    for xmol in moleculelist:
        jmol,min=0,float(100.0)
        mol02=[]
        translate_to_cm(xmol)
        ua,ub,uc=get_principal_rotation_inertias(xmol, 1, black_list)
        list=([ua, ub, uc], [-ua, -ub, uc])
        for x, y, z in list:
            a=np.array([x,y,z])
            mol01=copymol(xmol)
            for iatom in mol01.atoms:
                vector=np.array([iatom.xc, iatom.yc, iatom.zc])
                vr=np.matmul(a,vector)
                iatom.xc=vr[0]
                iatom.yc=vr[1]
                iatom.zc=vr[2]
            ne=get_principal_rotation_inertias(mol01, 2, black_list)
            mol02.extend([mol01])
            if (float(ne) < min): min, im=ne, jmol
            jmol=jmol+1
        xmol.atoms=mol02[im].atoms
    return moleculelist
#------------------------------------------------------------------------------------------
def rotational_constants(moleculein, black_list=[]):
    ia,ib,ic=get_principal_rotation_inertias(moleculein, 0, black_list)
    h=6.62606957
    pi=3.141592654
    amu2kg=1.660538921
    factor=h*10000.0/(8.0*pi*pi*amu2kg)
    #Rotational constants in GHz
    A=round(factor/ia,7) if (ia >= 0.001) else 0.0
    B=round(factor/ib,7) if (ib >= 0.001) else 0.0
    C=round(factor/ic,7) if (ic >= 0.001) else 0.0
    return A,B,C
#------------------------------------------------------------------------------------------
def molecular_stoichiometry(moleculein, opt):
    """ SAMPLE OF USE IN PYTHON:
    >>> from utils.libmoleculas import readxyzs, molecular_stoichiometry
    >>> mol0=readxyzs('Ethanol.xyz')
    >>> print molecular_stoichiometry(mol0[0],0) # [('O', 1), ('C', 2), ('H', 6)]
    >>> print molecular_stoichiometry(mol0[0],1) # ['O', 'C', 'H']
    >>> print molecular_stoichiometry(mol0[0],2) # [1, 2, 6]
    >>> print molecular_stoichiometry(mol0[0],3) # O1 C2 H6
    >>> print molecular_stoichiometry(mol0[0],4) # ['O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    """
    natoms=moleculein.n
    if natoms==0:
        return {
            0: [],
            1: (),
            2: (),
            3: 0,
            4: []
        }[opt]
    else:
        #====================================================================
        listsym=list([x.s for x in moleculein.atoms])
        liste=list(set(listsym))
        listn=[listsym.count(x) for x in liste]
        listm=[get_atomic_mass(x) for x in liste]
        s = zip(liste, listn, listm)
        t = sorted(s, key=lambda x: float(x[2]), reverse=True)
        allatoms=([(x[0], x[1]) for x in t])
        atoms=([x[0] for x in t])
        nofatoms=([x[1] for x in t])
        auxlist=[]
        for item in allatoms:
            sti=item[0] if item[1]==1 else item[0]+str(item[1]) 
            auxlist.extend(sti)
        clustername=''.join(auxlist)
        inatoms=[]
        for iii in allatoms:
            for jjj in range(iii[1]):
                inatoms.append(iii[0])
        #====================================================================
        return {
            0: allatoms[:],
            1: atoms,
            2: nofatoms,
            3: clustername,
            4: inatoms
        }[opt]
#------------------------------------------------------------------------------------------
def sort_by_stoichiometry(moleculein,opt=0):
    """ SAMPLE OF USE IN PYTHON:
    >>> from utils.libmoleculas import readxyzs, sort_by_stoichiometry, writexyzs
    >>> mol0=readxyzs('Ethanol.xyz')
    >>> mol1=sort_by_stoichiometry(mol0[0])
    >>> writexyzs([mol1],'w.xyz')
    """
    moleculeout=[]
    for imol in moleculein:
        listo=[]
        lists=list([imol.s for imol in imol.atoms])
        all_atoms=molecular_stoichiometry(imol,4)
        for symr in all_atoms:
            i=0
            for symu in lists:
                if (symr == symu) and (i not in listo): listo.append(i)
                i=i+1
        if opt==1: listo.reverse()
        xmol=Molecule(imol.i, imol.e, imol.m)
        for j in listo:
            xmol.add_atom(imol.atoms[j])
        moleculeout.extend([xmol])
    return moleculeout
#------------------------------------------------------------------------------------------
def scale_coords(moleculein, factor):
    moleculeout=copymol(moleculein)
    translate_to_cm(moleculeout)
    for iatom in moleculeout.atoms:
        iatom.xc=(iatom.xc)*factor
        iatom.yc=(iatom.yc)*factor
        iatom.zc=(iatom.zc)*factor
    return moleculeout
#------------------------------------------------------------------------------------------
def molecular_radius(moleculein):
    moleculeout=copymol(moleculein)
    translate_to_cm(moleculeout)
    r=[]
    for iatom in moleculeout.atoms:
        xyz=np.array([iatom.xc, iatom.yc, iatom.zc])
        ri=get_covalent_radius(iatom.s)
        dist=np.linalg.norm(xyz)+ri
        r.append(dist)
    r.sort()
    rmin=r[0]
    ravg=np.mean(r)
    rmax=r[-1]
    moleculeout.clear()
    return rmin,ravg,rmax
#------------------------------------------------------------------------------------------
def molecular_distance(moleculealone, moleculelist, opt=0):
    dmin=float(100.0)
    for iatom in moleculealone.atoms:
        vi=np.array([iatom.xc, iatom.yc, iatom.zc])
        ri=get_covalent_radius(iatom.s)
        for xmol in moleculelist:
            for jatom in xmol.atoms:
                vj=np.array([jatom.xc, jatom.yc, jatom.zc])
                rj=get_covalent_radius(jatom.s)
                dist=np.linalg.norm(vi-vj)
                dist=dist/(ri+rj) if opt==0 else dist
                if( dist < dmin): dmin=dist
    return dmin
#------------------------------------------------------------------------------------------
def atom_atom_dist(iatom, jatom):
    vi=np.array([iatom.xc, iatom.yc, iatom.zc])
    vj=np.array([jatom.xc, jatom.yc, jatom.zc])
    dist=np.linalg.norm(vi-vj)
    return dist
#------------------------------------------------------------------------------------------
def atom_molecule_dist(atomin, moleculein, opt=0):
    dmin=float(100.0)
    rr=get_covalent_radius(atomin.s)
    vv=np.array([atomin.xc, atomin.yc, atomin.zc])
    for iatom in moleculein.atoms:
        ri=get_covalent_radius(iatom.s)
        vi=np.array([iatom.xc, iatom.yc, iatom.zc])
        dist=np.linalg.norm(vi-vv)
        dist=dist/(ri+rr) if opt==0 else dist
        if( dist < dmin): dmin=dist
    return dmin
#------------------------------------------------------------------------------------------
def get_min_binding_distance(moleculein, opt=0):
    dmin=float(100.0)
    natoms=moleculein.n
    for iatom in range(natoms):
        si=moleculein.atoms[iatom].s
        ri=get_covalent_radius(si)
        xi=moleculein.atoms[iatom].xc
        yi=moleculein.atoms[iatom].yc
        zi=moleculein.atoms[iatom].zc
        vi=np.array([xi, yi, zi])
        for jatom in range(iatom+1,natoms):
            sj=moleculein.atoms[jatom].s
            rj=get_covalent_radius(sj)
            xj=moleculein.atoms[jatom].xc
            yj=moleculein.atoms[jatom].yc
            zj=moleculein.atoms[jatom].zc
            vj=np.array([xj, yj, zj])
            dist=np.linalg.norm(vj-vi)
            dist=dist/(ri+rj) if opt==0 else dist
            if( dist < dmin): dmin=dist
    return dmin
#------------------------------------------------------------------------------------------
def get_from_list(moleculein, list=[]):
    moleculeout=[]
    total_molecules=len(moleculein)
    for ii in list:
        if ii < total_molecules:
            moleculeout.extend([moleculein[ii]])
        else:
            print('Error in get_from_list: bad index')
            exit()
    return moleculeout
#------------------------------------------------------------------------------------------
def run_sample():
    from inout.molecules import write_molecule
    write_molecule('coenzymea')
    mol0=readxyzs('coenzymea.xyz')
    mol1=copymol(mol0[0])
    print(molecular_stoichiometry(mol1,0))
    print(molecular_stoichiometry(mol1,1))
    print(molecular_stoichiometry(mol1,2))
    print(molecular_stoichiometry(mol1,3))
    print(molecular_stoichiometry(mol1,4))
    mol2=sort_by_stoichiometry(mol0)
    writexyzs(mol2,'pelot.xyz')
    print(center_off_mass(mol2[0]))
    writexyzs(mol2,'pelou.xyz')
    vec=np.array([0, 3, 4])    
    translate(mol2[0], vec)
    writexyzs(mol2,'pelov.xyz')
    print(molecular_distance(mol1,mol2))
    print(get_min_binding_distance(mol2[0]))
    print(molecular_radius(mol2[0]))
    print(get_principal_rotation_inertias(mol2[0],0))
    print(get_principal_rotation_inertias(mol2[0],1))
    print(get_principal_rotation_inertias(mol2[0],2))
    mol3=rename_molecule(mol2,'hola',3)
    writexyzs(mol3,'pelow.xyz')
#run_sample()
#------------------------------------------------------------------------------------------
