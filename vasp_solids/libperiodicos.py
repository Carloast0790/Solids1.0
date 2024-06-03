import os.path
import numpy as np
from utils_solids.atomic import get_covalent_radius
from utils_solids.libmoleculas import Atom, Molecule, copymol, translate_to_cm, translate, align_all_inertia_axis_x
from inout_solids.getbilparam  import get_a_str, get_a_float
#------------------------------------------------------------------------------------------
def direct2cartesian(xd, yd, zd, matrix):
    #vector=np.array([xd, yd, zd])
    #vd=np.inner(matrix,vector)
    #xc, yc, zc = vd[0], vd[1], vd[2]
    a1x,a1y,a1z=matrix[0,0],matrix[0,1],matrix[0,2]
    a2x,a2y,a2z=matrix[1,0],matrix[1,1],matrix[1,2]
    a3x,a3y,a3z=matrix[2,0],matrix[2,1],matrix[2,2]
    xc=(a1x*xd+a2x*yd+a3x*zd)
    yc=(a1y*xd+a2y*yd+a3y*zd)
    zc=(a1z*xd+a2z*yd+a3z*zd)
    return xc, yc, zc
#------------------------------------------------------------------------------------------
def cartesian2direct(xc, yc, zc, matrix):
    #mi=np.linalg.inv(matrix)
    #vector=np.array([xc, yc, zc])
    #vd=np.inner(mi,vector)
    #xd, yd, zd = vd[0], vd[1], vd[2]
    a1x,a1y,a1z=matrix[0,0],matrix[0,1],matrix[0,2]
    a2x,a2y,a2z=matrix[1,0],matrix[1,1],matrix[1,2]
    a3x,a3y,a3z=matrix[2,0],matrix[2,1],matrix[2,2]
    det=(a1z*a2y*a3x - a1y*a2z*a3x - a1z*a2x*a3y + a1x*a2z*a3y + a1y*a2x*a3z - a1x*a2y*a3z)
    xd=(a2z*a3y*xc - a2y*a3z*xc - a2z*a3x*yc + a2x*a3z*yc + a2y*a3x*zc - a2x*a3y*zc)/det
    yd=(a1y*a3z*xc - a1z*a3y*xc + a1z*a3x*yc - a1x*a3z*yc - a1y*a3x*zc + a1x*a3y*zc)/det
    zd=(a1z*a2y*xc - a1y*a2z*xc - a1z*a2x*yc + a1x*a2z*yc + a1y*a2x*zc - a1x*a2y*zc)/det
    if zd < 0.0:
       zd=zd+1.0
    if zd > 1.0:
        zd=zd-1.0
    return xd, yd, zd
#------------------------------------------------------------------------------------------
def readposcars(filename):
    if not os.path.isfile(filename):
        print("The file %s does not exist." %(filename))
        exit()
    contcarfile=open(filename,'r')
    lines = contcarfile.readlines()
    last = lines[-1]
    contcarfile.close()
    contcarfile=open(filename,'r')
    line=contcarfile.readline()
    poscarout=[]
    while(line != '\n'):
        name=line.strip()
        #name=name.replace(" ", "")
        energy=float(0.0)
        line=contcarfile.readline()
        if str(line.split()[0])=='0.00000000E+00': break
        #-----------------------------------
        a0=float(line.split()[0])
        line=contcarfile.readline()
        a1x, a1y, a1z=map(float,line.split())
        line=contcarfile.readline()
        a2x, a2y, a2z=map(float,line.split())
        line=contcarfile.readline()
        a3x, a3y, a3z=map(float,line.split())
        matrix=np.array([[a0*a1x, a0*a1y, a0*a1z],[a0*a2x, a0*a2y, a0*a2z],[a0*a3x, a0*a3y, a0*a3z]])
        poscarx=Molecule(name, energy, matrix)
        #-----------------------------------
        line=contcarfile.readline()
        elements=line.split()
        line=contcarfile.readline()
        ocupnumchar=line.split()
        ocupnuminte=list(map(int, ocupnumchar))
        #chemformula=zip(elements,ocupnuminte)
        #-----------------------------------
        natom=sum(ocupnuminte)
        liste,kk=[],0
        for ii in ocupnuminte:
            for jj in range(ii):
                liste.append(elements[kk])
            kk=kk+1
        #-----------------------------------
        line=contcarfile.readline()
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
        poscarout.extend([poscarx])
        line='\n' if line == last else contcarfile.readline()
    contcarfile.close()
    return poscarout
#------------------------------------------------------------------------------------------
def writeposcars(poscarlist, file, opt='D'):
    fh=open(file,"w")
    for poscar in poscarlist:
        print(poscar.i, file=fh)
        print('1.0', file=fh)
        matrix=poscar.m
        print("%20.16f %20.16f %20.16f" %(matrix[0,0],matrix[0,1],matrix[0,2]), file=fh)
        print("%20.16f %20.16f %20.16f" %(matrix[1,0],matrix[1,1],matrix[1,2]), file=fh)
        print("%20.16f %20.16f %20.16f" %(matrix[2,0],matrix[2,1],matrix[2,2]), file=fh)
        ste=chem_formula(poscar)
        print(' '.join([str(item[0]) for item in ste]), file=fh)
        print(' '.join([str(item[1]) for item in ste]), file=fh)
        print('Selective dynamics', file=fh)
        if opt=='D':
            print('Direct', file=fh)
            ii=1
            for iatom in poscar.atoms:
                s,xc,yc,zc,xf,yf,zf=iatom.s,iatom.xc,iatom.yc,iatom.zc,iatom.xf,iatom.yf,iatom.zf
                xd, yd, zd=cartesian2direct(xc, yc, zc, poscar.m)
                print("%20.16f %20.16f %20.16f %3s %3s %3s !%s%d" %(xd,yd,zd,xf,yf,zf,s,ii), file=fh)
                ii=ii+1
        if opt=='C':
            print('Cartesian', file=fh)
            ii=1
            for iatom in poscar.atoms:
                s,xc,yc,zc,xf,yf,zf=iatom.s,iatom.xc,iatom.yc,iatom.zc,iatom.xf,iatom.yf,iatom.zf
                print("%20.16f %20.16f %20.16f %3s %3s %3s !%s%d" %(xc,yc,zc,xf,yf,zf,s,ii), file=fh)
                ii=ii+1
    fh.close()
    log_file =get_a_str('output_file','solids_out.txt')
    fopen = open(log_file,'a')
    print("Created: %s" %(file), file=fopen)
    fopen.close()
#------------------------------------------------------------------------------------------
def center_off_cell(singleposcar):
    matrix=singleposcar.m
    a1=np.array(matrix[0,:])
    a2=np.array(matrix[1,:])
    a3=np.array(matrix[2,:])
    cc=(a1+a2+a3)/2.0
    return cc
#------------------------------------------------------------------------------------------
def unfix_nfirst_slabs(poscarlist, nfixslab):
    for poscar in poscarlist:
        n=poscar.n
        zlist, ilist= [], []
        for j in range(n):
            zi=poscar.atoms[j].zc
            zlist.append(zi)
            ilist.append(j)
        s = zip(ilist, zlist)
        t = sorted(s, key=lambda x: float(x[1]))
        t.reverse()
        index=t[0][0]
        z0=poscar.atoms[index].zc
        nslab=1
        for j in range(n):
            index=t[j][0]
            zc=poscar.atoms[index].zc
            nslab=nslab if zc==z0 else nslab+1
            if nslab > nfixslab:
                poscar.atoms[index].xf='F'
                poscar.atoms[index].yf='F'
                poscar.atoms[index].zf='F'
            else:
                poscar.atoms[index].xf='T'
                poscar.atoms[index].yf='T'
                poscar.atoms[index].zf='T'
            if zc<z0: z0=zc
        print('Max number of slabs = %d' %(nslab))
    return poscarlist
#------------------------------------------------------------------------------------------
##CHECAR SI ES INDISPENSABLE chem_formula
#------------------------------------------------------------------------------------------
def chem_formula(poscar):
    listsym, xlist, ylist=[], [], []
    for iatom in poscar.atoms: listsym.append(iatom.s)
    ord=-1
    for ii in listsym:
        if ii not in xlist:
            xlist.append(ii)
            count=1
            ylist.append(count)
            ord=ord+1
        else:
            count=count+1
        ylist[ord]=count
    zlist = list(zip(xlist, ylist))
    return zlist[:]
#------------------------------------------------------------------------------------------
def expand_poscar(poscarlist, xff=2, yff=2, zff=2, xii=0, yii=0, zii=0):
    poscarout=[]
    for poscarx in poscarlist:
        matrix=poscarx.m
        a1, a2, a3=matrix[0,:], matrix[1,:], matrix[2,:]
        matrixf=np.array([(xff-xii)*a1,(yff-yii)*a2,(zff-zii)*a3])
        poscary=Molecule(poscarx.i, poscarx.e, matrixf)
        comment=[]
        for ii,iatom in enumerate(poscarx.atoms):
            for x in range(xii,xff,1):
                for y in range(yii,yff,1):
                    for z in range(zii,zff,1):
                        vt=float(x)*a1+float(y)*a2+float(z)*a3
                        xc, yc, zc=iatom.xc+vt[0], iatom.yc+vt[1], iatom.zc+vt[2]
                        #if x==1 and y==1 and z==1:
                        #    xf,yf,zf='T','T','T'
                        #else:
                        #    xf,yf,zf='F','F','F'
                        xf,yf,zf='T','T','T'
                        ai=Atom(iatom.s,xc,yc,zc,xf,yf,zf)
                        poscary.add_atom(ai)
                        comment.append(ii)
        poscary.c=comment
        poscarout.extend([poscary])
    return poscarout
#------------------------------------------------------------------------------------------
def expand_poscar3x3(singleposcar):
    matrix=singleposcar.m
    a1, a2, a3=matrix[0,:], matrix[1,:], matrix[2,:]
    matrixf=np.array([3*a1,3*a2,3*a3])
    poscarout=Molecule(singleposcar.i, singleposcar.e, matrixf)
    for ii,iatom in enumerate(singleposcar.atoms):
        for x in [-1,0,1]:
            for y in [-1,0,1]:
                for z in [-1,0,1]:
                    vt=float(x)*a1+float(y)*a2+float(z)*a3
                    xc, yc, zc=iatom.xc+vt[0], iatom.yc+vt[1], iatom.zc+vt[2]
                    if (x==0) and (y==0) and (z==0):
                        xf, yf, zf=iatom.xf, iatom.yf, iatom.zf
                    else:
                        xf, yf, zf='F','F','F'
                    ai=Atom(iatom.s,xc,yc,zc,xf,yf,zf)
                    poscarout.add_atom(ai)
    return poscarout
#------------------------------------------------------------------------------------------
def conventional(singleposcar, tol=0.01):
    matrix=np.copy(singleposcar.m)
    a1, a2, a3=matrix[0,:], matrix[1,:], matrix[2,:]
    mi=np.linalg.inv(matrix)
    poscarout=Molecule(singleposcar.i, singleposcar.e, matrix)
    xdmin,xdmax=0.0,1.0
    ydmin,ydmax=0.0,1.0
    zdmin,zdmax=0.0,1.0
    for ii,iatom in enumerate(singleposcar.atoms):
        for x in [-1,0,1]:
            for y in [-1,0,1]:
                for z in [-1,0,1]:
                    vt=float(x)*a1+float(y)*a2+float(z)*a3
                    xc, yc, zc=iatom.xc+vt[0], iatom.yc+vt[1], iatom.zc+vt[2]
                    xf, yf, zf=iatom.xf, iatom.yf, iatom.zf
                    vector=np.array([xc, yc, zc])
                    vd=np.matmul(mi,vector)
                    xd, yd, zd = vd[0], vd[1], vd[2]
                    suma=0
                    if (xd >= xdmin-tol) and (xd <= xdmax+tol): suma=suma+1 
                    if (yd >= ydmin-tol) and (yd <= ydmax+tol): suma=suma+1 
                    #if (zd >= zdmin-tol) and (zd <= zdmax+tol): suma=suma+1 
                    if (zd >= zdmin) and (zd <= zdmax): suma=suma+1 
                    if suma == 3:
                        ai=Atom(iatom.s,xc,yc,zc,xf,yf,zf)
                        poscarout.add_atom(ai)
    return poscarout
#------------------------------------------------------------------------------------------
def make_matrix(singlemoleculein):
    mol0=copymol(singlemoleculein)
    translate_to_cm(mol0)
    listx=[iatom.xc + get_covalent_radius(iatom.s) for iatom in mol0.atoms]
    listy=[iatom.yc + get_covalent_radius(iatom.s) for iatom in mol0.atoms]
    listz=[iatom.zc + get_covalent_radius(iatom.s) for iatom in mol0.atoms]
    latsp=get_a_float('latt_space', 5.0)
    diamx=np.abs(max(listx))+np.abs(min(listx))+latsp
    diamy=np.abs(max(listy))+np.abs(min(listy))+latsp
    diamz=np.abs(max(listz))+np.abs(min(listz))+latsp
    ##print(mol0.i, latsp, diamx, diamy, diamz)
    matrix=np.array([[diamx, float(0.0), float(0.0)],[float(0.0), diamy, float(0.0)],[float(0.0), float(0.0), diamz]])
    return matrix
#------------------------------------------------------------------------------------------
def molecule2poscar(singlemoleculein):
    singleposcar=copymol(singlemoleculein)
    align_all_inertia_axis_x([singleposcar])
    translate_to_cm(singleposcar)
    matrix=make_matrix(singleposcar)
    vc=np.array([matrix[0,0], matrix[1,1], matrix[2,2]])/float(2.0)
    translate(singleposcar,+vc)
    singleposcar.m=matrix
    return singleposcar
#------------------------------------------------------------------------------------------
def run_sample():
    from vasp.poscars import write_poscar_examples
    write_poscar_examples('Febcc')
    poscar1=readposcars('Febcc.vasp')
    poscar2=copymol(poscar1[0])
    unfix_nfirst_slabs([poscar2], 1)
    print(chem_formula(poscar2))
    poscar3=expand_poscar([poscar2],3,3,3)
    writeposcars(poscar3,'zposcar.vasp')
    poscar4=conventional(poscar2)
    poscar5=expand_poscar([poscar4],2,2,2)
    writeposcars(poscar5,'xposcar.vasp')
#run_sample()
#------------------------------------------------------------------------------------------
