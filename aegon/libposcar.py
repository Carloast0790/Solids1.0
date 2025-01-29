import os
import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from aegon.libutils import align
from collections import defaultdict
#------------------------------------------------------------------------------------------
def make_matrix(singlemoleculein, latsp):
    mol0=singlemoleculein.copy()
    mol0.translate(-mol0.get_center_of_mass())
    listx=[abs(iatom.position[0]) + covalent_radii[iatom.number] for iatom in mol0]
    listy=[abs(iatom.position[1]) + covalent_radii[iatom.number] for iatom in mol0]
    listz=[abs(iatom.position[2]) + covalent_radii[iatom.number] for iatom in mol0]
    diamx=max(listx)*2.0+latsp
    diamy=max(listy)*2.0+latsp
    diamz=max(listz)*2.0+latsp
    matrix=np.array([[diamx, 0.0, 0.0],[0.0, diamy, 0.0],[0.0, 0.0, diamz]])
    return matrix
#------------------------------------------------------------------------------------------
def molecule2poscar(atomslist, latsp=5.0):
    moleculeout=[]
    for iatoms in atomslist:
        singleposcar=iatoms.copy()
        align(singleposcar)
        singleposcar.translate(-singleposcar.get_center_of_mass())
        matrix=make_matrix(singleposcar, latsp)
        vc=np.array([matrix[0,0], matrix[1,1], matrix[2,2]])/float(2.0)
        singleposcar.translate(+vc)
        singleposcar.set_cell(matrix)
        singleposcar.set_pbc(True)
        moleculeout.extend([singleposcar])
    return moleculeout
#------------------------------------------------------------------------------------------
def readposcars(filename):
    if not os.path.isfile(filename):
        print("The file %s does not exist." %(filename))
        exit()
    contcarfile=open(filename,'r')
    poscarout=[]
    for line in contcarfile:
        header=line.split()
        name=header[0]
        energy=float(header[1]) if len(header)>1 else 0.0
        line=contcarfile.readline()
        if str(line.split()[0])=='0.00000000E+00': break
        #-----------------------------------
        scale=float(line.split()[0])
        line=contcarfile.readline()
        a1x, a1y, a1z=map(float,line.split())
        line=contcarfile.readline()
        a2x, a2y, a2z=map(float,line.split())
        line=contcarfile.readline()
        a3x, a3y, a3z=map(float,line.split())
        lattice_vectors=np.array([[a1x, a1y, a1z],[a2x, a2y, a2z],[a3x, a3y, a3z]])*scale
        #-----------------------------------
        line=contcarfile.readline()
        elements=line.split()
        line=contcarfile.readline()
        ocupnumchar=line.split()
        ocupnuminte=list(map(int, ocupnumchar))
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
        symbols=[]
        positions=[]
        if 'Direct' in line:
            for iatom in range(natom):
                line=contcarfile.readline()
                vecxyz=line.split()
                s=liste[iatom]
                symbols.append(s)
                xd=float(vecxyz[0])
                yd=float(vecxyz[1])
                zd=float(vecxyz[2])
                direct_coords = np.array([xd, yd, zd])
                cart_coords = np.dot(direct_coords, lattice_vectors)
                positions.append(cart_coords)
        if 'Cartesian' in line:
            for iatom in range(natom):
                s=liste[iatom]
                symbols.append(s)
                line=contcarfile.readline()
                cart_coords =np.array(line.split()[0:3])
                positions.append(cart_coords)
        positions=np.array(positions)
        iposcar = Atoms(symbols=symbols, positions=positions, cell=lattice_vectors, pbc=True)
        iposcar.info['i']=name
        iposcar.info['e']=energy
        poscarout.extend([iposcar])
    contcarfile.close()
    return poscarout
#------------------------------------------------------------------------------------------
def writeposcars(poscarlist, file, opt='D', in_log=0):
    f=open(file,"w")
    for atoms in poscarlist:
        print("%s %12.8f" %(atoms.info['i'], atoms.info['e']), file=f)
        print('1.0', file=f)
        matrix=atoms.cell
        print("%20.16f %20.16f %20.16f" %(matrix[0,0],matrix[0,1],matrix[0,2]), file=f)
        print("%20.16f %20.16f %20.16f" %(matrix[1,0],matrix[1,1],matrix[1,2]), file=f)
        print("%20.16f %20.16f %20.16f" %(matrix[2,0],matrix[2,1],matrix[2,2]), file=f)
        element_count = defaultdict(int)
        for symbol in atoms.symbols: element_count[symbol] += 1
        print(' '.join([str(item) for item in element_count.keys()]), file=f)
        print(' '.join([str(item) for item in element_count.values()]), file=f)
        element_count = defaultdict(int)
        if opt=='C':
            print('Cartesian', file=f)
            for atom in atoms:
                symbol=atom.symbol
                element_count[symbol] += 1
                xc, yc, zc = atom.position
                print("%20.16f %20.16f %20.16f   !%s%d" %(xc,yc,zc,symbol,element_count[symbol]), file=f)
        if opt=='D':
            print('Direct', file=f)
            mi=np.linalg.inv(matrix)
            for atom in atoms:
                symbol=atom.symbol
                element_count[symbol] += 1
                cart_coords=np.array(atom.position)
                ##CHECAR PARA SOLIDS. LA T.
                vd=np.dot(mi.T, cart_coords)
                xd, yd, zd = vd[0], vd[1], vd[2]
                print("%20.16f %20.16f %20.16f   !%s%d" %(xd,yd,zd,symbol,element_count[symbol]), file=f)
    f.close()
    if in_log==0: print("Writing %s" %(file))
#------------------------------------------------------------------------------------------
def writecfg(poscarlist, file):
    f=open(file,"w")
    for atoms in poscarlist:
        print("BEGIN_CFG", file=f)
        print(" Size", file=f)
        natoms=len(atoms)
        print("    %d" %(natoms), file=f)
        print(" Supercell", file=f)
        matrix=atoms.cell
        mi=np.linalg.inv(matrix)
        print("       %10.6f    %10.6f    %10.6f" %(matrix[0,0],matrix[0,1],matrix[0,2]), file=f)
        print("       %10.6f    %10.6f    %10.6f" %(matrix[1,0],matrix[1,1],matrix[1,2]), file=f)
        print("       %10.6f    %10.6f    %10.6f" %(matrix[2,0],matrix[2,1],matrix[2,2]), file=f)
        print(" AtomData:  id type       cartes_x      cartes_y      cartes_z", file=f)
        for k, atom in enumerate(atoms):
            xc, yc, zc = atom.position
            print("             0    %d     %10.6f    %10.6f    %10.6f" %(k+1,xc,yc,zc), file=f)
        print("END_CFG", file=f)
    f.close()
#------------------------------------------------------------------------------------------
def sort_poscar(atoms):
    atoms_copy = atoms.copy()
    species_positions = defaultdict(list)
    for position, species in zip(atoms.positions, atoms.symbols):
        species_positions[species].append(position)
    k=0
    for specie in species_positions.keys():
        for position in species_positions[specie]:
            atoms_copy[k].symbol=specie
            atoms_copy[k].position=position
            k=k+1
    return atoms_copy
#------------------------------------------------------------------------------------------
def conventional(atoms, tol=0.01):
    lattice_vectors=np.copy(atoms.cell)
    xdmin,xdmax=0.0,1.0
    ydmin,ydmax=0.0,1.0
    zdmin,zdmax=0.0,1.0
    direct_coords = atoms.get_scaled_positions()
    atom_symbols = atoms.get_chemical_symbols()
    symbols, positions=[], []
    for i, icoord in enumerate(direct_coords):
        xd, yd, zd=icoord
        for x in [-1,0,1]:
            for y in [-1,0,1]:
                for z in [-1,0,1]:
                    xd2=xd+x
                    yd2=yd+y
                    zd2=zd+z
                    cond1=(xd2 >= xdmin-tol) and (xd2 <= xdmax+tol)
                    cond2=(yd2 >= ydmin-tol) and (yd2 <= ydmax+tol)
                    cond3=(zd2 >= zdmin-tol) and (zd2 <= zdmax+tol) 
                    if cond1 and cond2 and cond3:
                        symbols.append(atom_symbols[i])
                        vectord=np.array([xd2, yd2, zd2])
                        cart_coords=np.dot(vectord, lattice_vectors)
                        positions.append(cart_coords)
    poscarout = Atoms(symbols=symbols, positions=positions, cell=lattice_vectors, pbc=True)
    #poscarout.info['i']=atoms.info['i']
    #poscarout.info['e']=atoms.info['e']
    return poscarout
#------------------------------------------------------------------------------------------