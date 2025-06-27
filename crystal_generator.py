import random
import os.path
from pyxtal import pyxtal
from pyxtal.symmetry import Group
from pyxtal.tolerance import Tol_matrix
from pyxtal.lattice import Lattice
from ase.data import covalent_radii, atomic_numbers
from ase.io import write
from aegon.libstdio import read_main_input
import warnings
warnings.filterwarnings("ignore")

file = 'INPUT.txt'
log_file = 'solids_out.txt'

df = read_main_input(file)
df_composition = df.get_comp(key='COMPOSITION')
species = [k[0] for k in df_composition.comp]
formula_units = df.get_int('formula_units',2)
atms_per_specie = [k[1]*formula_units for k in df_composition.comp]
dimension = df.get_int('dimension',3)
vol_factor = df.get_float('volume_factor',1.0)
revisit_syms = df.get_int('revisit_syms',1)

#------------------------------------------------------------------------------------------------
'''
Unit cell constrictions 
'''
def uc_restriction():
    '''Obtains the UC restrictions imposed by the user. Format a b c alpha beta gamma.

    out: 
    restr_uc (Lattice); Pyxtal object for the lattice
    '''
    if os.path.isfile(file):
        restr_lattice = False
        f=open(file,"r")
        for line in f:
            if not line.startswith('#') and 'fixed_lattice' in line:
                line = line.split()
                a,b,c,alpha,beta,gamma = float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6])
                restr_lattice = Lattice.from_para(a,b,c,alpha,beta,gamma)
                break
        f.close()
    return restr_lattice
#------------------------------------------------------------------------------------------------
'''
Interatomic distances
'''
def get_percent_tolerances():
    '''Gets the default tolerances for each pair of atoms in the structure
    
    out: 
    solids_tolerances (list), List containing tuples with each min int dist, [(s1,s2,d1),(s1,s3,d2),...]
    pyxtal_tolerances (Tol_matrix), PyXtal object used for atomic tolerance in generation of structures
    '''
    pyxtal_tolerances = Tol_matrix()
    solids_tolerances = []
    species_number = [atomic_numbers[s] for s in species]
    if len(species) == 1:
        s = species[0]
        n = species_number[0]
        r = covalent_radii[n]
        tv = r*tolerance_percent*2
        tv = round(tv,2)
        solids_tolerances.append((s,s,tv))
        pyxtal_tolerances.set_tol(s,s,tv)
    else:
        for i in range(len(species)):
            s1 = species[i]
            r1 = covalent_radii[species_number[i]]
            tv = r1*tolerance_percent
            tv = round(tv,2)
            solids_tolerances.append((s1,s1,tv))
            pyxtal_tolerances.set_tol(s1,s1,tv)
            for j in range(i+1,len(species)):
                s2 = species[j]
                r2 = covalent_radii[species_number[j]]
                tv_mix = (r1+r2)*tolerance_percent
                tv_mix = round(tv_mix,2)
                solids_tolerances.append((s1,s2,tv_mix))
                pyxtal_tolerances.set_tol(s1,s2,tv_mix)
    return solids_tolerances,pyxtal_tolerances

#------------------------------------------------------------------------------------------------
def interatom_restriction():
    pyxtal_tolerances = Tol_matrix()
    solids_tolerances = []
    if os.path.isfile(file):
        xfile = open(file,"r")
        for line in xfile:
            if not line.startswith('#') and 'custom_tolerances' in line:
                readline = line.split()
                readline = readline[1:]
                for i in readline:
                    x = i.split(',')
                    tupla = (x[0],x[1],float(x[2]))
                    solids_tolerances.append(tupla)
                    pyxtal_tolerances.set_tol(x[0],x[1],float(x[2]))
            elif not line.startswith('#') and 'tol_atomic_overlap' in line:
                solids_tolerances, pyxtal_tolerances = get_percent_tolerances()
        xfile.close()
    return solids_tolerances, pyxtal_tolerances

#------------------------------------------------------------------------------------------------
'''
Restrictions on symmetry
'''
def get_symmetry_constrains(str_range, dimension=3):
    ''' This routine extracts a desired range of integers to be used as SGs in the construction of
    structures. The result is presented in list format, eg. range 1-5, range_list = [1,2,3,4,5]. If 
    the restriction is not provided, the list ranges from 2-80 for 2D structures and from 2-230 for 3D.

    in: str_range (str), flag to locate the desired range of integers
        dimension (int), list of all numbers within the desired range
    out: range_list (list), list of all numbers within the desired range
    '''
    import os.path
    file = 'INPUT.txt'
    range_list = []
    if os.path.isfile(file):
        f = open(file,'r')
        flag = False
        for line in f:
            if not line.startswith('#') and str_range in line:
                line = line.lstrip('\t\n\r')
                line = line.split()
                readline = line[1].split('-')
                bottom, top = int(readline[0])-1, int(readline[1])
                range_list = [s+1 for s in range(bottom,top)]
                flag = True
                break
        f.close()
        if flag == False and dimension == 2:
            range_list = [i for i in range(2,81)]
        elif flag == False and dimension == 3:
            range_list = [i for i in range(2,231)]
    return range_list

#------------------------------------------------------------------------------------------------
def random_crystal_gen_SM(species,atoms_per_specie,dimension,volume_factor,pyxtal_mtx_tolerance,uc_restr,revisit_syms):
    '''Generates a specified number of random crystals. It considers symmetry and dimension constrictions.
    If "revisit_syms" is provided, it rebuilds all random structures revisiting the same space groups. 

    in:
    species (list); This list contains each atomic symbol of the chemical composition 
    atoms_per_specie (list); This list contains the number of atoms corresponding to each species
    dimension (int); Dimension of the structure, right now only 3D
    volume_factor (float); Escalling factor for the unit cell
    pyxtal_mtx_tolerance (Tol_matrix); pyxtal.tolerance Tol_matrix object that contains the atomico overlap permited in 
        the construction fo the structure
    uc_restr (Pyxtal's Lattice or bool); If the user provided cell restriction (a,b,c.alpha,beta,gamma),
        this is Pyxtal's Lattice. False, otherwise.
    revisit_syms (int); The number of times space gorups are used to build random structures    

    out:
    xtalist_out (list); This list will contain all crystal structures as Molecule object
    '''
    xtalist_out = []
    fopen = open(log_file,'a')
    xc = 1
    sym_list = get_symmetry_constrains('symmetries', dimension)
    for i in range(revisit_syms):
        for sym in sym_list:
            xtal = pyxtal()
            if dimension == 2:
                try:
                    xtal.from_random(dimension,sym,species,atoms_per_specie,thickness=0.0)
                except:
                    continue
                else:
                    sg = Group (sym)
                    sg_symbol = str(sg.symbol)
                    ase_xtal = xtal.to_ase()
                    print('random_000_'+str(xc).zfill(3)+' ---> SG_'+str(sg_symbol)+"_("+str(sym)+")",file=fopen)
                    xtalist_out.append(ase_xtal)
                    xc = xc + 1
            elif dimension == 3:
                try:
                    if uc_restr:
                        xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,pyxtal_mtx_tolerance, lattice=uc_restr)
                    else:
                        xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,pyxtal_mtx_tolerance)
                except:
                    continue
                else:
                    sg = Group (sym)
                    sg_symbol = str(sg.symbol)
                    ase_xtal = xtal.to_ase()
                    print('random_000_'+str(xc).zfill(3)+' ---> SG_'+str(sg_symbol)+"_("+str(sym)+")",file=fopen)
                    xtalist_out.append(ase_xtal)
                    xc = xc + 1
    fopen.close()
    return xtalist_out

tolerance_percent = df.get_float(key='tol_atomic_overlap',default=0.97)
solids_mtx_tolerance, pyxtal_mtx_tolerance = interatom_restriction()
unit_cell_restriction = uc_restriction()

x = random_crystal_gen_SM(species,atms_per_specie,dimension,vol_factor,pyxtal_mtx_tolerance,unit_cell_restriction,revisit_syms)
write('POSCAR',x[0],vasp5=True, direct=True, append=False)