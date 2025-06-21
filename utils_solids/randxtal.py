import random
from pyxtal import pyxtal
from pyxtal.symmetry import Group
from utils_solids.libmoleculas import sort_by_stoichiometry, rename_molecule
from utils_solids.miscellaneous import pyxtal2solids, rescale_str, unit_cell_non_negative_coordinates,get_symmetry_constrains, vol_restriction
from inout_solids.getbilparam import get_a_str, get_a_int
import warnings
warnings.filterwarnings("ignore")
#------------------------------------------------------------------------------------------------
# Variables
log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------------
def random_crystal_gen_SM(revisit_syms,species,atoms_per_specie,p_list,formula_units=1,dimension=3,volume_factor=1.0, uc_restr=False):
    '''Generates a specified number of random crystals. It considers symmetry and dimension constrictions.
    If "revisit_syms" is provided, it rebuilds all random structures revisiting the same space groups. 

    in:
    revisit_syms (int); The number of times space gorups are used to build random structures
    species (list); This list contains each atomic symbol of the chemical composition 
    atoms_per_specie (list); This list contains the number of atoms corresponding to each species
    p_list (Tol_matrix); pyxtal.tolerance Tol_matrix object that contains the atomico overlap permited in 
        the construction fo the structure
    formula_units (int); This number is used to repeat the number of atoms in the chemical composition
    dimension (int); Dimension of the structure, right now only 3D
    volume_factor (float); Escalling factor for the unit cell
    
    uc_restr (Pyxtal's Lattice or bool); If the user provided cell restriction (a,b,c.alpha,beta,gamma),
        this is Pyxtal's Lattice. False, otherwise.

    out:
    xtalist_out (list); This list will contain all crystal structures as Molecule object
    '''
    vrestr = vol_restriction()
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
                    s_xtal = pyxtal2solids(xtal,dimension)
                    s_xtal = unit_cell_non_negative_coordinates(s_xtal)
                    s_xtal.c.append(0)
                    print('random_000_'+str(xc).zfill(3)+' ---> SG_'+str(sg_symbol)+"_("+str(sym)+")",file=fopen)
                    xtalist_out.append(s_xtal)
                    xc = xc + 1
            elif dimension == 3:
                try:
                    if uc_restr:
                        xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list, lattice=uc_restr)
                    else:
                        xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list)
                except:
                    continue
                else:
                    sg = Group (sym)
                    sg_symbol = str(sg.symbol)
                    s_xtal = pyxtal2solids(xtal,dimension)
                    s_xtal = unit_cell_non_negative_coordinates(s_xtal)
                    s_xtal.c.append(0)
                    print('random_000_'+str(xc).zfill(3)+' ---> SG_'+str(sg_symbol)+"_("+str(sym)+")",file=fopen)
                    xtalist_out.append(s_xtal)
                    xc = xc + 1
    fopen.close()
    xtalist_out = sort_by_stoichiometry(xtalist_out)
    return xtalist_out
#------------------------------------------------------------------------------------------------
def random_crystal_gen_GA(total_of_xtals,species,atoms_per_specie,p_list,formula_units=1,dimension=3,volume_factor=1.0,uc_restr=False):
    '''Generates a specified number of crystal structures

    in:
    total_of_xtals (int); The total number of crystals to be built
    species (list); This list contains each atomic symbol of the chemical composition 
    atoms_per_specie (list); This list contains the number of atoms corresponding to each species
    p_list (Tol_matrix); pyxtal.tolerance Tol_matrix object that contains the atomico overlap permited in 
        the construction fo the structure
    formula_units (int); This number is used to repeat the number of atoms in the chemical composition
    dimension (int); Dimension of the structure, right now only 3D
    volume_factor (float); Escalling factor for the unit cell
    
    uc_restr (float); There's two kinds of volume restriction: Using the lattice vectors or the value of
        the volume. If any of those is provided the structure's volume will be reecaled to it.

    out:
    xtalist_out (list); This list will contain all crystal structures as Molecule object
    '''
    xtalist_out = []
    xc = 1
    fopen = open(log_file,'a')
    sym_list = get_symmetry_constrains('symmetries', dimension)
    while xc <= total_of_xtals:
        xtal = pyxtal()
        sym = random.choice(sym_list)
        if dimension == 2:
            try:
                xtal.from_random(dimension,sym,species,atoms_per_specie,thickness=0.0)
            except:
                continue
            else:
                sg = Group (sym)
                sg_symbol = str(sg.symbol)
                s_xtal = pyxtal2solids(xtal,dimension)
                s_xtal = unit_cell_non_negative_coordinates(s_xtal)
                s_xtal.c.append(0)
                print('random_000_'+str(xc).zfill(3)+' ---> SG_'+str(sg_symbol)+"_("+str(sym)+")",file=fopen)
                xtalist_out.append(s_xtal)
                xc = xc + 1
        elif dimension == 3:
            try:
                if uc_restr:
                    xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list, lattice=uc_restr)
                else:
                    xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list)
            except:
                continue
            else:
                sg = Group (sym)
                sg_symbol = str(sg.symbol)
                s_xtal = pyxtal2solids(xtal,dimension)
                s_xtal = unit_cell_non_negative_coordinates(s_xtal)
                # if vol_restr:
                #     s_xtal = rescale_str(s_xtal,vol_restr)
                s_xtal.c.append(0)
                print('random_000_'+str(xc).zfill(3)+' ---> SG_'+str(sg_symbol)+"_("+str(sym)+")",file=fopen)
                xtalist_out.append(s_xtal)
                xc = xc + 1
    fopen.close()
    xtalist_out = sort_by_stoichiometry(xtalist_out)
    return xtalist_out
#------------------------------------------------------------------------------------------------
def popgen_fresh_random_gen(number_of_random,species,atoms_per_specie,volume_factor,p_list,dimension,index,uc_restr=False):
    '''Generates a specified number of random crystal structures

    in:
    number_of_random (int); The total number of crystals to be built
    species (list); This list contains each atomic symbol of the chemical composition 
    atoms_per_specie (list); This list contains the number of atoms corresponding to each species
    dimension (int); Dimension of the structure, right now only 3D
    
    uc_restr

    out:
    xtalist_out (list); This list will contain all crystal structures as Molecule object
    '''
    if number_of_random == 0:
        poscarout=[]
        return poscarout
    logfile = open(log_file,'a')
    print("-------------------------------------------------------------------", file=logfile)
    print("------------------------ RANDOM STRUCTURES ------------------------", file=logfile)
    xtalist_out = []
    xc = 1
    sym_list = get_symmetry_constrains('symmetries', dimension)
    basename = 'random_' + str(index).zfill(3) + '_'
    while xc <= number_of_random:
        sym = random.choice(sym_list)
        xtal = pyxtal()
        if dimension == 2:
            try:
                xtal.from_random(dimension,sym,species,atoms_per_specie,thickness=0.0)
            except:
                continue
            else:
                sg = Group (sym)
                sg_symbol = str(sg.symbol)
                s_xtal = pyxtal2solids(xtal,dimension)
                s_xtal = unit_cell_non_negative_coordinates(s_xtal)
                s_xtal.c.append(0)
                print(basename + str(xc).zfill(3) +'---> SG_'+str(sg_symbol)+"_("+str(sym)+")",file=logfile)
                xtalist_out.append(s_xtal)
                xc = xc + 1
        elif dimension == 3:
            try:
                if uc_restr:
                    xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list, lattice=uc_restr)
                else:
                    xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list)
            except:
                continue
            else:
                sg = Group (sym)
                sg_symbol = str(sg.symbol)
                s_xtal = pyxtal2solids(xtal,dimension)
                s_xtal = unit_cell_non_negative_coordinates(s_xtal)
                s_xtal.c.append(0)
                print(basename + str(xc).zfill(3) + ' ---> SG_'+str(sg_symbol)+"_("+str(sym)+")",file=logfile)
                xtalist_out.append(s_xtal) 
                xc = xc + 1
    print("We have %d POSCAR type RANDOM from %d solicited" %(len(xtalist_out), number_of_random), file=logfile)
    logfile.close()
    xtalist_out = sort_by_stoichiometry(xtalist_out)
    xtalist_out  = rename_molecule(xtalist_out, basename, 3)
    return xtalist_out