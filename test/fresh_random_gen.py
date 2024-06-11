import random
from pyxtal import pyxtal
from utils_solids.libmoleculas import sort_by_stoichiometry
from inout_solids.getbilparam import get_a_str, get_a_int
from utils_solids.miscellaneous import pyxtal2xyz, rescale_str, unit_cell_non_negative_coordinates
#------------------------------------------------------------------------------------------------
# Variables
log_file = get_a_str('output_file','solids_out.txt')
num_of_random = get_a_int('num_of_random',10)
#------------------------------------------------------------------------------------------------
def fresh_random_gen(total_of_xtals,species,atoms_per_specie,formula_units=1,dimension=3):
    '''Generates a specified number of random crystal structures

    in:
    total_of_xtals (int); The total number of crystals to be built
    species (list); This list contains each atomic symbol of the chemical composition 
    atoms_per_specie (list); This list contains the number of atoms corresponding to each species
    formula_units (int); This number is used to repeat the number of atoms in the chemical composition
    dimension (int); Dimension of the structure, right now only 3D

    out:
    xtal_list (list); This list will contain all crystal structures as Molecule object
    '''
    xtal_list = []
    z = len(atoms_per_specie)
    for i in range(z):
        atoms_per_specie[i] = atoms_per_specie[i] * formula_units
    xc = 0
    if dimension == 2:
        topsym = 80
    else:
        topsym = 230
    while xc <= total_of_xtals:
        if dimension == 2:
            sym = random.randint(2,80)
            xtal = pyxtal()
            xtal.from_random(dimension,sym,species,atoms_per_specie,thickness=0.0, force_pass=True)
        elif dimension == 3:
            sym = random.randint(2,230)
            xtal = pyxtal()
            xtal.from_random(dimension, sym, species, atoms_per_specie,force_pass=True)        
        if xtal.valid:
            xc = xc + 1
            glomos_xtal = pyxtal2xyz(xtal)
            glomos_xtal = unit_cell_non_negative_coordinates(glomos_xtal)
            glomos_xtal.i = 'fresh_rand'+str(xc).zfill(3)+'_sym_'+str(sym).zfill(3)
            glomos_xtal.c.append(0)
            xtal_list.append(glomos_xtal)
            xc = xc + 1
        if xc == total_of_xtals:
            break
    xtal_list = sort_by_stoichiometry(xtal_list)
    return xtal_list
