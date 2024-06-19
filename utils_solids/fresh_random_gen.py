import random
from utils_solids.libmoleculas import sort_by_stoichiometry
from inout_solids.getbilparam import get_a_str
from utils_solids.miscellaneous import pyxtal2xyz, rescale_str, unit_cell_non_negative_coordinates
import warnings
warnings.filterwarnings("ignore")
#------------------------------------------------------------------------------------------------
# Variables
log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------------
def popgen_fresh_random_gen(number_of_random,species,atoms_per_specie,generation,formula_units,dimension):
    '''Generates a specified number of random crystal structures

    in:
    number_of_random (int); The total number of crystals to be built
    species (list); This list contains each atomic symbol of the chemical composition 
    atoms_per_specie (list); This list contains the number of atoms corresponding to each species
    formula_units (int); This number is used to repeat the number of atoms in the chemical composition
    dimension (int); Dimension of the structure, right now only 3D

    out:
    xtal_list (list); This list will contain all crystal structures as Molecule object
    '''
    if number_of_random==0:
        poscarout=[]
        return poscarout
    from pyxtal import pyxtal
    logfile = open(log_file,'a')
    print("-------------------------------------------------------------------", file=logfile)
    print("------------------------ RANDOM STRUCTURES ------------------------", file=logfile)
    xtal_list = []
    xc = 0
    while xc <= number_of_random:
        xtal = pyxtal()
        if dimension == 2:
            try:
                sym = random.randint(2,80)
                xtal.from_random(dimension,sym,species,atoms_per_specie,thickness=0.0)
            except:
                continue
        elif dimension == 3:
            try:
                sym = random.randint(2,230)
                xtal.from_random(dimension, sym, species, atoms_per_specie)
            except:
                continue
            else:
                xc = xc + 1
                solids_xtal = pyxtal2xyz(xtal)
                solids_xtal = unit_cell_non_negative_coordinates(solids_xtal)
                name = 'random_'+str(generation).zfill(3)+'_'+str(xc).zfill(3)
                solids_xtal.i = name
                xtal_list.append(solids_xtal)
                print(name + ' ---> sym_' + str(sym), file=logfile)
        if xc == number_of_random:
            break
    print("We have %d POSCAR type RANDOM from %d solicited" %(len(xtal_list), number_of_random), file=logfile)
    logfile.close()
    xtal_list = sort_by_stoichiometry(xtal_list)
    return xtal_list