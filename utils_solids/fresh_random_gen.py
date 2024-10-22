import random
from pyxtal import pyxtal
from utils_solids.libmoleculas import sort_by_stoichiometry,rename_molecule
from inout_solids.getbilparam import get_a_str
from utils_solids.miscellaneous import pyxtal2solids, rescale_str, unit_cell_non_negative_coordinates
import warnings
warnings.filterwarnings("ignore")
#------------------------------------------------------------------------------------------------
# Variables
log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------------
def popgen_fresh_random_gen(number_of_random,species,atoms_per_specie,volume_factor,p_list,dimension,index):
    '''Generates a specified number of random crystal structures

    in:
    number_of_random (int); The total number of crystals to be built
    species (list); This list contains each atomic symbol of the chemical composition 
    atoms_per_specie (list); This list contains the number of atoms corresponding to each species
    dimension (int); Dimension of the structure, right now only 3D

    out:
    xtalist_out (list); This list will contain all crystal structures as Molecule object
    '''
    if number_of_random==0:
        poscarout=[]
        return poscarout
    logfile = open(log_file,'a')
    print("-------------------------------------------------------------------", file=logfile)
    print("------------------------ RANDOM STRUCTURES ------------------------", file=logfile)
    xtalist_out = []
    xc = 0
    basename = 'random_' + str(index+1).zfill(3) + '_'
    while xc <= number_of_random:
        xtal = pyxtal()
        if dimension == 2:
            try:
                sym = random.randint(2,80)
                xtal.from_random(dimension,sym,species,atoms_per_specie,thickness=0.0)
            except:
                continue
            else:
                xc = xc + 1
                s_xtal = pyxtal2solids(xtal,dimension)
                s_xtal = unit_cell_non_negative_coordinates(s_xtal)
                s_xtal.c.append(0)
                print(basename + str(xc).zfill(3) +'---> sym_' + str(sym).zfill(3),file=logfile)
                xtalist_out.append(s_xtal)
        elif dimension == 3:
            try:
                sym = random.randint(2,230)
                xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list)
            except:
                continue
            else:
                xc = xc + 1
                s_xtal = pyxtal2solids(xtal,dimension)
                s_xtal = unit_cell_non_negative_coordinates(s_xtal)
                s_xtal.c.append(0)
                print(basename + str(xc).zfill(3) + ' ---> sym_' + str(sym).zfill(3),file=logfile)
                xtalist_out.append(s_xtal)        
        if xc == number_of_random:
            break
    print("We have %d POSCAR type RANDOM from %d solicited" %(len(xtalist_out), number_of_random), file=logfile)
    logfile.close()
    xtalist_out = sort_by_stoichiometry(xtalist_out)
    xtalist_out  = rename_molecule(xtalist_out, basename, 3)
    return xtalist_out

# popgen_fresh_random_gen(5,['Ti','O'],[1,2],3)
