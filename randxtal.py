import random
from utils.libmoleculas import sort_by_stoichiometry
from inout.getbilparam import get_a_str
from pyxtal import pyxtal

from miscellaneous import pyxtal2xyz, possible_sym, rescale_str

log_file = get_a_str('output_file','glomos_out.txt')
#------------------------------------------------------------------------------------------------
def random_crystal_gen(total_of_xtals,species,atoms_per_specie,formula_units=1,dimension=3,volume_factor=1.1,vol_restr=False):
    '''
    This function generates the initial population for the algorithm using Pyxtal to build random structures
    and later translate them into GLOMOS Molecule object. It takes into consideration the restrictions imposed 
    by the user, if requested.

    in: total_of_xtals (int), the requested number of random structures
        species (list), the species of all the atoms in the unit cell ['Mg','O']
        atoms_per_specie (list), it concatenates with 'species' according with composition [1,2]
        formula_units (int), the number of formula units within the unit cell
        dimension (int), so far only 3D
        volumen_factor (float), a fixed value to expand or contract the unit cell
        vol_restr (Lattice), Pyxtal's Lattice object if requested, False otherwise

    out: xtal_list (list), This list contains all the generated structures
    '''
    xtal_list = []
    possym = possible_sym(atoms_per_specie)
    z = len(atoms_per_specie)
    # Expand the chemical formula to match the amount of formula units
    for i in range(z):
        atoms_per_specie[i] = atoms_per_specie[i] * formula_units
    c = 0
    for ii in range(total_of_xtals):
        # select one of the possible symmetries
        if c < len(possym):
            random_sym = possym[c]
        else:
            random.shuffle(possym)
            random_sym = random.choice(possym)
        c = c + 1
        # build the crystal with volume restrictions if requested
        xtal = pyxtal()
        xtal.from_random(dimension, random_sym, species, atoms_per_specie, volume_factor)
        # translate the crystal to GLOMOS and add them to xtal_list
        glomos_xtal = pyxtal2xyz(xtal)
        if vol_restr:
            glomos_xtal = rescale_str(glomos_xtal,vol_restr)
        glomos_xtal.i = 'random' + str(ii+1).zfill(3) + '_sym_' + random_sym.zfill(3)
        glomos_xtal.c.append(0)
        fopen = open(log_file,'a')
        print('random'+str(ii+1).zfill(3)+'_sym_'+random_sym.zfill(3),file=fopen)
        fopen.close()
        xtal_list.append(glomos_xtal)
    xtal_list = sort_by_stoichiometry(xtal_list)
    return xtal_list

#------------------------------------------------------------------------------------------------

def run_sample():
    from inout.getbilparam import get_a_int, get_a_float
    from inout.readbil import read_var_composition
    from vasp.libperiodicos import writeposcars
    from miscellaneous import get_xcomp
    composition = read_var_composition('composition')
    atms_specie, atms_per_specie = get_xcomp(composition)
    total_structures = get_a_int('total_structures', 4)
    formula_units = get_a_int('formula_units',1)
    dimension = get_a_int('dimension',3)
    volume_factor = get_a_float('volume_factor', 1.1)
    poscarlist=random_crystal_gen(total_structures,atms_specie,atms_per_specie, formula_units,dimension,volume_factor)
    writeposcars(poscarlist,'test.vasp','D')
# run_sample()
