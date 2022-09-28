import random
from utils.libmoleculas import Molecule, Atom, sort_by_stoichiometry
from inout.getbilparam import get_a_str
from pyxtal import pyxtal
from ws_check import possible_sym

log_file=get_a_str('output_file','glomos_out.txt')
#------------------------------------------------------------------------------------------------
def pyxtal2xyz(xtal_from_pyxtal):
    """
    This function transforms the Pyxtal object to a Molecule object

    in: xtal_from_pyxtal (pyxtal), the object to be transformed
    out: glomos_xtal (Molecule), the transformed object
    """
    energy = float(0.0)
    name = 'pyxtal2glomos'
    # Getting the coordinates, species and lattice vectors
    coordinates, species = xtal_from_pyxtal._get_coords_and_species(True)
    lattice_vectors = xtal_from_pyxtal.lattice.get_matrix()
    # Creating the object Molecule in GLOMOS
    glomos_xtal = Molecule(name, energy, lattice_vectors)
    total_atms = len(species)
    for ii in range(total_atms):
        si = species[ii]
        xc = coordinates[ii][0]
        yc = coordinates[ii][1]
        zc = coordinates[ii][2]
        iatm = Atom(si,xc,yc,zc)
        glomos_xtal.add_atom(iatm)
    return glomos_xtal

#------------------------------------------------------------------------------------------------
def random_crystal_gen(total_of_xtals,species,atoms_per_specie,formula_units=1,dimension=3,volume_factor=1.1):
    '''
    This function generates the initial population for the algorithm using Pyxtal to build the unit cell 
    for random structures.

    in: total_of_xtals (int), the requested amount of random structures
        species (list), the species of all the atoms in the unit cell ['Mg','O',...]
        atoms_per_specie (list), it concatenates with 'species' according with the chemical formula [1,2,...]
        formula_units (int), the amount of desired formula units in the unit cell
        dimension (int), so far only 3D
        volumen_factor (float), a fixed value to expand or contract the unit cell

    out: xtal_list (list), This list contains all the generated Molecule objects
    '''
    xtal_list = []
    possym = possible_sym(atoms_per_specie)
    random.shuffle(possym)
    z = len(atoms_per_specie)
    # Expand the chemical formula to match the amount of formula units
    for i in range(z):
        atoms_per_specie[i] = atoms_per_specie[i] * formula_units
    # Find the possible symmetries
    c = 0
    for ii in range(total_of_xtals):
        # select one of the possible symmetries and generate the random crystal
        if c < len(possym):
            random_sym = possym[c]
        else:
            random_sym = random.choice(possym)
        c = c + 1
        xtal = pyxtal()
        xtal.from_random(dimension, random_sym, species, atoms_per_specie, volume_factor)
        # translate the crystal to GLOMOS and add them to xtal_list
        glomos_xtal = pyxtal2xyz(xtal)
        glomos_xtal.i = 'random' + str(ii+1).zfill(3) + '_sym_' + random_sym.zfill(3)
        glomos_xtal.c.append(0)
        fopen = open(log_file,'a')
        print('random'+str(ii+1).zfill(3)+'_sym_'+random_sym.zfill(3),file=fopen)
        fopen.close()
        xtal_list.append(glomos_xtal)
    xtal_list = sort_by_stoichiometry(xtal_list)
    return xtal_list

#------------------------------------------------------------------------------------------------
def get_xcomp(composition):
    x = len(composition[0])
    species = []
    atms_per_specie = []
    for ii in range(x):
        s = composition[0][ii][0]
        species.append(s)
        n = composition[0][ii][1]
        atms_per_specie.append(n)
        return species, atms_per_specie

def run_sample():
    from inout.getbilparam import get_a_int, get_a_float
    from inout.readbil import read_var_composition
    from vasp.libperiodicos import writeposcars
    composition = read_var_composition('composition')
    atms_specie, atms_per_specie = get_xcomp(composition)
    total_structures = get_a_int('total_structures', 4)
    formula_units = get_a_int('formula_units',1)
    dimension = get_a_int('dimension',3)
    volume_factor = get_a_float('volume_factor', 1.1)
    poscarlist=random_crystal_gen(total_structures,atms_specie,atms_per_specie, formula_units,dimension,volume_factor)
    writeposcars(poscarlist,'initial.vasp','D')
# run_sample()
