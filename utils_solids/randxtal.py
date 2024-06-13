from pyxtal import pyxtal
from utils_solids.libmoleculas import sort_by_stoichiometry
from inout_solids.getbilparam import get_a_str, get_a_int
from utils_solids.miscellaneous import pyxtal2xyz, rescale_str, unit_cell_non_negative_coordinates

# import warnings
# warnings.filterwarnings("ignore")
#------------------------------------------------------------------------------------------------
# Variables
log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------------
# def random_crystal_gen(total_of_xtals,species,atoms_per_specie,p_list,formula_units=1,dimension=3,volume_factor=1.1,vol_restr=False):
#     '''Generates a specified number of crystal structures

#     in:
#     total_of_xtals (int); The total number of crystals to be built
#     species (list); This list contains each atomic symbol of the chemical composition 
#     atoms_per_specie (list); This list contains the number of atoms corresponding to each species
#     p_list (Tol_matrix); pyxtal.tolerance Tol_matrix object that contains the atomic overlap permited in 
#         the construction fo the structure
#     formula_units (int); This number is used to repeat the number of atoms in the chemical composition
#     dimension (int); Dimension of the structure, right now only 3D
#     volume_factor (float); Escalling factor for the unit cell
#     vol_restr (float); There's two kinds of volume restriction: Using the lattice vectors or the value of
#         the volume. If any of those is provided the structure's volume will be reecaled to it.

#     out:
#     xtal_list (list); This list will contain all crystal structures as Molecule object
#     '''
#     import random
#     xtal_list = []
#     xc = 0
#     if dimension == 2:
#         topsym = 80
#     else:
#         topsym = 230
#     while xc <= total_of_xtals:
#         if dimension == 2:
#             sym = random.randint(2,80)
#             xtal = pyxtal()
#             xtal.from_random(dimension,sym,species,atoms_per_specie,thickness=0.0, force_pass=True)
#         elif dimension == 3:
#             sym = random.randint(2,230)
#             xtal = pyxtal()
#             xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list,force_pass=True)        
#         if xtal.valid:
#             xc = xc + 1
#             glomos_xtal = pyxtal2xyz(xtal)
#             glomos_xtal = unit_cell_non_negative_coordinates(glomos_xtal)
#             if vol_restr:
#                 glomos_xtal = rescale_str(glomos_xtal,vol_restr)
#             glomos_xtal.i = 'rand'+str(xc).zfill(3)+'_sym_'+str(sym).zfill(3)
#             glomos_xtal.c.append(0)
#             fopen = open(log_file,'a')
#             print('rand'+str(xc).zfill(3)+'_sym_'+str(sym).zfill(3),file=fopen)
#             fopen.close()
#             xtal_list.append(glomos_xtal)
#             xc = xc + 1
#         if xc == total_of_xtals:
#             break
#     xtal_list = sort_by_stoichiometry(xtal_list)
#     return xtal_list

def random_crystal_gen(total_of_xtals,species,atoms_per_specie,p_list,formula_units=1,dimension=3,volume_factor=1.1,vol_restr=False):
    '''Generates a specified number of crystal structures

    in:
    total_of_xtals (int); The total number of crystals to be built
    species (list); This list contains each atomic symbol of the chemical composition 
    atoms_per_specie (list); This list contains the number of atoms corresponding to each species
    p_list (Tol_matrix); pyxtal.tolerance Tol_matrix object that contains the atomic overlap permited in 
        the construction fo the structure
    formula_units (int); This number is used to repeat the number of atoms in the chemical composition
    dimension (int); Dimension of the structure, right now only 3D
    volume_factor (float); Escalling factor for the unit cell
    vol_restr (float); There's two kinds of volume restriction: Using the lattice vectors or the value of
        the volume. If any of those is provided the structure's volume will be reecaled to it.

    out:
    xtal_list (list); This list will contain all crystal structures as Molecule object
    '''
    import random
    xtal_list = []
    xc = 0
    fopen = open(log_file,'a')
    while xc <= total_of_xtals:
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
                xtal.from_random(dimension, sym, species, atoms_per_specie,volume_factor,p_list)
            except:
                continue
            else:
                xc = xc + 1
                glomos_xtal = pyxtal2xyz(xtal)
                glomos_xtal = unit_cell_non_negative_coordinates(glomos_xtal)
                if vol_restr:
                    glomos_xtal = rescale_str(glomos_xtal,vol_restr)
                glomos_xtal.i = 'rand'+str(xc).zfill(3)+'_sym_'+str(sym).zfill(3)
                glomos_xtal.c.append(0)
                print('rand'+str(xc).zfill(3)+'_sym_'+str(sym).zfill(3),file=fopen)
                xtal_list.append(glomos_xtal)
                xc = xc + 1
        if xc == total_of_xtals:
            break
    fopen.close()
    xtal_list = sort_by_stoichiometry(xtal_list)
    return xtal_list
#------------------------------------------------------------------------------------------------
def run_sample():
    # from inout.getbilparam import get_a_int, get_a_float
    # from inout.readbil import read_var_composition
    # from vasp.libperiodicos import writeposcars
    # from miscellaneous import get_xcomp,get_tolerances
    # composition = read_var_composition('composition')
    # atms_specie, atms_per_specie = get_xcomp(composition)
    # l_tol,p_tol = get_tolerances(atms_specie)
    # total_structures = get_a_int('total_structures', 4)
    # formula_units = get_a_int('formula_units',1)
    # dimension = get_a_int('dimension',3)
    # volume_factor = get_a_float('volume_factor', 1.1)
    from vasp.libperiodicos import writeposcars 
    from pyxtal.lattice import Lattice
    from miscellaneous import pyxtal2xyz
    #cell = Lattice.from_para(3,3,20,90,90,90)
    total_structures=10
    atms_specie='C'
    atms_per_specie=[4]
    dimension=2
    xtal = pyxtal()
    xtal.from_random(2,12,['C'],[4],thickness=0.0, factor=0.8)#lattice=cell)
    glomos_xtal = pyxtal2xyz(xtal)
    writeposcars([glomos_xtal],'test.vasp','D')
# run_sample()