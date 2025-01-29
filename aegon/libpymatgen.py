from pymatgen.core import Molecule
#from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
#------------------------------------------------------------------------------------------
def molase2pymatgen(moleculease):
    symbols = moleculease.get_chemical_symbols()
    positions = moleculease.get_positions()
    moleculepymatgen=Molecule(symbols, positions)
    return moleculepymatgen
#------------------------------------------------------------------------------------------
def inequivalent_finder(moleculease, tolerance=0.3, eigen_tolerance=0.008, matrix_tolerance=0.1):
    molpmg=molase2pymatgen(moleculease)
    #molpmg=AseAtomsAdaptor.get_structure(moleculease)
    molx=PointGroupAnalyzer(molpmg, tolerance, eigen_tolerance, matrix_tolerance)
    dictionary=molx.get_equivalent_atoms()
    dict=dictionary['eq_sets']
    inequiv=[ix for ix in dict.keys()]
    #for ix in dict.keys(): print(list(dict[ix]))
    #print('')
    inequiv.sort()
    return inequiv
#------------------------------------------------------------------------------------------
