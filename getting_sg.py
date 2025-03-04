#!/home/carlos0790/miniconda3/bin/python -u
from pymatgen.io.vasp import Poscar
import spglib

poscar = Poscar.from_file("c8_graphite.vasp")
structure = poscar.structure
lattice = structure.lattice.matrix # 3x3 array of lattice vectors
positions = structure.frac_coords # Fractional coordinates of atoms
atomic_numbers = structure.atomic_numbers # List of atomic numbers
cell = (lattice, positions, atomic_numbers)
spacegroup = spglib.get_spacegroup(cell)
print("Space group:", spacegroup)
