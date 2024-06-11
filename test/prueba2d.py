from pyxtal import pyxtal
from vasp_solids.libperiodicos import readposcars, writeposcars
import numpy as np
from utils_solids.libmoleculas import get_covalent_radius
import random

def pyxtal2xyz_2D(xtal_from_pyxtal, vacuum=6):
    """Transforms Pyxtal 2D object to Molecule object adding vacuum 

    in: 
    xtal_from_pyxtal (pyxtal); the object to be transformed

    out: 
    xtal (Molecule), Molecule object
    """
    from utils_solids.libmoleculas import Molecule, Atom
    from utils_solids.libmoleculas import get_covalent_radius
    coordinates, species = xtal_from_pyxtal._get_coords_and_species(True)
    lattice_vectors = xtal_from_pyxtal.lattice.get_matrix()

    all_species = list(dict.fromkeys(species))
    print(all_species)
    vol_sum = 0
    for s in all_species:
        c = species.count(s)
        r = get_covalent_radius(s)
        atmvol = (np.pi*pow(r,2))
        print('radio',r,'conteo',c,'volumen de atomo',atmvol)
        vol_sum = vol_sum + atmvol

    area = np.cross(lattice_vectors[0],lattice_vectors[1])
    area = np.linalg.norm(area)
    atoms = len(species)
    div= atoms/area
    print('area',area,'div',div,'\n')
    red_vol = 0.8

    if div <= 0.1:
        lattice_vectors[0]=lattice_vectors[0]*red_vol
        lattice_vectors[1]=lattice_vectors[1]*red_vol
    lattice_vectors[2] = lattice_vectors[2] * vacuum
    xtal = Molecule('pyxtal2molecule', 0.0, lattice_vectors)
    for ii in range(atoms):
        si = species[ii]
        if div <= 0.1:
            xc = coordinates[ii][0]*red_vol
            yc = coordinates[ii][1]*red_vol
        else:
            xc = coordinates[ii][0]
            yc = coordinates[ii][1]
        zc = coordinates[ii][2] * vacuum
        iatm = Atom(si,xc,yc,zc)
        xtal.add_atom(iatm)
    return xtal

# x =readposcars('initial.vasp')[0]
# l,m=get_min_binding_distance_xy(x)
xtals = []
for x in range(3):
    sym = random.randint(2,80)
    print(sym)
    xtal = pyxtal()
    xtal.from_random(2,sym,['C'],[6],thickness=0.0,force_pass=True)
    if xtal.valid:
        ase_struc = xtal.to_ase()
        ase_struc.write('ase'+str(x)+f'.vasp', format='vasp', vasp5=True, direct=True)
        glomos_xtal = pyxtal2xyz_2D(xtal)
        xtals.append(glomos_xtal)
        writeposcars(xtals,'strs'+str(x)+'.vasp','D')
writeposcars(xtals,'initial.vasp','D')