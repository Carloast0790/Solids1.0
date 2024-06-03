from pyxtal import pyxtal
from pyxtal.tolerance import Tol_matrix
from vasp_solids.libperiodicos import readposcars, writeposcars
# from utils_solids.miscellaneous import pyxtal2xyz

def pyxtal2xyz_2D(xtal_from_pyxtal, vacuum=6):
    """Transforms Pyxtal 2D object to Molecule object

    in: 
    xtal_from_pyxtal (pyxtal); the object to be transformed

    out: 
    xtal (Molecule), Molecule object
    """
    from utils_solids.libmoleculas import Molecule, Atom
    coordinates, species = xtal_from_pyxtal._get_coords_and_species(True)
    print(coordinates,species)
    lattice_vectors = xtal_from_pyxtal.lattice.get_matrix()
    print(lattice_vectors)
    print(lattice_vectors[2])
    lattice_vectors[2] = lattice_vectors[2] * vacuum
    print(lattice_vectors)
    xtal = Molecule('pyxtal2molecule', 0.0, lattice_vectors)
    total_atms = len(species)
    for ii in range(total_atms):
        si = species[ii]
        xc = coordinates[ii][0]
        yc = coordinates[ii][1]
        zc = coordinates[ii][2] * vacuum
        iatm = Atom(si,xc,yc,zc)
        xtal.add_atom(iatm)
    return xtal


# tol_m_1 = Tol_matrix(prototype="molecular", factor=1.0)
for x in range(10):
	xtal = pyxtal()
	# xtal.from_random(2,25,['C'],[6],2,thickness=0.0,tm=tol_m_1)
	xtal.from_random(dim=2,group=25,species=['C'],numIons=[8],thickness=0.0)
	ase_struc = xtal.to_ase()
	ase_struc.write('ase_'+str(x)+'.vasp',format='vasp',vasp5=True, direct=True)
	# solids
	glomos_xtal = pyxtal2xyz_2D(xtal)
	writeposcars([glomos_xtal],'solids_'+str(x)+'.vasp','D')
	break