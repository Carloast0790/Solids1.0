from ase.data import covalent_radii, atomic_numbers
import numpy as np
import py3Dmol
#------------------------------------------------------------------------------------------
def viewmol_RDKit(atoms, width=300, height=300):
    view = py3Dmol.view(width=width, height=height)
    xyz_str=atoms2xyz_str(atoms)
    view.addModel(xyz_str, 'xyz')
    view.setStyle({'sphere': {'scale': 0.3}, 'stick': {'radius': 0.2}})
    view.setBackgroundColor("white")
    view.zoomTo()
    view.show()
#------------------------------------------------------------------------------------------
def viewmolecule(atoms, width=300, height=300):
    view = py3Dmol.view(width=width, height=height)
    xyz_str = f"{len(atoms)}\n"
    xyz_str += f"\n"
    for atom in atoms:
        symbol = atom.symbol
        x, y, z = atom.position
        xyz_str += f"{symbol:2s} {x:14.9f} {y:16.9f} {z:16.9f}\n"
    view.addModel(xyz_str, 'xyz')
    view.setStyle({'sphere': {'scale': 0.3}, 'stick': {'radius': 0.2}})
    view.setBackgroundColor("white")
    view.zoomTo()
    view.show()
#------------------------------------------------------------------------------------------
def celda(a1, a2, a3, view):
    vertices = [np.array([0, 0, 0]), a1, a2, a3, a1 + a2, a1 + a3, a2 + a3, a1 + a2 + a3]
    vertices_dict = [{'x': float(v[0]), 'y': float(v[1]), 'z': float(v[2])} for v in vertices]
    aristas = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]
    for arista in aristas:
        start = vertices_dict[arista[0]]
        end = vertices_dict[arista[1]]
        view.addCylinder({
            'start': start,'end': end,'radius': 0.05,'color': 'black','fromCap': True,'toCap': True})
#------------------------------------------------------------------------------------------
def viewcrystal(atoms, width=300, height=300):
    view = py3Dmol.view(width=width, height=height)
    xyz_str = f"{len(atoms)}\n"
    xyz_str += f"\n"
    for atom in atoms:
        symbol = atom.symbol
        x, y, z = atom.position
        #view.addLabel(f"{atom.symbol}",{'position':{'x':x,'y':y,'z':z}})
        xyz_str += f"{symbol:2s} {x:14.9f} {y:16.9f} {z:16.9f}\n"
    view.addModel(xyz_str, 'xyz')
    for atom in atoms:
      radio = covalent_radii[atomic_numbers[atom.symbol]]
      view.setStyle(
           { "elem": atom.symbol},
           {'sphere': {'scale': radio*0.3}, 'stick': {'radius': 0.2}})
    a1,a2,a3=atoms.cell
    celda(a1, a2, a3 , view)
    view.setBackgroundColor("white")
    view.zoomTo()
    view.show()
#------------------------------------------------------------------------------------------
def viewmol_ASE(atoms, width=300, height=300):
    if any(atoms.pbc):
        viewmolecule(atoms, width, height)
    else:
        viewcrystal(atoms, width, height)
#------------------------------------------------------------------------------------------
