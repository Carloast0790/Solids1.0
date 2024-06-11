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
#     xtal_list = []
#     z = len(atoms_per_specie)
#     for i in range(z):
#         atoms_per_specie[i] = atoms_per_specie[i] * formula_units
#     xc = 0
#     if dimension == 2:
#         topsym = 80
#     else:
#         topsym = 230
#     for sym in range (2,topsym):
#         if xc <= total_of_xtals:
#             if dimension == 2:
#                 xtal = pyxtal()
#                 xtal.from_random(dimension,sym,species,atoms_per_specie,thickness=0.0,force_pass=True)#, factor=0.8
#             else:
#                 xtal = pyxtal()
#                 xtal.from_random(dimension,sym,species,atoms_per_specie,volume_factor,p_list,force_pass=True)
#             if xtal.valid:
#                 glomos_xtal = pyxtal2xyz(xtal)
#                 glomos_xtal = unit_cell_non_negative_coordinates(glomos_xtal)
#                 if vol_restr:
#                     glomos_xtal = rescale_str(glomos_xtal,vol_restr)
#                 glomos_xtal.i = 'rand'+str(xc+1).zfill(3)+'_sym_'+str(sym).zfill(3)
#                 glomos_xtal.c.append(0)
#                 fopen = open(log_file,'a')
#                 print('rand'+str(xc+1).zfill(3)+'_sym_'+str(sym).zfill(3),file=fopen)
#                 fopen.close()
#                 xtal_list.append(glomos_xtal)
#                 xc = xc + 1
#         else:
#             break
#     xtal_list = sort_by_stoichiometry(xtal_list)
#     return xtal_list
