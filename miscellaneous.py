import numpy as np
#------------------------------------------------------------------------------------------------
#-------------------------First Random Generations-----------------------------------------------
#------------------------------------------------------------------------------------------------
def check_compatible(group, numIons):
    from copy import deepcopy
    """
    Checks if the number of atoms is compatible with the Wyckoff
    positions. Considers the number of degrees of freedom for each Wyckoff
    position, and makes sure at least one valid combination of WP's exists.
    NOTE Comprhys: Is degrees of freedom used symnomously with multiplicity?
    perhaps standardising to multiplicity would be clearer?
    """
    # Store whether or not at least one degree of freedom exists
    has_freedom = False
    # Store the wp's already used that don't have any freedom
    used_indices = []
    # Loop over species
    for numIon in numIons:
        # Get lists of multiplicity, maxn and freedom
        l_mult0 = []
        l_maxn0 = []
        l_free0 = []
        indices0 = []
        for i_wp, wp in enumerate(group):
            indices0.append(i_wp)
            l_mult0.append(len(wp))
            l_maxn0.append(numIon // len(wp))
            if np.allclose(wp[0].rotation_matrix, np.zeros([3, 3])):
                l_free0.append(False)
            else:
                l_free0.append(True)
        # Remove redundant multiplicities:
        l_mult = []
        l_maxn = []
        l_free = []
        indices = []
        for mult, maxn, free, i_wp in zip(l_mult0, l_maxn0, l_free0, indices0):
            if free is True:
                if mult not in l_mult:
                    l_mult.append(mult)
                    l_maxn.append(maxn)
                    l_free.append(True)
                    indices.append(i_wp)
            elif free is False and i_wp not in used_indices:
                l_mult.append(mult)
                indices.append(i_wp)
                if mult <= numIon:
                    l_maxn.append(1)
                elif mult > numIon:
                    l_maxn.append(0)
                    l_free.append(False)
        # Loop over possible combinations
        p = 0  # Create pointer variable to move through lists
        # Store the number of each WP, used across possible WP combinations
        n0 = [0] * len(l_mult)
        n = deepcopy(n0)
        for i, mult in enumerate(l_mult):
            if l_maxn[i] != 0:
                p = i
                n[i] = l_maxn[i]
                break
        p2 = p
        if n == n0:
            return False
        while True:
            num = np.dot(n, l_mult)
            dobackwards = False
            # The combination works: move to next species
            if num == numIon:
                # Check if at least one degree of freedom exists
                for val, free, i_wp in zip(n, l_free, indices):
                    if val > 0:
                        if free is True:
                            has_freedom = True
                        elif free is False:
                            used_indices.append(i_wp)
                break
            # All combinations failed: return False
            if n == n0 and p >= len(l_mult) - 1:
                return False
            # Too few atoms
            if num < numIon:
                # Forwards routine
                # Move p to the right and max out
                if p < len(l_mult) - 1:
                    p = p + 1
                    n[p] = min((numIon - num) // l_mult[p], l_maxn[p])
                elif p == len(l_mult) - 1:
                    # p is already at last position: trigger backwards routine
                    dobackwards = True
            # Too many atoms
            if num > numIon or dobackwards is True:
                # Backwards routine
                # Set n[p] to 0, move p backwards to non-zero, and decrease by 1
                n[p] = 0
                while p > 0 and p > p2:
                    p -= 1
                    if n[p] != 0:
                        n[p] -= 1
                        if n[p] == 0 and p == p2:
                            p2 = p + 1
                        break
    if has_freedom:
        # All species passed: return True
        return True
    else:
        # All species passed, but no degrees of freedom: return 0
        return 0

#------------------------------------------------------------------------------------------------
def possible_sym(total_atms):
    """
    This function returns a list of the possible symmetry groups for a given stoichiometry

    in: total_atms (int), the total amount of atoms in the stoicheometry
    out: pos_sym (list), a list containing the possible symmetry for the crystal
    """
    from pyxtal.symmetry import Group
    pos_sym = []
    for ii in range(1,230):
        sg = Group (ii)
        sg_symbol = str(sg.symbol)
        x = check_compatible(sg,total_atms)
        if x == True:
            pos_sym.append(sg_symbol)
    return pos_sym

#------------------------------------------------------------------------------------------------
def get_xcomp(composition):
    '''
    This function translates the composition that is returned from read_var_composition function in
    GLOMOS into a PyXtal readible format

    in: composition (list), the variable composition as read by GLOMOS
    out: species (list), a list with the species of the atoms in the crystal
         atms_per_specie (list), this list concatenates with species, matching the amount of atoms of each kind
    '''
    x = len(composition[0])
    species = []
    atms_per_specie = []
    for ii in range(x):
        s = composition[0][ii][0]
        species.append(s)
        n = composition[0][ii][1]
        atms_per_specie.append(n)
    return species, atms_per_specie

#------------------------------------------------------------------------------------------------
def pyxtal2xyz(xtal_from_pyxtal):
    """
    This function transforms the Pyxtal object to a Molecule object

    in: xtal_from_pyxtal (pyxtal), the object to be transformed
    out: glomos_xtal (Molecule), the transformed object
    """
    from utils.libmoleculas import Molecule, Atom
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
#---------------------------------General Constrains---------------------------------------------
#------------------------------------------------------------------------------------------------
def uc_restriction():
    '''
    This functions obtains the restrictions imposed by the user to the unit cell

    out: restr_uc (Lattice), the lattice object from Pyxtal if found in the INPUT.txt file
         or False if not found.
         ### we're testing wether if just reescaling the structures works better than creating 
         them from scratch with the lattice parameters
    '''
    import os.path
    #from pyxtal.lattice import Lattice
    file = 'INPUT.txt'
    if os.path.isfile(file):
        restr_v = False
        a,b,c = 0,0,0
        f=open(file,"r")
        for line in f:
            if '#' in line:
                continue
            if 'fixed_UC' in line:
                line = line.split()
                a,b,c = float(line[1]),float(line[2]),float(line[3])
                restr_v = a*b*c
                # restr_uc = Lattice.from_para(a,b,c,A,B,G)
                break
            elif 'fixed_vol' in line:
                line = line.split()
                restr_v = float(line[1])
                break
        f.close()
    return restr_v

#------------------------------------------------------------------------------------------------
def rescale_str(xtal_in, reference_volume):
    '''
    This function reescales the volume of a given cell to match a given parameter

    in: xtal_in (Molecule), the structure to be reescaled
        reference_volume (float), the desired volume
    out: xtal_out (Molecule), the reescaled structure 
    '''
    from utils.libmoleculas import Molecule, Atom
    from vasp.libperiodicos import cartesian2direct, direct2cartesian
    volume = lambda v0,v1,v2: abs(np.dot(np.cross(v0,v1),v2))
    mtx = xtal_in.m
    org_vol = volume(mtx[0],mtx[1],mtx[2])
    r = reference_volume / org_vol
    dv = np.cbrt([r])
    mtx[0],mtx[1],mtx[2] = dv*mtx[0],dv*mtx[1],dv*mtx[2]
    xtal_out = Molecule(xtal_in.i,xtal_in.e,mtx)
    for a in xtal_in.atoms:
        xc,yc,zc = cartesian2direct(a.xc,a.yc,a.zc,xtal_in.m)
        xn, yn,zn = direct2cartesian(xc,yc,zc,mtx)
        natm = Atom(a.s,xn,yn,zn)
        xtal_out.add_atom(natm)
    return xtal_out

#------------------------------------------------------------------------------------------------
def get_default_tolerances(species,scale_value):
    '''
    This function gets the default tolerances of each pair of atoms in the structure.
    
    in: species (list), a list of strings, with each species found in the crystal
        scale_value (float), a scaling value for the sum of each species' covalent radius.
        The scaled sum will be used as the minimum interatomic distance.
    out: tolerances (list), a list containing tuples with each min int dist, [(s1,s2,d1),(s1,s3,d2),...]
         py_tol (Tol_matrix), a PyXtal object used for the initial generation of structures
    '''
    from utils.libmoleculas import get_covalent_radius
    from pyxtal.tolerance import Tol_matrix
    py_tol = Tol_matrix()
    tolerances = []
    if len(species) == 1:
        s = species[0]
        r = get_covalent_radius(s)
        tv = r*scale_value*2
        tv = round(tv,2)
        tolerances.append((s,s,tv))
        py_tol.set_tol(s,s,tv)
        return tolerances,py_tol
    else:
        for i,s1 in enumerate(species):
            r1 = get_covalent_radius(s1)
            tv = r1*scale_value*2
            tv = round(tv,2)
            tolerances.append((s1,s1,tv))
            py_tol.set_tol(s1,s1,tv)
            for j in range(i+1,len(species)):
                s2 = species[j]
                r2 = get_covalent_radius(s2)
                tv_mix = (r1+r2)*scale_value
                tv_mix = round(tv_mix,2)
                tolerances.append((s1,s2,tv_mix))
                py_tol.set_tol(s1,s2,tv_mix)
        return tolerances,py_tol 

#------------------------------------------------------------------------------------------------
def get_custom_tolerances():
    '''
    This function obtains a list of the tuples containing each of the min interatomic distances 
    provided by the user.
    
    out: tolerances (list), a list containing tuples with each min int dist, [(s1,s2,d1),(s1,s3,d2),...]
         py_tol (Tol_matrix), a PyXtal object used for the initial generation of structures
    '''
    import os.path
    from pyxtal.tolerance import Tol_matrix
    file = 'INPUT.txt'
    py_tol = Tol_matrix()
    if os.path.isfile(file):
        bilfile = open(file,"r")
        for line in bilfile:
            if not line.startswith('#') and 'custom_tolerances' in line:
                tolerances = []
                line = line.lstrip('\n')
                readline = line.split()
                readline = readline[1:]
                for i in readline:
                    x = i.split(',')
                    tupla = (x[0],x[1],float(x[2]))
                    tolerances.append(tupla)
                    py_tol.set_tol(x[0],x[1],float(x[2]))
        bilfile.close()
    return tolerances,py_tol

#------------------------------------------------------------------------------------------------
def get_tolerances(species):
    from inout.getbilparam import get_a_float
    tv = get_a_float('interatom_scale_value',False)
    if tv:
        dtv,p_dtv = get_default_tolerances(species,tv)
        return dtv,p_dtv
    else:
        ctv, p_ctv = get_custom_tolerances()
        return ctv,p_ctv

#------------------------------------------------------------------------------------------------
# def virtual_expansion(xtal_in):
#     '''
#     This function takes the original structure and expands it in all dimensions in order to later check
#     the overlapping between atoms. The cell that will be of interest is the one in the middle of them  

#     in: xtal_in (Molecule), the structure that will be expanded
#     out: xtal_out (Molecule), the expanded structure 
#     '''
#     from utils.libmoleculas import copymol, Molecule, Atom
#     cp_xtal = copymol(xtal_in)
#     mtx = cp_xtal.m
#     a1, a2, a3 = mtx[0,:], mtx[1,:], mtx[2,:]
#     exp_mtx = np.array([3*a1,3*a2,3*a3])
#     xtal_out = Molecule(cp_xtal.i, cp_xtal.e, exp_mtx)
#     for ii,iatom in enumerate(cp_xtal.atoms):
#         for x in [-1,0,1]:
#             for y in [-1,0,1]:
#                 for z in [-1,0,1]:
#                     vt = float(x)*a1 + float(y)*a2 + float(z)*a3
#                     xc, yc, zc = iatom.xc + vt[0], iatom.yc + vt[1], iatom.zc + vt[2]
#                     if x==1 and y==1 and z==1:
#                         xf, yf, zf = iatom.xf, iatom.yf, iatom.zf
#                     else:
#                         xf, yf, zf = 'F','F','F'
#                     ai = Atom(iatom.s,xc,yc,zc,xf,yf,zf)
#                     xtal_out.add_atom(ai)
#     return xtal_out

def virtual_expansion(singleposcar, tol=0.01):
    import numpy as np
    from utils.libmoleculas import Atom, Molecule
    matrix=np.copy(singleposcar.m)
    a1, a2, a3=matrix[0,:], matrix[1,:], matrix[2,:]
    mi=np.linalg.inv(matrix)
    poscarout=Molecule(singleposcar.i, singleposcar.e, matrix)
    xdmin,xdmax=0.0,1.0
    ydmin,ydmax=0.0,1.0
    zdmin,zdmax=0.0,1.0
    for ii,iatom in enumerate(singleposcar.atoms):
        for x in [-1,0,1]:
            for y in [-1,0,1]:
                for z in [-1,0,1]:
                    vt=float(x)*a1+float(y)*a2+float(z)*a3
                    xc, yc, zc=iatom.xc+vt[0], iatom.yc+vt[1], iatom.zc+vt[2]
                    xf, yf, zf=iatom.xf, iatom.yf, iatom.zf
                    vector=np.array([xc, yc, zc])
                    vd=np.matmul(mi,vector)
                    xd, yd, zd = vd[0], vd[1], vd[2]
                    suma=0
                    if (xd >= xdmin-tol) and (xd <= xdmax+tol): suma=suma+1 
                    if (yd >= ydmin-tol) and (yd <= ydmax+tol): suma=suma+1 
                    #if (zd >= zdmin-tol) and (zd <= zdmax+tol): suma=suma+1 
                    if (zd >= zdmin) and (zd <= zdmax): suma=suma+1 
                    if suma == 3:
                        ai=Atom(iatom.s,xc,yc,zc,xf,yf,zf)
                        poscarout.add_atom(ai)
    return poscarout


#------------------------------------------------------------------------------------------------
def distance_atom_atom(iatom, jatom):
    vi=np.array([iatom.xc, iatom.yc, iatom.zc])
    vj=np.array([jatom.xc, jatom.yc, jatom.zc])
    c_dist = np.linalg.norm(vi-vj)
    return c_dist

#------------------------------------------------------------------------------------------------
def overlap_check(xatom, xtal_host,ref_dist):
    from utils.libmoleculas import copymol
    import copy
    # from vasp.libperiodicos import writeposcars
    new_xtal = copymol(xtal_host)
    cp_atom = copy.copy(xatom)
    for a in new_xtal.atoms:
        a.xf, a.yf, a.zf='F','F','F'
    cp_atom.xf,cp_atom.yf,cp_atom.zf='T','T','T'
    new_xtal.add_atom(cp_atom)
    extended_xtal = virtual_expansion(new_xtal,0.3)
    # writeposcars([extended_xtal],'ext_check.vasp','D')
    natoms = extended_xtal.n
    flag = False
    for i, ia in enumerate(extended_xtal.atoms):
        if ia.xf == 'T':
            for j in range(i+1,natoms):
                ja = extended_xtal.atoms[j]
                auxt = [ia.s,ja.s]
                auxt.sort()
                for t in ref_dist:
                    refsym = [t[0],t[1]]
                    refsym.sort()
                    if refsym == auxt:
                        dmin = float(t[2])
                        break
                d = distance_atom_atom(ia,ja)
                if d < dmin:
                    # print('traslape detectado entre',ia.s,ja.s)
                    flag = True
                    break
                # else:
                #     print('todo bien')
            if flag == True:
                break
    return flag


# def overlap_check(xatom, xtal_host,ref_dist):
#     from utils.libmoleculas import copymol
#     # from vasp.libperiodicos import writeposcars
#     new_xtal = copymol(xtal_host)
#     cp_atom = xatom
#     for a in new_xtal.atoms:
#         a.xf, a.yf, a.zf='F','F','F'
#     cp_atom.xf,cp_atom.yf,cp_atom.zf='T','T','T'
#     new_xtal.add_atom(cp_atom)
#     extended_xtal = virtual_expansion(new_xtal)
#     # writeposcars([extended_xtal],'ext_check.vasp','D')
#     # poscarhlp = unit_cell_non_negative_coordinates(poscarhlp)
#     natoms = extended_xtal.n
#     flag = False
#     for i, ia in enumerate(extended_xtal.atoms):
#         if ia.xf == 'T':
#             for j in range(i+1,natoms):
#                 ja = extended_xtal.atoms[j]
#                 auxt = [ia.s,ja.s]
#                 auxt.sort()
#                 for t in ref_dist:
#                     if t[0] == auxt:
#                         dmin = t[1]
#                     break
#                 d = distance_atom_atom(ia,ja)
#                 if d < dmin:
#                     # print('traslape detectado entre',ia.s,ja.s)
#                     flag = True
#                     break
#                 # else:
#                 #     print('todo bien')
#             if flag == True:
#                 break
#     return flag

#----------------------------------------------------------------------------------------------------------
def unit_cell_non_negative_coordinates(xtal_in):
    '''
    This function transforms every negative coordinate in the crystal structure to the equivalent positive one
    in order to correctly calculate the distances between atoms whithin the UC.
    '''
    from vasp.libperiodicos import cartesian2direct, direct2cartesian
    from utils.libmoleculas import copymol
    xtal_out = copymol(xtal_in)
    for a in xtal_out.atoms:
        x,y,z = cartesian2direct(a.xc,a.yc,a.zc,xtal_out.m)
        l = [x,y,z]
        for i in range(3):
            if l[i] < 0:
                l[i] = l[i] + 1
            elif l[i] > 1:
                l[i] = l[i] - 1
        a.xc, a.yc, a.zc = l[0], l[1], l[2]
        a.xc,a.yc,a.zc = direct2cartesian(a.xc,a.yc,a.zc,xtal_out.m)
    return xtal_out


# from vasp.libperiodicos import readposcars,writeposcars
# o = readposcars('test.vasp')[0]
# o = overlap_check(o)
# writeposcars([o],'no_neg.vasp','D')