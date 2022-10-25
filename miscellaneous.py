import numpy as np
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
def uc_restriction():
    '''
    This functions obtains the restrictions imposed by the user to the unit cell

    out: restr_uc (Lattice), the lattice object from Pyxtal if found in the INPUT.txt file
         or False if not found.
         ### we're testing wether if just reescaling the structures works better than creating 
         them from scratch with the lattice parameters
    '''
    import os.path
    from pyxtal.lattice import Lattice
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

# from vasp.libperiodicos import readposcars,writeposcars
# o = readposcars('ours.vasp')[0]
# # vol org = 123.32073085200001
# n = reescale_str(o,122)
# writeposcars([n],'compression.vasp','D')