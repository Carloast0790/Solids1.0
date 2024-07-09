import numpy as np
#------------------------------------------------------------------------------------------------
#----------------------------------First Random Generations--------------------------------------
#------------------------------------------------------------------------------------------------
def get_xcomp(composition):
    '''Translates the composition that is returned from read_var_composition function into a 
    PyXtal readible format

    in: 
    composition (list); the variable composition as read by read_var_composition
    
    out: 
    species (list); list with the species of the atoms in the crystal
    atms_per_specie (list); list concatenating with species, matching the amount of atoms of each kind
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
    """Transforms Pyxtal object to Molecule object

    in: 
    xtal_from_pyxtal (pyxtal); the object to be transformed

    out: 
    xtal (Molecule), Molecule object
    """
    from utils_solids.libmoleculas import Molecule, Atom
    coordinates, species = xtal_from_pyxtal._get_coords_and_species(True)
    lattice_vectors = xtal_from_pyxtal.lattice.get_matrix()
    xtal = Molecule('pyxtal2molecule', 0.0, lattice_vectors)
    total_atms = len(species)
    for ii in range(total_atms):
        si = species[ii]
        xc = coordinates[ii][0]
        yc = coordinates[ii][1]
        zc = coordinates[ii][2]
        iatm = Atom(si,xc,yc,zc)
        xtal.add_atom(iatm)
    return xtal

def pyxtal2solids(xtal_from_pyxtal,dimension):
    from utils_solids.libmoleculas import Molecule, Atom
    coordinates, species = xtal_from_pyxtal._get_coords_and_species(True)
    lattice_vectors = xtal_from_pyxtal.lattice.get_matrix()
    if dimension == 2:
        lattice_vectors[2] = lattice_vectors[2] * 6
    xtal = Molecule('pyxtal2molecule', 0.0, lattice_vectors)
    total_atms = len(species)
    zavg = sum([coordinates[i][2] for i in range(total_atms)])/total_atms
    for ii in range(total_atms):
        si = species[ii]
        xc = coordinates[ii][0]
        yc = coordinates[ii][1]
        if dimension == 2:
            zc = coordinates[ii][2] - zavg + lattice_vectors[2][2]/2
        else:
            zc = coordinates[ii][2]
        iatm = Atom(si,xc,yc,zc)
        xtal.add_atom(iatm)
    return xtal

#------------------------------------------------------------------------------------------------
#-------------------------------------General Constrains-----------------------------------------
#------------------------------------------------------------------------------------------------
def uc_restriction():
    '''Obtains the volume for the restrictions imposed to the unit cell

    out: 
    restr_uc (float); Desired volume (it is under evaluation if Pyxtal tools are necessary for this)
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
            if 'fixed_uc' in line:
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
    '''Reescales the volume of the unit cell to match a given parameter

    in: 
    xtal_in (Molecule); Structure to be reescaled
    reference_volume (float); Desired volume
    
    out: 
    xtal_out (Molecule); Reescaled structure 
    '''
    from utils_solids.libmoleculas import Molecule, Atom
    from vasp_solids.libperiodicos import cartesian2direct, direct2cartesian
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
    '''Gets the default tolerances for each pair of atoms in the structure
    
    in: 
    species (list); List with each species found in the structure
    scale_value (float); Scaling value for the sum of each species' covalent radius
    The scaled sum will be used as the minimum interatomic distance.
    
    out: 
    tolerances (list), List containing tuples with each min int dist, [(s1,s2,d1),(s1,s3,d2),...]
    py_tol (Tol_matrix), PyXtal object used for atomic tolerance in generation of structures
    '''
    from utils_solids.libmoleculas import get_covalent_radius
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
            tv = r1*scale_value
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
    '''Obtains a list of tuples containing each min interatomic distance provided by the user
    
    out: 
    tolerances (list); List containing tuples with each min int dist [(s1,s2,d1),(s1,s3,d2),...]
    py_tol (Tol_matrix); PyXtal object used for atomic tolerance in generation of structures
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
    '''Builds the tolerance pyxtal object based on the pair of species of each atom. It also returns a 
    list of tuples

    in:
    species (str); Atomic species

    out:
    dtv, ctv (list); list of tuples, dtv if default, ctv otherwise
    p_dtv, p_ctv (Tol_matrix); Pyxtal tolerance object, p_dtv if default, p_ctv otherwise
    '''
    from inout_solids.getbilparam import get_a_float
    tv = get_a_float('interatom_scale_value',False)
    if tv:
        dtv,p_dtv = get_default_tolerances(species,tv)
        return dtv,p_dtv
    else:
        ctv, p_ctv = get_custom_tolerances()
        return ctv,p_ctv

#------------------------------------------------------------------------------------------------
#--------------------------------------General Tools---------------------------------------------
#------------------------------------------------------------------------------------------------
def virtual_expansion(singleposcar, tol=0.01):
    '''Expands to a certain degree the size of the unit cell in order to search atomic overlaps

    in:
    singleposcar (Molecule); Structure at hand
    tol (float); Degree of expansion of the unit cell in all directions

    out:
    poscarout (Molecule); Expanded structure
    '''
    import numpy as np
    from utils_solids.libmoleculas import Atom, Molecule
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
    '''Calculates and returns c_dist, the distance between iatom and jatom
    '''
    vi=np.array([iatom.xc, iatom.yc, iatom.zc])
    vj=np.array([jatom.xc, jatom.yc, jatom.zc])
    c_dist = np.linalg.norm(vi-vj)
    return c_dist

#------------------------------------------------------------------------------------------------
def overlap_check(xatom, xtal_host,ref_dist):
    '''Receives a structure, expands it to a certain degree and checks for atomic overlapp

    in:
    xatom (Atom); To-be-checked and added atom, its .xf,.yf,.zf are marked as 'T' to distinguish
    xtal_host (Molecule); Structure that will receive the atom
    ref_dist (list); List containing tuples with minimum interatomic distances

    out:
    flag (bool); If the interatomic distance between xatom and other in the structure, flag=True,
        flag=False any other case   
    '''
    from utils_solids.libmoleculas import copymol
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
                    flag = True
                    break
            if flag == True:
                break
    return flag

#----------------------------------------------------------------------------------------------------------
def unit_cell_non_negative_coordinates(xtal_in):
    '''To correctly calculate the distances between atoms, it transforms every negative atomic coordinate 
        into its positive-valued equivalent

    in:
    xtal_in (Molecule); To-be-corrected structure

    out:
    xtal_out (Molecule); Corrected structure
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