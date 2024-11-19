import random
import numpy as np
from utils_solids.libmoleculas import Molecule, Atom, copymol, rename_molecule
from vasp_solids.libperiodicos import cartesian2direct, direct2cartesian
from inout_solids.getbilparam import get_a_int, get_a_str
from gega_solids.solids_roulette import get_roulette_wheel_selection
from utils_solids.miscellaneous import rescale_str

#number_of_mutants = get_a_int('number_of_mutants',5)
number_of_atmxchange = get_a_int('number_of_xchange',5)
number_of_lattstr = get_a_int('number_of_strains',5)
log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
def lattice_correction(xtal_in):
    xtal_out = copymol(xtal_in)
    vs = [xtal_out.m[0],xtal_out.m[1],xtal_out.m[2]]
    nvs = []
    for i,vi in enumerate(vs):
        for j in range(i+1,3):
            vj = vs[j]
            mvi,mvj = np.linalg.norm(vi), np.linalg.norm(vj)
            vdot = np.dot(vi,vj)
            vdot = round(vdot,6)
            # print('vdot',vdot)
            div = vdot/(mvi*mvj)
            mag = abs(np.arccos(div) - np.pi)
            if mag > np.pi/3 and mvi >= mvj:
                # print(vi,'needs change its mag',mag,'and',mvi,'greater than',mvj)
                div2 = abs(vdot/(mvj*mvj))
                # print('div',div2)
                sign = np.sign(vdot)
                # print('sign',sign)
                vi = vi - np.ceil(div2)*sign*vj
                # print('new vector',nv)
            nvs.append(vi)
            # print(nvs)
    xtal_out.m[0],xtal_out.m[1],xtal_out.m[2] = nvs[0],nvs[1],nvs[2]
    return xtal_out

#------------------------------------------------------------------------------------------
def strained_lattice_unrestricted(original_lattice):
    '''
    This functions multiplies a random strain matrix by each of the lattice vectors of a crystal unit cell

    in: original_lattice (numpy array)
    out: mutated_lattice (numpy array) 
    '''
    aux = []
    cont = 0 
    while cont <= 5:
        x = random.gauss(0,1)
        if abs(x) <= 0.5:
            aux.append(x)
            cont = cont + 1        
    # random.shuffle(aux)
    e11, e22, e33 = 1 + aux[0], 1 + aux[1], 1 + aux[2]
    e12, e13, e23 = aux[3] / 2, aux[4] / 2, aux[5] / 2 
    strain_matrix = np.array([[e11,e12,e13],[e12,e22,e23],[e13,e23,e33]])
    mutated_lattice = np.array([np.dot(original_lattice[ii],strain_matrix) for ii in range(3)])
    return mutated_lattice

#------------------------------------------------------------------------------------------
def lattice_mutation(xtal_in):
    '''
    This function mutates the crystal lattice by appling a stain matrix which entries are taken 
    randomly from a gaussian distribution. The matrix is then applied to all lattice vectors 
    of the cell and the cell is reescaled to have the original volume.

    in: xtal_in (Molecule), the structure which UC will be mutated
    out: xtal_out (Molecule), the structure with the mutated UC
    '''
    
    o_lat = xtal_in.m
    volume = lambda v0,v1,v2: abs(np.dot(np.cross(v0,v1),v2))
    org_vol = volume(o_lat[0],o_lat[1],o_lat[2])
    str_lat = strained_lattice_unrestricted(o_lat)
    name = '_lattice_strain'
    xtal_out = Molecule(xtal_in.i+name,xtal_in.e,str_lat)
    for a in xtal_in.atoms:
        ox,oy,oz = cartesian2direct(a.xc,a.yc,a.zc,o_lat)
        nxc,nyc,nzc = direct2cartesian(ox,oy,oz,str_lat)
        n_atm = Atom(a.s,nxc,nyc,nzc)
        xtal_out.add_atom(n_atm)
    xtal_out = rescale_str(xtal_out,org_vol)
    # xtal_out = lattice_correction(xtal_out)
    return xtal_out

#------------------------------------------------------------------------------------------
def atom_exchange(xtal_in,rounds):
    xtal_out = copymol(xtal_in)
    l = len(xtal_out.atoms)
    for _ in range(rounds):
        atm1 = random.choice(xtal_out.atoms)
        s1,x1,y1,z1 = atm1.s,atm1.xc,atm1.yc,atm1.zc
        flag = 0
        while flag == 0:
            atm2 = random.choice(xtal_out.atoms)
            s2,x2,y2,z2 = atm2.s,atm2.xc,atm2.yc,atm2.zc
            if s2 != s1:
                atm1.xc = x2
                atm1.yc = y2
                atm1.zc = z2
                atm2.xc = x1
                atm2.yc = y1
                atm2.zc = z1
                flag = 1
                break
    xtal_out.i = xtal_in.i + '_atom_exchange'
    return xtal_out

#--------------------------------------------------------------------------------------------
def make_mutants(the_chosen_ones_list):
    xtal_out = []
    xchange_list = random.choices(the_chosen_ones_list, k = number_of_atmxchange)
    strain_list = random.choices(the_chosen_ones_list, k = number_of_lattstr)
    s_list = [a.s for a in the_chosen_ones_list[0].atoms]
    s_list = list(dict.fromkeys(s_list))
    if len(s_list) == 1:
        for s_str in strain_list:
            mut_s = lattice_mutation(s_str)
            xtal_out.append(mut_s)
    else:
        for s_str in strain_list:
            mut_s = lattice_mutation(s_str)
            xtal_out.append(mut_s)        
        for x_str in xchange_list:
            r = random.randint(1,4)
            mut_x = atom_exchange(x_str,r)
            xtal_out.append(mut_x)
    return xtal_out

#--------------------------------------------------------------------------------------------
def popgen_mutants(xtal_list, generation):
    if number_of_lattstr == 0 and number_of_atmxchange == 0:
        xtal_out = []
        return xtal_out
    number_of_mutants = number_of_lattstr + number_of_atmxchange
    logfile = open(log_file,'a')
    print ("-------------------------------------------------------------------", file=logfile)
    print ("-----------------------GENERATOR OF MUTANTS------------------------", file=logfile)
    logfile.close()
    the_chosen_ones = get_roulette_wheel_selection(xtal_list, number_of_mutants)
    xtal_out = make_mutants(the_chosen_ones)
    # xtal_out = xtal_out[0:number_of_mutants]
    logfile = open(log_file,'a')
    cont = 1
    for x in xtal_out:
        aux_name = 'mutant_' + str(generation).zfill(3) + '_' + str(cont).zfill(3)
        print ('%s ---> %s' %(aux_name,x.i), file=logfile)
        cont = cont + 1
    print ("We have %d POSCAR type MUTANT from %d solicited" %(len(xtal_out), number_of_mutants), file=logfile)
    logfile.close()
    name = 'mutant_' + str(generation).zfill(3)
    xtal_out = rename_molecule(xtal_out, name +'_', 3)
    return xtal_out

# from vasp_solids.libperiodicos import readposcars, writeposcars 

# x = readposcars('same1.vasp')[0]
# y = lattice_mutation(x)
# writeposcars([y],'mutated_lattice.vasp','D')
# # z = lattice_correction(y)
# # writeposcars([z],'corrected.vasp','D')


# def make_mutants(the_chosen_ones_list):
#     xtal_out = []
#     s_list = []
#     [s_list.append(a.s) for a in the_chosen_ones_list[0].atoms if a.s not in s_list]
#     l = len(s_list)
#     for xtal in the_chosen_ones_list:
#         tmp_xtal = copymol(xtal)
#         mut_type = random.gauss(0,1)
#         if mut_type < 0 and l > 1:
#             r = random.randint(1,4)
#             muty = atom_exchange(xtal,r)
#         else: 
#             muty = lattice_mutation(xtal)
#         xtal_out.append(muty)
#     return xtal_out

# def strained_lattice_unrestricted(original_lattice):
#     # Nota ahorita es para pruebas, originalmente es el strained_lattice_restricted
#     '''
#     This functions multiplies a random strain matrix by all the lattice vectors of a crystal unit cell

#     in: original_lattice (numpy array)
#     out: mutated_lattice (numpy array) 
#     '''
#     flag = False
#     stop = 0
#     while flag == False:
#         aux = []
#         cont = 0 
#         while cont <= 5:
#             x = random.gauss(0,1)
#             if abs(x) <= 1:
#                 aux.append(x)
#                 cont = cont + 1        
#         random.shuffle(aux)
#         e11, e22, e33 = 1 + aux[0], 1 + aux[1], 1 + aux[2]
#         e12, e13, e23 = aux[3] / 2, aux[4] / 2, aux[5] / 2 
#         strain_matrix = np.array([[e11,e12,e13],[e12,e22,e23],[e13,e23,e33]])
#         mutated_lattice = np.array([np.dot(original_lattice[ii],strain_matrix) for ii in range(3)])
#         for i in range(len(mutated_lattice)):
#             vi = mutated_lattice[i]
#             mi = np.linalg.norm(vi)
#             for j in range(i+1,len(mutated_lattice)):
#                 vj = mutated_lattice[j]
#                 mj = np.linalg.norm(vj)
#                 div = np.dot(vi,vj)/(mi*mj)
#                 theta = np.arccos(div)
#                 theta = round(theta,2)
#                 if theta > 1.05 and theta < 2.10:
#                     flag = True
#                 else:
#                     flag = False
#                     break
#             if flag == False:
#                 break
#     return mutated_lattice

# def lattice_mutation(xtal_in):
#     '''
#     This function mutates the crystal lattice by appling a stain matrix which entries are taken 
#     randomly from a gaussian distribution. The matrix is then applied to all lattice vectors 
#     of the cell and the cell is reescaled to have the original volume.

#     in: xtal_in (Molecule), the structure which UC will be mutated
#     out: xtal_out (Molecule), the structure with the mutated UC
#     '''
    
#     o_lat = xtal_in.m
#     str_lat = strained_lattice_unrestricted(o_lat)
#     name = '_lattice_strain'
#     xtal_out = Molecule(xtal_in.i+name,xtal_in.e,str_lat)
#     for a in xtal_in.atoms:
#         ox,oy,oz = cartesian2direct(a.xc,a.yc,a.zc,o_lat)
#         nxc,nyc,nzc = direct2cartesian(ox,oy,oz,str_lat)
#         n_atm = Atom(a.s,nxc,nyc,nzc)
#         xtal_out.add_atom(n_atm)
#     # xtal_out = lattice_correction(xtal_out)
#     return xtal_out