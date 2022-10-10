import random
import numpy as np
from utils.libmoleculas import Molecule, Atom, copymol, rename_molecule, get_min_binding_distance
from vasp.libperiodicos import cartesian2direct, direct2cartesian
from inout.getbilparam import get_a_int, get_a_str
from solids_roulette import get_roulette_wheel_selection

number_of_mutants = get_a_int('number_of_mutants',5)
log_file = get_a_str('output_file','glomos_out.txt')
#------------------------------------------------------------------------------------------
def strained_lattice_restricted(original_lattice):
    '''
    This functions multiplies a random strain matrix by all the lattice vectors of a crystal cell

    in: original_lattice (numpy array)
    out: mutated_lattice (numpy array) 
    '''
    flag = False
    stop = 0
    while flag == False:
        aux = []
        cont = 0 
        while cont <= 5:
            x = random.gauss(0,0.2)
            if abs(x) <= 1:
                aux.append(x)
                cont = cont + 1        
        random.shuffle(aux)
        e11, e22, e33 = 1 + aux[0], 1 + aux[1], 1 + aux[2]
        e12, e13, e23 = aux[3] / 2, aux[4] / 2, aux[5] / 2 
        strain_matrix = np.array([[e11,e12,e13],[e12,e22,e23],[e13,e23,e33]])
        mutated_lattice = np.array([np.dot(original_lattice[ii],strain_matrix) for ii in range(3)])
        for i in range(len(mutated_lattice)):
            vi = mutated_lattice[i]
            mi = np.linalg.norm(vi)
            for j in range(i+1,len(mutated_lattice)):
                vj = mutated_lattice[j]
                mj = np.linalg.norm(vj)
                div = np.dot(vi,vj)/(mi*mj)
                theta = np.arccos(div)
                theta = round(theta,2)
                if theta > 1.05 and theta < 2.10:
                    flag = True
                else:
                    flag = False
                    break
            if flag == False:
                break
    return mutated_lattice

#------------------------------------------------------------------------------------------
def strained_lattice_unrestricted(original_lattice):
    '''
    This functions multiplies a random strain matrix by all the lattice vectors of a crystal cell

    in: original_lattice (numpy array)
    out: mutated_lattice (numpy array) 
    '''
    aux = []
    cont = 0 
    while cont <= 5:
        x = random.gauss(0,0.2)
        if abs(x) <= 1:
            aux.append(x)
            cont = cont + 1        
    random.shuffle(aux)
    e11, e22, e33 = 1 + aux[0], 1 + aux[1], 1 + aux[2]
    e12, e13, e23 = aux[3] / 2, aux[4] / 2, aux[5] / 2 
    strain_matrix = np.array([[e11,e12,e13],[e12,e22,e23],[e13,e23,e33]])
    mutated_lattice = np.array([np.dot(original_lattice[ii],strain_matrix) for ii in range(3)])
    return mutated_lattice
#------------------------------------------------------------------------------------------
# def lattice_mutation(xtal_in):
#     """
#     This function mutates the crystal lattice by appling a stain matrix which entries
#     are taken randomly from a gaussian distribution. The matrix is applied to all lattice vectors 
#     of the cell and the cell is reescaled to have the original volume.

#     in: xtal_in (Molecule), the structure whose UC will be mutated
#     out: xtal_out (Molecule), the structure with the mutated UC
#     """ 
#     o_lat = xtal_in.m
#     volume = lambda v0,v1,v2: abs(np.dot(np.cross(v0,v1),v2))  
#     org_vol = volume(o_lat[0],o_lat[1],o_lat[2])
#     # ref_vol_b = org_vol*0.9
#     # ref_vol_t = org_vol*1.1
#     # new_vol = 0.0
#     # while new_vol < ref_vol_b or new_vol > ref_vol_t  :
#     #     str_lat = strained_lattice(o_lat)
#     #     new_vol = volume(str_lat[0],str_lat[1],str_lat[2])
#     str_lat = strained_lattice(o_lat)
#     str_volume = volume(str_lat[0],str_lat[1],str_lat[2])
#     if str_volume < org_vol:
#         str_lat[0],str_lat[1],str_lat[2] = str_lat[0]*1.1,str_lat[1]*1.1,str_lat[2]*1.1 
#     elif str_volume > org_vol:
#         str_lat[0],str_lat[1],str_lat[2] = str_lat[0]*0.9,str_lat[1]*0.9,str_lat[2]*0.9 
#     xtal_out = Molecule(xtal_in.i,xtal_in.e,str_lat)
#     for a in xtal_in.atoms:
#         ox,oy,oz = cartesian2direct(a.xc,a.yc,a.zc,o_lat)
#         nxc,nyc,nzc = direct2cartesian(ox,oy,oz,str_lat)
#         n_atm = Atom(a.s,nxc,nyc,nzc)
#         xtal_out.add_atom(n_atm)
#     return xtal_out

def lattice_mutation(xtal_in):
    """
    This function mutates the crystal lattice by appling a stain matrix which entries
    are taken randomly from a gaussian distribution. The matrix is applied to all lattice vectors 
    of the cell and the cell is reescaled to have the original volume.

    in: xtal_in (Molecule), the structure whose UC will be mutated
    out: xtal_out (Molecule), the structure with the mutated UC
    """ 
    o_lat = xtal_in.m
    str_lat = strained_lattice_unrestricted(o_lat)
    name = '_lattice_strain'
    xtal_out = Molecule(xtal_in.i+name,xtal_in.e,str_lat)
    for a in xtal_in.atoms:
        ox,oy,oz = cartesian2direct(a.xc,a.yc,a.zc,o_lat)
        nxc,nyc,nzc = direct2cartesian(ox,oy,oz,str_lat)
        n_atm = Atom(a.s,nxc,nyc,nzc)
        xtal_out.add_atom(n_atm)
    return xtal_out

#------------------------------------------------------------------------------------------
def atom_exchange(xtal_in,rounds):
    xtal_out = copymol(xtal_in)
    xtal_out.m = xtal_out.m * 1.1
    l = len(xtal_out.atoms)
    for _ in range(rounds):
        atm1 = random.choice(xtal_out.atoms)
        s1,x1,y1,z1 = atm1.s,atm1.xc,atm1.yc,atm1.zc
        cont = 0
        while cont <= l:
          atm2 = random.choice(xtal_out.atoms)
          s2,x2,y2,z2 = atm2.s,atm2.xc,atm2.yc,atm2.zc
          if s2 != s1:
             atm1.xc = x2
             atm1.yc = y2
             atm1.zc = z2
             atm2.xc = x1
             atm2.yc = y1
             atm2.zc = z1
             break
          cont = cont + 1
    return xtal_out

#----------------------------------------------------------------------------------------------------------
def make_mutants(the_chosen_ones_list):
    xtal_out = []
    for xtal in the_chosen_ones_list:
        tmp_xtal = copymol(xtal)
        mut_type = random.gauss(0,1)
        if mut_type > 0: 
            muty = lattice_mutation(xtal)
        if mut_type < 0: 
            muty = atom_exchange(xtal,3)
            muty.i = xtal.i + '_atom_exchange'
        xtal_out.append(muty)
    return xtal_out

#----------------------------------------------------------------------------------------------------------
def popgen_mutants(xtal_list, generation):
    if number_of_mutants==0:
        xtal_out=[]
        return xtal_out
    logfile = open(log_file,'a')
    print ("-------------------------------------------------------------------", file=logfile)
    print ("-----------------------GENERATOR OF MUTANTS------------------------", file=logfile)
    logfile.close()
    the_chosen_ones = get_roulette_wheel_selection(xtal_list, number_of_mutants)
    xtal_out = make_mutants(the_chosen_ones)
    xtal_out = xtal_out[0:number_of_mutants]
    logfile = open(log_file,'a')
    cont = 1
    for x in xtal_out:
        aux_name = 'mutant_' + str(generation).zfill(3) + '_' + str(cont).zfill(3)
        print ('%s ---> %s' %(aux_name,x.i), file=logfile)
        cont = cont + 1
    print ("We have %d POSCAR type MUTANT from %d solicited" %(len(xtal_out), number_of_mutants), file=logfile)
    logfile.close()
    name = 'mutant_' + str(generation).zfill(3)
    xtal_out = rename_molecule(xtal_out, name +'_', 4)
    return xtal_out

###### Nota mental, agregar desplazamientos aleatorios de los átomos dentro de la celda unitaria y extensión de la misma en una sola dirección