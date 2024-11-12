import spglib
import numpy as np
from vasp_solids.libperiodicos import readposcars, writeposcars
from vasp_solids.libperiodicos import direct2cartesian, cartesian2direct
from utils_solids.atomic import get_atomic_number
from utils_solids.libmoleculas import sort_by_energy
from inout_solids.getbilparam import get_a_str, get_a_float
#----------------------------------------------------------------------------------------------------------
#Variables
log_file = get_a_str('output_file','solids_out.txt')
min_dE = get_a_float('min_energy_difference',0.01)
min_dV = get_a_float('min_volume_difference',0.01) 
#----------------------------------------------------------------------------------------------------------
def find_spacegroup(xtal_in):
	'''Uses spglib to find the space group of the crystal structures
	
	in: 
	xtal_in (Molecule); the structure of interest
	
	out: 
	spc_grp (str); the symbol for the space group of the structure
	'''
	lattice = xtal_in.m
	atms_position = []
	atms_specie = []
	for a in xtal_in.atoms:
		x,y,z = cartesian2direct(a.xc,a.yc,a.zc,xtal_in.m)
		x,y,z = round(x,2),round(y,2),round(z,2)
		atms_position.append([x,y,z])
		n = get_atomic_number(a.s)
		atms_specie.append(n)
	atms_specie.sort()
	x = (lattice,atms_position,atms_specie)
	spc_grp = spglib.get_spacegroup(x)
	return spc_grp

#----------------------------------------------------------------------------------------------------------
def get_volume(xtal_in):
	i,j,k = xtal_in.m[0],xtal_in.m[1],xtal_in.m[2]
	jxk = np.cross(j,k)
	vol = np.dot(i,jxk)
	return vol

#----------------------------------------------------------------------------------------------------------
def get_fingerprint(xtal_in):
	'''Builts the fingerprint of each structure. According to Zurek et.al. this comprises a list
	like this: [energy,spacial_group,volume]

	in:
	xtal_in (Molecule); The structures whose fingerprint is desired

	out:
	fp (list); [energy,spacial_group,volume]   
	'''
	e = xtal_in.e
	sg = find_spacegroup(xtal_in)
	v = get_volume(xtal_in)
	fp = [e,sg,v]
	return fp

#----------------------------------------------------------------------------------------------------------
def compare_fingerprints(xtal_a,xtal_b,vol_restr):
	'''Compares the fingerprints of two different structures, but takes into consideration if a volume restriction
	is in place. 

	in:
	xtal_a, xtal_b (Molecule); Two structures to be compared
	vol_restr; float, if there is a volume restriction; false, otherwise

	out:
	equal (bool); True if the structures have the same fingerprint, False otherwise
	edif (float); Energy difference between structures
	vdif (float); Volume difference between structures
	sym_a, sym_b (str); The symmetry found for each structure 
	'''
	t_atma, t_atmb = len(xtal_a.atoms), len(xtal_b.atoms)
	fp_a = get_fingerprint(xtal_a)
	fp_b = get_fingerprint(xtal_b)
	edif = abs(fp_a[0] - fp_b[0])/t_atma
	edif = round(edif,5)
	vdif = abs(fp_a[2] - fp_b[2])/t_atma
	vdif = round(vdif,5)
	sym_a, sym_b = fp_a[1],fp_b[1]
	equal = True
	if vol_restr:
		if edif <= min_dE and sym_a == sym_b:
			return equal, edif, vdif, sym_a, sym_b
		else:
			equal = False
			return equal, edif, vdif, sym_a, sym_b
	else:
		if edif <= min_dE and vdif <= min_dV:
			if sym_a == sym_b:
				return equal, edif, vdif, sym_a, sym_b
			elif edif == 0:
				return equal, edif, vdif, sym_a, sym_b
			else:
				equal = False
				return equal, edif, vdif, sym_a, sym_b
		else:
			equal = False
			return equal, edif, vdif, sym_a, sym_b
		
#---------------------------------------------------------------------------------------------------------
def discriminate_calculated(xtalist_in, vol_restr):
	'''Compares each structure with the remaining fo the list taking into consideration volume restriction

	in: 
	xtalist_in (list); List with all to-be-compared structures 
	vol_restr; Float if volume restriction, False otherwise

	out:
	xtalist_out (list); Curated list with only different structures
	'''
	fopen = open(log_file,'a')
	print('\n-------------------------------------------------------------------',file=fopen)
	print('--------------- Duplicates Removal in Generation ------------------',file=fopen)
	print('\nTolerance Values Given: Maximum dE = '+str(min_dE)+', Maximum dV = '+str(min_dV),file=fopen)
	size_in = len(xtalist_in)
	disc_count = 0
	xtalist_out = []
	for i,str_a in enumerate(xtalist_in):
		stop_flag = False
		for j in range(i+1,size_in):
			str_b = xtalist_in[j]
			e, dE, dV, sym_a, sym_b = compare_fingerprints(str_a,str_b,vol_restr)
			if e == True:
				print('%s PG:%18s removed, similar to %s PG:%18s, dE=%.5f dV=%.5f' %(str_a.i,sym_a,str_b.i,sym_b,dE,dV),file=fopen)
				disc_count = disc_count + 1
				stop_flag = True
				break			
		if stop_flag == False:
			xtalist_out.append(str_a)
	print('\n'+str(disc_count)+' structures removed by similarity in generation comparison',file=fopen)
	fopen.close()
	return xtalist_out

#---------------------------------------------------------------------------------------------------------
def discriminate_calculated_vs_pool(calulation_list, pool_list, vol_restr):
	'''Compares the list obtained in discriminate_calculated with the general pool of structures to prevent 
	repeating structures in the final result

	in: 
	calulation_list (list); The curated list with only different structures
	pool_list (list); List with every structure in the pool
	vol_restr; Float if volume restriction, False otherwise

	out:
	xtalist_out (list); Final list of curated structures, this will become the new pool
	'''
	fopen = open(log_file,'a')
	print('\n-------------------------------------------------------------------',file=fopen)
	print('---------------- Duplicates Removal Gen vs Pool -------------------\n',file=fopen)
	xtalist_out = []
	diff_elements = []
	disc_count = 0 
	for calc_str in calulation_list:
		stop_flag = False
		for pool_str in pool_list:
			e, dE, dV, sym_p, sym_c = compare_fingerprints(pool_str,calc_str,vol_restr)
			if e == True:
				print('%s PG: %18s removed, similar to %s PG: %18s, dE=%.5f dV=%.5f' %(calc_str.i,sym_c,pool_str.i,sym_p,dE,dV),file=fopen)
				disc_count = disc_count + 1
				stop_flag = True
				break
		if stop_flag == False:
			diff_elements.append(calc_str)
	if len(xtalist_out) == len(pool_list):
		print(str(disc_count)+' structures removed by similarity in Gen vs Pool comparison'+'\n',file=fopen)
	else:
		xtalist_out.extend(diff_elements)
		xtalist_out = sort_by_energy(xtalist_out,1)
		print('\n'+str(disc_count)+' structures removed by similarity in Gen vs Pool comparison''\n',file=fopen)
	fopen.close()
	return xtalist_out