import random
import numpy as np
from utils.atomic  import get_covalent_radius
from utils.libmoleculas import copymol, Molecule, Atom, sort_by_stoichiometry, molecular_stoichiometry, rename_molecule
from vasp.libperiodicos import direct2cartesian, cartesian2direct

from solids_roulette import get_roulette_wheel_selection
from inout.getbilparam    import get_a_int, get_a_str
#----------------------------------------------------------------------------------------------------------
number_of_childs = get_a_int('number_of_matings',10)
log_file = get_a_str('output_file','glomos_out.txt')
#----------------------------------------------------------------------------------------------------------
def bounded_random_point(r_max,r_min):
	'''
	This function returns a random number (random_point) such as r_min < random_point < r_max
	Note: 0 < r_max & r_min < 1 
	'''
	random_point = 1
	while random_point > r_max or random_point < r_min:
		if r_max > 1 or r_min < 0:
			random_point = 0.5
			break
		random_point = random.random()
	random_point = round(random_point,2)
	return random_point

#----------------------------------------------------------------------------------------------------------
def combined_unit_cell(xtal_a, xtal_b, weight_h=0.6, weight_s=0.4):
	'''
	This function averages the matrices that comprise the unit cells of the crystal structures. The weighted
	averaged unit cell is returned.

	in: xtal_a, xtal_b (Molecule), the two structures whose UC will be averaged
	    weight_h (float), the higher weight to be given to the vectors
	    weight_s (float), the lower weigth to be given to the vectors
	out: new_mtx (Numpy Array), the weighted averaged UC
	'''
	mtx_a,mtx_b = xtal_a.m, xtal_b.m
	new_mtx=[]
	for i in range(3):
		ai = mtx_a[i]
		bi = mtx_b[i]
		mag_ai = np.linalg.norm(ai)
		mag_bi = np.linalg.norm(bi)
		if mag_ai > mag_bi:
			ci = bi*weight_s + ai*weight_h
		elif mag_ai < mag_bi:
			ci = ai*weight_s + bi*weight_h
		else:
			ci = ai*0.5 + bi*0.5
		new_mtx.append(ci)
	return np.array(new_mtx)

#----------------------------------------------------------------------------------------------------------
def unit_cell_non_negative_coordinates(xtal_in):
	'''
	This function transforms every negative coordinate in the crystal structure to the equivalent positive one
	in order to correctly calculate the distances between atoms whithin the UC.
	'''
	xtal_out = copymol(xtal_in)
	for a in xtal_out.atoms:
		x,y,z = cartesian2direct(a.xc,a.yc,a.zc,xtal_out.m)
		l = [x,y,z]
		for i in range(3):
			if l[i] < 0:
				l[i] = l[i] + 1
		a.xc, a.yc, a.zc = l[0], l[1], l[2]
		a.xc,a.yc,a.zc = direct2cartesian(a.xc,a.yc,a.zc,xtal_out.m)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def parent_cut_bellow(xtal_in, axis, cutting_point):
	'''
	This function copies into xtal_out the atoms of xtal_in that are under cutting_point in the selected axis

	in: xtal_in (Molecule), the cristal to be cut
	    axis (str), a, b or c correspondig to the axis in which the cu will be performed
	    cutting_point (float), under this point (0,1), the atoms will be saved in xtal_out
	out: xtal_out (Molecule), the structure that contains the atoms under cuttin_point
	'''
	cp_xtal = copymol(xtal_in)
	name = cp_xtal.i + str('_cutB_') + str(cutting_point) + axis 
	xtal_out = Molecule(name, cp_xtal.e, cp_xtal.m)
	for a in cp_xtal.atoms:
		x,y,z = cartesian2direct(a.xc,a.yc,a.zc,cp_xtal.m)
		s = a.s
		if axis == 'a':
			if x < cutting_point:
				xc,yc,zc = direct2cartesian(x,y,z,cp_xtal.m)
				atom = Atom(s,xc,yc,zc)
				xtal_out.add_atom(atom)
		if axis == 'b':
			if y < cutting_point:
				xc,yc,zc=direct2cartesian(x,y,z,cp_xtal.m)
				atom = Atom(s,xc,yc,zc)
				xtal_out.add_atom(atom)
		if axis == 'c':
			if z < cutting_point:
				xc,yc,zc=direct2cartesian(x,y,z,cp_xtal.m)
				atom = Atom(s,xc,yc,zc)
				xtal_out.add_atom(atom)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def parent_cut_above(xtal_in, axis, cutting_point):
	'''
	This function is exactly the same than parent_cut_bellow but cutting above the cutting_point
	'''
	cp_xtal = copymol(xtal_in)
	name = cp_xtal.i + str('_cutA_') + str(cutting_point) + axis 
	new_xtal = Molecule(name, cp_xtal.e, cp_xtal.m)
	for a in cp_xtal.atoms:
		x,y,z = cartesian2direct(a.xc,a.yc,a.zc,cp_xtal.m)
		s = a.s
		if axis == 'a':
			if x > cutting_point:
				xc,yc,zc = direct2cartesian(x,y,z,cp_xtal.m)
				atom = Atom(s,xc,yc,zc)
				new_xtal.add_atom(atom)
		if axis == 'b':
			if y > cutting_point:
				xc,yc,zc = direct2cartesian(x,y,z,cp_xtal.m)
				atom = Atom(s,xc,yc,zc)
				new_xtal.add_atom(atom)
		if axis == 'c':
			if z > cutting_point:
				xc,yc,zc = direct2cartesian(x,y,z,cp_xtal.m)
				atom = Atom(s,xc,yc,zc)
				new_xtal.add_atom(atom)
	return new_xtal
#----------------------------------------------------------------------------------------------------------
def translating_to_avg_uc(xtal_in,avg_uc):
	'''
	This function translates xtal_in to the average UC, that is, the same atomic position but in the avg UC

	in: xtal_in (Molecule), the original crystal structure that will be translated
	    avg_uc (Numpy array), th UC upon wich the translation will be done.
	out: xtal_out (Molecule), the new translated crystal structure 
	'''
	cp_xtal = copymol(xtal_in)
	name = cp_xtal.i + str('_avg_uc') 
	xtal_out = Molecule(name,cp_xtal.e,avg_uc)
	for iatom in cp_xtal.atoms:
		x,y,z = cartesian2direct(iatom.xc,iatom.yc,iatom.zc,cp_xtal.m)
		s = iatom.s
		xn,yn,zn = direct2cartesian(x,y,z,avg_uc)
		new_atom = Atom(s,xn,yn,zn)
		xtal_out.add_atom(new_atom) 
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def missing_atoms_identifier(original_stoich,new_stoich):
	'''
	This function compares the original composition of the structure and gets the amount of 
	missing atoms and their species.

	in: original_stoich (list[tuples]), a list that has all tuples with (symbol,#atoms)
	    new_stoich (list[tuples]), the composition after the cut(symbol,#atoms)
	out: missing_atoms[list[tuples]], the amount of atoms missing on each symbol (symbol,#atoms) 
	'''
	full_stoich = original_stoich + new_stoich
	cont = 0
	counted_atoms = []
	missing_atoms = []
	for ii in range(len(original_stoich)):
		symbol = full_stoich[ii][0]
		for jtup in full_stoich:
			if symbol == jtup[0]:
				cont = cont + jtup[1]
				tup_aux = (symbol,cont)
		counted_atoms.append(tup_aux)	
		cont = 0
	for itup,jtup in zip(original_stoich,counted_atoms):
		x,y = itup[1], jtup[1]
		dif = y - 2*x
		if dif < 0:
			tup_aux = (itup[0],abs(dif))
			missing_atoms.append(tup_aux)
	return missing_atoms

#----------------------------------------------------------------------------------------------------------
def virtual_expansion(xtal_in):
	'''
	This function takes the original structure and expands it in all dimensions in order to later check
	the overlapping between atoms. The cell that will be of interest is the one in the middle of them  

	in: xtal_in (Molecule), the structure that will be expanded
	out: xtal_out (Molecule), the expanded structure 
	'''
	cp_xtal = copymol(xtal_in)
	mtx = cp_xtal.m
	a1, a2, a3 = mtx[0,:], mtx[1,:], mtx[2,:]
	exp_mtx = np.array([3*a1,3*a2,3*a3])
	xtal_out = Molecule(cp_xtal.i, cp_xtal.e, exp_mtx)
	for ii,iatom in enumerate(cp_xtal.atoms):
		for x in [-1,0,1]:
			for y in [-1,0,1]:
				for z in [-1,0,1]:
					vt = float(x)*a1 + float(y)*a2 + float(z)*a3
					xc, yc, zc = iatom.xc + vt[0], iatom.yc + vt[1], iatom.zc + vt[2]
					if x==1 and y==1 and z==1:
						xf, yf, zf = iatom.xf, iatom.yf, iatom.zf
					else:
						xf, yf, zf = 'F','F','F'
					ai = Atom(iatom.s,xc,yc,zc,xf,yf,zf)
					xtal_out.add_atom(ai)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def min_dist_atm_xtal(xatom, xtal_host):
	'''
	This function calculates the distance between the added atom, marked as T, and the rest of them,
	marked as F, in the same expanded structure

	in: xatom (Atom), the recently added atom
	    xtal_host (Molecule), the structure that will receive the atom and be expanded
	out: dmin (float), the minimal distance found between the added atom and everything else
	'''
	new_xtal = copymol(xtal_host)
	cp_atom = xatom
	for a in new_xtal.atoms:
		a.xf, a.yf, a.zf='F','F','F'
	cp_atom.xf,cp_atom.yf,cp_atom.zf='T','T','T'
	new_xtal.add_atom(cp_atom)
	poscarhlp = virtual_expansion(new_xtal)
	poscarhlp = unit_cell_non_negative_coordinates(poscarhlp)
	natoms = poscarhlp.n
	dmin = 100.0
	for iatom in range(natoms):
		qi = poscarhlp.atoms[iatom].xf
		if qi == 'T':
			si = poscarhlp.atoms[iatom].s
			ri = get_covalent_radius(si)
			xi = poscarhlp.atoms[iatom].xc
			yi = poscarhlp.atoms[iatom].yc
			zi = poscarhlp.atoms[iatom].zc
			vi = np.array([xi, yi, zi])
			for jatom in range(natoms):
				if jatom != iatom:
					qj = poscarhlp.atoms[jatom].xf
					sj = poscarhlp.atoms[jatom].s
					rj = get_covalent_radius(sj)
					xj = poscarhlp.atoms[jatom].xc
					yj = poscarhlp.atoms[jatom].yc
					zj = poscarhlp.atoms[jatom].zc
					vj = np.array([xj, yj, zj])
					dist = np.linalg.norm(vj-vi)
					dist = dist /  (ri + rj)
					if dist < dmin:
						dmin = dist
	return dmin

#----------------------------------------------------------------------------------------------------------
def crossover(base_xtal,complement_xtal):
	'''
	This function gets two structures, randomly transfers each unit cell in 3 dimensions, cuts them 
	in a random point located on a random vector, and unites the halves into one single structure

	in: base_xtal, complement_xtal (Molecule), the two structures that will be attempted to crossed
	out: xtal_out, (Molecule or False), the crossover between the structures if passible, False otherwise
	'''
	from translation_rotation import translation_3D,rotation2D
	# copy the original xtals and get their stoichiometry
	base = copymol(base_xtal)
	base = unit_cell_non_negative_coordinates(base)
	org_stoich = molecular_stoichiometry(base,0)
	comp = copymol(complement_xtal)
	comp = unit_cell_non_negative_coordinates(comp)
	# find the averaged unit cell
	avg_uc = combined_unit_cell(base, comp)
	# translate them into this new uc
	avg_base = translating_to_avg_uc(base,avg_uc)
	avg_comp = translating_to_avg_uc(comp,avg_uc)
	stop = 0
	wflag = False
	# Note: ref_dist was selected randomly and this value seems to work effectively
	ref_dist = 0.5
	while wflag == False:
		iflag = True
		# find the random vector and point on which the cut will be performed
		r_vect = random.choice(['a','b','c'])
		r_point = bounded_random_point(0.6,0.2)
		# cut the structures
		xtal_out = parent_cut_bellow(avg_base,r_vect,r_point)
		cut_comp = parent_cut_above(avg_comp,r_vect,r_point)
		# writeposcars([xtal_out],'cuta.vasp','D')
		# writeposcars([cut_comp],'cutb.vasp','D')
		xtal_out = translation_3D(xtal_out)
		cut_comp = translation_3D(cut_comp)
		# find the stoichiometry of the structures after the cut
		new_stoich = molecular_stoichiometry(xtal_out,0)
		comp_stoich = molecular_stoichiometry(cut_comp,0)
		# find how many atoms of each kind are missing
		m_atm = missing_atoms_identifier(org_stoich,new_stoich)
		# for each tuple in m_atm, get the symbol and amount of missing atoms
		for t in m_atm:
			# iflag tells if there are enough atoms in comp to perform the union
			if iflag == False:
				break
			else:
				# to know iflag, explore the available stoichiometry
				ms,ma = t[0],t[1]
				for tt in comp_stoich:
					avs,ava = tt[0], tt[1]
					# if the amount of available atoms is less than the required one break the cycle 
					if avs == ms and ava < ma:
						iflag = False
						break
		# the search for overlaps will occur only if the structures are complementary
		if iflag == True:
			# for each tuple in m_atm, get the symbol and amount of missing atoms
			for t in m_atm:
				ms, ma = t[0],t[1]
				# search for the symbol in the comp structure
				for a in cut_comp.atoms:
					s = a.s
					if s == ms:
						# once found, check for overlaps
						if ma == 0:
							continue
						try_atom = Atom(s,a.xc,a.yc,a.zc)
						min_dist = min_dist_atm_xtal(try_atom, xtal_out)
						# if there are non overlaps, add the atom to the structure
						if min_dist > ref_dist:
							xtal_out.add_atom(try_atom)
							ma = ma - 1
			# organize the resulting crossover
			xtal_out = sort_by_stoichiometry([xtal_out])[0]
			# writeposcars([xtal_out],'out.vasp','D')
		# get its stoichiometry and compare it with the original 
		out_stoich = molecular_stoichiometry(xtal_out,0)
		# if there was succes,break the cycle, if there weren't try again
		if out_stoich == org_stoich:
			wflag = True
			xtal_out.i = str(base_xtal.i) + '_x_' + str(complement_xtal.i)
		else:
			stop = stop + 1
		# if after n attempts there were no luck, break the cycle and return False
		if stop == 10:
			xtal_out = False
			break
	return xtal_out

#----------------------------------------------------------------------------------------------------------
# From this section on, there is another way to make the crossovers; that is, by using only half 
# of the chemical formula and randomly selecting the atoms  of each structure.
#----------------------------------------------------------------------------------------------------------
def chem_formula_bisect(xtal_in):
	'''
	This functions splits the chemical formula into two lists of tuples, where each tuple is (simbol,num_atoms)

	in: xtal_in (Molecule), the original structures whose formula is to be split
	out: chem_for
	'''
	big_formula = molecular_stoichiometry(xtal_in,0)
	f_half, s_half = [], []
	for ii in big_formula:
		n1 = ii[1]//2
		n2 = ii[1] - n1
		if n1 != 0 : f_half.append((ii[0], n1))
		if n2 != 0: s_half.append((ii[0], n2))
	return f_half, s_half


#------------------------------------------------------------------------------------------
def build_str_half_comp(xtal_in,chemformula):
	'''
	This function builds a structure, based on the input one, by randomly selecting atoms
	until the desired chemical formula is achieved

	in: xtal_in (Molecule), the original structure that will be reduced in its composition
	    chemformula (list), a list of tuples (symbol,atm_number) with the desired composition
	out: xtal_out (Molecule), the new structure with the desired composition
	'''
	xtal = copymol(xtal_in)
	xtal_out = Molecule(xtal.i, 0.0, xtal.m)
	inatoms=[]
	for ii in chemformula:
		for _ in range(ii[1]):
			inatoms.append(ii[0])
	random.shuffle(inatoms)
	atomlist = xtal.atoms.copy()
	random.shuffle(atomlist)
	for si in inatoms:
		for iatom in atomlist:
			if si==iatom.s:
				xtal_out.add_atom(iatom)
				atomlist.remove(iatom)
				break
	xtal_out=sort_by_stoichiometry([xtal_out])[0]
	return xtal_out

#------------------------------------------------------------------------------------------
def get_complement(xtal_receipt,xtal_donor,comp_stoich,original_stoich):
	'''
	This function builts the crossover between two structures using a different idea that the 
	one in crossover function. Under this scheme, the receiver will contain the half of the 
	composition but the atoms that comprais it will be chosen randomly. The full structure then 
	is built by adding random atoms from the donor.

	in: xtal_receipt (Molecule), the structure that will receive the atoms
	    xtal_donor (Molecule), the donor of random atoms
	    comp_stoich (list), a list of tuples of the missing composition of xtal_receipt
	    original_stoich(list), a list of tuples of the original stoichiometry
	out: new_xtal (Molecule), if the structure was succesfuly built, then it contains it. If it 
	     wasn't, then returns a False value.
	'''
	donor = copymol(xtal_donor)
	new_xtal = copymol(xtal_receipt)
	inatoms=[]
	for iii in comp_stoich:
		for _ in range(iii[1]):
			inatoms.append(iii[0])
			random.shuffle(inatoms)
			atomlist = donor.atoms.copy()
			random.shuffle(atomlist)
	for si in inatoms:
		for iatom in atomlist:
			if si == iatom.s:
				di = min_dist_atm_xtal(iatom,new_xtal)
				if di > 0.5:
					new_xtal.add_atom(iatom)
					atomlist.remove(iatom)
					break
	new_stoich = molecular_stoichiometry(new_xtal,0)
	if new_stoich == original_stoich:
		new_xtal=sort_by_stoichiometry([new_xtal])[0]
	else:
		new_xtal = False
	return new_xtal

#------------------------------------------------------------------------------------------
def new_crossover(base_xtal,complement_xtal):
	base = copymol(base_xtal)
	comp = copymol(complement_xtal)
	org_stoich = molecular_stoichiometry(base,0)
	avg_uc = combined_unit_cell(base,comp)
	avg_base = translating_to_avg_uc(base,avg_uc)
	avg_comp = translating_to_avg_uc(comp,avg_uc)
	f_half,s_half = chem_formula_bisect(avg_base)
	receptor = build_str_half_comp(avg_base,f_half)
	combined = get_complement(receptor,avg_comp,s_half,org_stoich)
	if combined:
		combined.i = base.i + '_XC_' + comp.i
	return combined

#----------------------------------------------------------------------------------------------------------
def many_crossovers(m_list,f_list):
	all_cross = []
	for i, ix in enumerate(m_list):
		for j, jx in enumerate(f_list):
			if i == j and ix != jx:
				r = random.gauss(0,2)
				if r > 0:
					child = crossover(ix,jx)
				else:
					child = new_crossover(ix,jx)
				if child:
					all_cross.append(child)
	return all_cross
#----------------------------------------------------------------------------------------------------------

def popgen_childs(poscarlist, index):
	if number_of_childs==0:
		poscarout=[]
		return poscarout
	logfile = open(log_file,'a')
	print("-------------------------------------------------------------------", file=logfile)
	print("------------------------ MATING: CROSSOVER ------------------------", file=logfile)
	print("CONSTRUCTION OF THE PARENT LIST THROUGH ROULETTE", file=logfile)
	logfile.close()
	mom = get_roulette_wheel_selection(poscarlist, number_of_childs*2)
	pop = get_roulette_wheel_selection(poscarlist, number_of_childs*2)
	poscarout = many_crossovers(mom,pop)
	poscarout = poscarout[0:number_of_childs]
	logfile = open(log_file,'a')
	cont = 1 
	for x in poscarout:
		aux_name = 'mating_' + str(index).zfill(3) + '_' + str(cont).zfill(3)
		print('%s ---> %s' %(aux_name,x.i), file=logfile)
		cont = cont + 1
	print("We have %d POSCAR files type MATING from %d solicited" %(len(poscarout), number_of_childs), file=logfile)
	logfile.close()
	basename = 'mating_' + str(index).zfill(3) + '_'
	poscarout  = rename_molecule(poscarout, basename, 4)
	return poscarout   
