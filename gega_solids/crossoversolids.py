import random
import numpy as np
from utils_solids.libmoleculas import copymol, Molecule, Atom, sort_by_stoichiometry, molecular_stoichiometry, rename_molecule
from utils_solids.miscellaneous import unit_cell_non_negative_coordinates, overlap_check
from utils_solids.translation_rotation import translation_3D
from vasp_solids.libperiodicos import direct2cartesian, cartesian2direct
from gega_solids.solids_roulette import get_roulette_wheel_selection
from inout_solids.getbilparam import get_a_int, get_a_str
#----------------------------------------------------------------------------------------------------------
number_of_childs = get_a_int('number_of_matings',10)
log_file = get_a_str('output_file','solids_out.txt')
#----------------------------------------------------------------------------------------------------------
def bounded_random_point(r_max,r_min):
	'''Finds random_point such as r_min < random_point < r_max

	in: 
	r_max, r_min (float); Top and Lower limits for the random number

	out:
	random_point (float); r_min < random_point < r_max 
	'''
	random_point = 1
	if r_min > r_max:
		random_point = 0.4
		return random_point
	while random_point > r_max or random_point < r_min:
		random_point = random.random()
	random_point = round(random_point,2)
	return random_point

#----------------------------------------------------------------------------------------------------------
def rot(xtal_in):
	'''Randomly exchange the unit cell vectors. ------------------Does it really work?

	in:
	xtal_in (Molecule); The structure whose vectors will be exchanged

	out:
	xtal_out (Molecule); The new structure
	'''
	xtal_out = copymol(xtal_in)
	org = list(xtal_out.m)
	random.shuffle(org)
	new = np.array(org)
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
	xtal_out = Molecule(cp_xtal.i, cp_xtal.e, cp_xtal.m)
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
def atms_correction_above_cut(xtal_in):
	xtal_out = copymol(xtal_in)
	for a in xtal_out.atoms:
		a.xc,a.yc,a.zc = cartesian2direct(a.xc,a.yc,a.zc,xtal_out.m)
		if a.xc == 0:
			a.xc = 1
		if a.yc == 0:
			a.yc = 1 
		if a.zc == 0:
			a.zc = 1
		a.xc,a.yc,a.zc = direct2cartesian(a.xc,a.yc,a.zc,xtal_out.m)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def parent_cut_above(xtal_in, axis, cutting_point):
	'''
	This function is exactly the same than parent_cut_bellow but cutting above the cutting_point
	'''
	cp_xtal = copymol(xtal_in)
	cp_xtal = atms_correction_above_cut(cp_xtal)
	new_xtal = Molecule(cp_xtal.i, cp_xtal.e, cp_xtal.m)
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
def combined_unit_cell(xtal_a, xtal_b, weight_h=0.6, weight_s=0.4):
	'''Averages the unit cells of parent structures, a weighted-averaged unit cell is returned.
	
	in: 
	xtal_a, xtal_b (Molecule); two structures whose UC will be averaged
	weight_h (float); the higher weight to be given to the vectors
	weight_s (float); the lower weigth to be given to the vectors
	
	out: 
	new_mtx (Numpy Array); weighted averaged UC
	'''
	mtx_a, mtx_b = xtal_a.m, xtal_b.m
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
def translating_to_avg_uc(xtal_in,avg_uc):
	'''Translates a given structure to the corresponding one in the weighted-averaged UC

	in: 
	xtal_in (Molecule); the original crystal structure that will be translated
	avg_uc (Numpy array); th UC upon wich the translation will be done.
	
	out: 
	xtal_out (Molecule); the new translated crystal structure 
	'''
	cp_xtal = copymol(xtal_in)
	xtal_out = Molecule(cp_xtal.i,cp_xtal.e,avg_uc)
	for iatom in cp_xtal.atoms:
		x,y,z = cartesian2direct(iatom.xc,iatom.yc,iatom.zc,cp_xtal.m)
		s = iatom.s
		xn,yn,zn = direct2cartesian(x,y,z,avg_uc)
		new_atom = Atom(s,xn,yn,zn)
		xtal_out.add_atom(new_atom) 
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def missing_atoms_identifier(original_stoich,new_stoich):
	'''Compares the original composition of the structure and gets the amount of missing atoms and their species.

	in: 
	original_stoich (list[tuples]); a list that has all tuples with (symbol,#atoms)
	new_stoich (list[tuples]); the composition after the cut(symbol,#atoms)
	
	out: 
	missing_atoms[list[tuples]]; the amount of atoms missing on each symbol (symbol,#atoms) 
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
def crossover(base_xtal,complement_xtal,ref_d):
	'''Transfer the genetic information of two parent structures to a new one, called child

	in: 
	base_xtal, complement_xtal (Molecule); Two parent structures for gene transmition

	out: 
	xtal_out (Molecule/False); Child structure, False if crossover is not possible
	'''
	#from vasp_solids.libperiodicos import readposcars, writeposcars
	base = copymol(base_xtal)
	base = unit_cell_non_negative_coordinates(base)	
	org_stoich = molecular_stoichiometry(base,0)
	comp = copymol(complement_xtal)
	comp = unit_cell_non_negative_coordinates(comp)
	rw1 = random.random()
	rw2 = 1 - rw1
	avg_uc = combined_unit_cell(base, comp,rw1,rw2)
	avg_base = translating_to_avg_uc(base,avg_uc)
	
	#writeposcars([avg_base],'base.vasp','D')

	avg_comp = translating_to_avg_uc(comp,avg_uc)
	
	#writeposcars([avg_comp],'comp.vasp','D')

	general_stop = 0
	wflag = False
	while wflag == False:
		r_vect = random.choice(['a','b','c'])
		r_point = bounded_random_point(0.7,0.4)
		avg_base = translation_3D(avg_base)

		#writeposcars([avg_base],'base_trans.vasp','D')

		avg_comp = translation_3D(avg_comp)

		#writeposcars([avg_comp],'comp_trans.vasp','D')


		xtal_out = parent_cut_bellow(avg_base,r_vect,r_point)

		#writeposcars([xtal_out],'cut_base.vasp','D')

		cut_comp = parent_cut_above(avg_comp,r_vect,r_point)

		#writeposcars([cut_comp],'cut_comp.vasp','D')

		new_stoich = molecular_stoichiometry(xtal_out,0)
		comp_stoich = molecular_stoichiometry(cut_comp,0)
		m_atm = missing_atoms_identifier(org_stoich,new_stoich)
		av_atms, ms_atms = 0,0
		for ta,tm in zip(comp_stoich,m_atm):
			av_atms = av_atms + ta[1]
			ms_atms = ms_atms + tm[1]
		if av_atms >= ms_atms:
			aux_list = cut_comp.atoms.copy()
			random.shuffle(aux_list)
			for t in m_atm:
				ms, ma = t[0],t[1]
				for i in range(ma):
					if aux_list:
						r_atm = random.choice(aux_list)
						x,y,z = r_atm.xc,r_atm.yc, r_atm.zc 
						natm = Atom(ms,x,y,z)
						oc = overlap_check(natm,xtal_out,ref_d)
						if oc == False:
							xtal_out.add_atom(natm)
							aux_list.remove(r_atm)
		else:
			continue
		xtal_out = sort_by_stoichiometry([xtal_out])[0]
		out_stoich = molecular_stoichiometry(xtal_out,0)
		if out_stoich == org_stoich:
			wflag = True
			xtal_out.i = str(base_xtal.i) + '_x_' + str(complement_xtal.i)
		else:
			general_stop = general_stop + 1
		if general_stop == 3:
			xtal_out = False
			break

		#writeposcars([xtal_out],'result.vasp','D')

	return xtal_out

# from vasp_solids.libperiodicos import readposcars, writeposcars

# x = readposcars('1.vasp')
# y = readposcars('2.vasp')
# z = crossover(x,y,[(Mg,Mg,1.0),(Mg,Al,1.0),(Mg,O,1.0),(Al,Al,1.0),(Al,O,1.0),(O,O,1.0)])

#----------------------------------------------------------------------------------------------------------
def many_crossovers(m_list,f_list,ref_d):
	all_cross = []
	for i, ix in enumerate(m_list):
		for j, jx in enumerate(f_list):
			if i == j:
				child = crossover(ix,jx,ref_d)
				if child:
					all_cross.append(child)
	return all_cross

#----------------------------------------------------------------------------------------------------------
def popgen_childs(poscarlist, ref_d, index):
	if number_of_childs==0:
		poscarout=[]
		return poscarout
	logfile = open(log_file,'a')
	print("-------------------------------------------------------------------", file=logfile)
	print("------------------------ MATING: CROSSOVER ------------------------", file=logfile)
	print("CONSTRUCTION OF THE PARENT LIST THROUGH ROULETTE", file=logfile)
	logfile.close()
	c,c2 = 0,0
	while c < number_of_childs:
		mom = get_roulette_wheel_selection(poscarlist, number_of_childs)
		pop = get_roulette_wheel_selection(poscarlist, number_of_childs)
		poscarout = many_crossovers(mom,pop,ref_d)
		c = c + len(poscarout)
		c2 = c2 + 1
		if c2 > 15:
			break
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
	poscarout  = rename_molecule(poscarout, basename, 3)
	return poscarout
