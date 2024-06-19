import random
import numpy as np
from utils_solids.libmoleculas import copymol
from vasp_solids.libperiodicos import direct2cartesian, cartesian2direct
#----------------------------------------------------------------------------------------------------------
def translation_x(xtal_in,delta_x):
	"""This function as well as 'translation_y' and 'translation_z' translates the unit cell across 
	the x/y/z-axis a 'delta_x/y/z' distance.

	in: 
	xtal_in (Molecule); the crystal from wich the unit cell will be translated
	delta_x/y/z (float); the amount the unit cell will be translated

	out: 
	xtal_out (Molecule); the new crystal with the unit cell translated
	"""
	xtal_out = copymol(xtal_in)
	for iatom in xtal_out.atoms:
		iatom.xc, iatom.yc, iatom.zc = cartesian2direct(iatom.xc, iatom.yc, iatom.zc,xtal_out.m)
		if abs(iatom.xc) == 0:
			iatom.xc = 1 - delta_x
		elif delta_x == iatom.xc:
			iatom.xc = 0
		elif delta_x > iatom.xc:
			iatom.xc = 1 - delta_x + iatom.xc
		else:
			iatom.xc = iatom.xc - delta_x
	for iatom in xtal_out.atoms:
		iatom.xc,iatom.yc,iatom.zc = direct2cartesian(iatom.xc,iatom.yc,iatom.zc,xtal_out.m)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def translation_y(xtal_in, delta_y):
	xtal_out = copymol(xtal_in)
	for iatom in xtal_out.atoms:
		iatom.xc, iatom.yc, iatom.zc = cartesian2direct(iatom.xc, iatom.yc, iatom.zc,xtal_out.m)
		if abs(iatom.yc) == 0:
			iatom.yc = 1 - delta_y
		elif delta_y == iatom.yc:
			iatom.yc = 0
		elif delta_y > iatom.yc:
			iatom.yc = 1 - delta_y + iatom.yc
		else:
			iatom.yc = iatom.yc - delta_y
	for iatom in xtal_out.atoms:
		iatom.xc,iatom.yc,iatom.zc = direct2cartesian(iatom.xc,iatom.yc,iatom.zc,xtal_out.m)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def translation_z(xtal_in, delta_z):
	xtal_out = copymol(xtal_in)
	for iatom in xtal_out.atoms:
		iatom.xc, iatom.yc, iatom.zc = cartesian2direct(iatom.xc, iatom.yc, iatom.zc,xtal_out.m)
		if abs(iatom.zc) == 0:
			iatom.zc = 1 - delta_z
		elif delta_z == iatom.zc:
			iatom.zc = 0
		elif delta_z > iatom.zc:
			iatom.zc = 1 - delta_z + iatom.zc
		else:
			iatom.zc = iatom.zc - delta_z
	for iatom in xtal_out.atoms:
		iatom.xc,iatom.yc,iatom.zc = direct2cartesian(iatom.xc,iatom.yc,iatom.zc,xtal_out.m)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def translation_3D(xtal_in):
	"""This function randomly translates the unit cell in the same amount in all three dimensions of crystal

	in:
	xtal_in (Molecule); the crystal from wich the unit cell will be translated

	out: 
	xtal_out (Molecule); the new crystal with the unit cell translated  
	"""
	xtal_out0 = copymol(xtal_in)
	delta = random.random()
	xtal_out1 = translation_x(xtal_out0,delta)
	xtal_out2 = translation_y(xtal_out1,delta)
	xtal_out = translation_z(xtal_out2,delta)
	return xtal_out


