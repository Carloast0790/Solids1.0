import random
import numpy as np
from utils.libmoleculas import copymol
from vasp.libperiodicos import writeposcars, direct2cartesian, cartesian2direct, readposcars

#----------------------------------------------------------------------------------------------------------
def translation_x(xtal_in,delta_x):
	"""
	This function as well as 'translation_y' and 'translation_z' translates the unit cell across 
	the x/y/z-axis a 'delta_x/y/z' distance.

	in: xtal_in (Molecule object), the crystal from wich the unit cell will be translated
	    delta_x/y/z (float), the amount the unit cell will be translated

	out: xtal_out (Molecule object), the new crystal with the unit cell translated
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
	"""
	This function randomly translates the unit cell in the same amount in all three dimensions of crystal

	in: xtal_in (Molecule object), the crystal from wich the unit cell will be translated

	out: xtal_out (Molecule object), the new crystal with the unit cell translated  
	"""
	xtal_out0 = copymol(xtal_in)
	delta = random.random()
	xtal_out1 = translation_x(xtal_out0,delta)
	xtal_out2 = translation_y(xtal_out1,delta)
	xtal_out = translation_z(xtal_out2,delta)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def rotation_matrix_transformation(lattice_matrix,rand_degree):
	cos, sin = np.cos(rand_degree), np.sin(rand_degree)
	A = np.array([[cos, -sin, 0],
		  [sin, cos, 0],
		  [0,0,1]])
	v = np.array([lattice_matrix[0],lattice_matrix[1],lattice_matrix[2]])
	v_inv = np.linalg.inv(v)
	m = np.dot(v_inv,A)
	transformed_rotation_matrix = np.dot(m,v)
	return(transformed_rotation_matrix)

#----------------------------------------------------------------------------------------------------------
def rotation_z(xtal_in):
	xtal_out = copymol(xtal_in)
	xtal_out = translation_3D(xtal_out)
	r = random.random()
	print('se rotó',r)
	r = np.pi * r
	rotation_matrix = rotation_matrix_transformation(xtal_out.m,r)
	aux = 0
	for iatom in xtal_out.atoms:
		p_direct = np.array(cartesian2direct(iatom.xc,iatom.yc,iatom.zc,xtal_out.m))
		n_vector = np.dot(rotation_matrix,p_direct)
		if n_vector[0] > 1 or n_vector[1] > 1 or n_vector[2] > 1:
			print(iatom.s,'está fuera de la celda en x o y o z > 1')
			del xtal_out.atoms[aux]
		aux = aux + 1
		iatom.xc, iatom.yc, iatom.zc = direct2cartesian(n_vector[0],n_vector[1],n_vector[2],xtal_out.m)
	return xtal_out

#----------------------------------------------------------------------------------------------------------
def rotation2D(xtal_in):
	xtal_out = copymol(xtal_in)
	r = random.randint(2,10)
	deg = np.pi / r
	cos, sin = np.cos(deg), np.sin(deg)
	m = np.array([[cos,-sin],[sin,cos]])
	for a in xtal_out.atoms:
		xd, yd = a.xc,a.yc
		v = np.array([xd,yd])
		nv = np.matmul(v,m)
		nx,ny = nv[0],nv[1]
		a.xc, a.yc = nx,ny	
	return xtal_out

# from vasp.libperiodicos import readposcars, writeposcars
# x = readposcars('rutilo.vasp')[0]
# y = translation_3D(x)
# writeposcars([y],'translation.vasp','D')
# z = rotation2D(y)
# z.i = 'rot01'
# writeposcars([z],'rotation.vasp','D')
