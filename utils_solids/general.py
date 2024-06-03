import random
import numpy as np
#------------------------------------------------------------------------------------------
def randunitvector(option=1):
#>>> from utils.general import randunitvector
    if option == 1:
        ##SPHERE DISTRIBUTION
        phi=float(random.uniform(0.0, 2.0*(np.pi)))
        theta=float(random.uniform(0.0,(np.pi)))
    if option == 2:
        ##DISC DISTRIBUTION
        phi=float(random.uniform(0.0, 2.0*(np.pi)))
        theta=np.pi/2.0
    if option == 3:
        ##LINEAR DISTRIBUTION
        phi=float(random.uniform(0.0, 2.0*(np.pi)))
        theta=float(0.0)
    xu=np.sin(theta) * np.cos(phi)
    yu=np.sin(theta) * np.sin(phi)
    zu=np.cos(theta)
    vector = np.array([xu, yu, zu])
    return vector
#------------------------------------------------------------------------------------------
def randvector(r, option=1):
    ru=float(random.uniform(0.0, r))
    vector = ru*randunitvector(option)
    return vector
#------------------------------------------------------------------------------------------
def atomrandunitvector(moleculein):
#>>> from utils.general import atomrandunitvector
    lista=range(moleculein.n)
    random.shuffle(lista)
    iatom=lista.pop()
    si=moleculein.atoms.s
    xi=moleculein.atoms.xc
    yi=moleculein.atoms.yc
    zi=moleculein.atoms.zc
    v=np.array([xi, yi, zi])
    a=np.linalg.norm(v)
    u=v/a
    return si, v, u
#------------------------------------------------------------------------------------------
def euler_matrix(qdegx, qdegy, qdegz):
    qradx=float(qdegx)*(np.pi)/180.0
    qrady=float(qdegy)*(np.pi)/180.0
    qradz=float(qdegz)*(np.pi)/180.0
    ##qradz=np.deg2rad(qdegz)
    cx,cy,cz=np.cos(qradx),np.cos(qrady),np.cos(qradz)
    sx,sy,sz=np.sin(qradx),np.sin(qrady),np.sin(qradz)
    row1=[cy*cz,-cy*sz,sy]
    row2=[cx*sz+cz*sx*sy,cx*cz-sx*sy*sz,-cy*sx]
    row3=[sx*sz-cx*cz*sy,cz*sx+cx*sy*sz, cx*cy]
    eulerm=np.array([row1,row2,row3])
    return eulerm
#------------------------------------------------------------------------------------------
def rodrigues_rotation_matrix(kvector, qdeg):
    qrad=float(qdeg)*np.pi/180.0
    kvec=np.array(kvector)
    kuv=kvec/np.linalg.norm(kvec)
    kmat1 = np.array([[0.0, -kuv[2], kuv[1]], [kuv[2], 0.0, -kuv[0]], [-kuv[1], kuv[0], 0.0]])
    kmat2 = np.matmul(kmat1,kmat1)
    rodriguesrm = np.eye(3) + np.sin(qrad)*kmat1 + (1.0 - np.cos(qrad))*kmat2
    return rodriguesrm
#------------------------------------------------------------------------------------------
###Calculate Rotation Matrix to align Vector A to Vector B in 3d
def rotation_matrix_that_rotates_uva_onto_uvb(unit_vector_a, unit_vector_b):
    v = np.cross(unit_vector_a, unit_vector_b)
    s = np.linalg.norm(v)
    c = np.dot(unit_vector_a, unit_vector_b)
    vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    vxvx = np.matmul(vx,vx)
    r = np.eye(3) + vx + vxvx*(1.0-c)/(s*s)
    return r
#------------------------------------------------------------------------------------------
def run_sample():
    v1=randunitvector()
    v2=randunitvector()
    rmatrix=rotation_matrix_that_rotates_uva_onto_uvb(v1,v2)
    print(v1)
    print(v2)
    v3=np.matmul(rmatrix,v1)
    print(v3)
#run_sample()
#------------------------------------------------------------------------------------------
