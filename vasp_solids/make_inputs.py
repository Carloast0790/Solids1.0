import os.path
from inout_solids.getbilparam import get_a_int, get_str_list, get_a_str
from utils_solids.libmoleculas import copymol
from vasp_solids.libperiodicos import writeposcars, molecule2poscar
#------------------------------------------------------------------------------------------
restricted=get_str_list('restricted_atoms',[])
firstatoms=get_a_int('restricted_first',0)
ufirstatoms=get_a_int('unrestricted_first',0)
incar_list=get_str_list('incar_files'  ,[])
kpoints_list=get_str_list('kpoints_files',[])
#------------------------------------------------------------------------------------------
def make_a_poscar(moleculein, path, opt):
    nameinp='POSCAR'
    posx=molecule2poscar(moleculein) if len(moleculein.m) == 0 else copymol(moleculein)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if len(restricted)==0 and firstatoms==0 and ufirstatoms==0:
        for iatom in posx.atoms:        iatom.xf, iatom.yf, iatom.zf='T', 'T', 'T'
    elif len(restricted)>0:
        for iatom in posx.atoms:
            if (iatom.s in restricted): iatom.xf, iatom.yf, iatom.zf='F', 'F', 'F'
            else:                       iatom.xf, iatom.yf, iatom.zf='T', 'T', 'T'
    elif firstatoms>0:
        for ii, iatom in enumerate(posx.atoms):
            if (ii < firstatoms):       iatom.xf, iatom.yf, iatom.zf='F', 'F', 'F'
            else:                       iatom.xf, iatom.yf, iatom.zf='T', 'T', 'T'
    elif ufirstatoms>0:
        for ii, iatom in enumerate(posx.atoms):
            if (ii>=ufirstatoms):       iatom.xf, iatom.yf, iatom.zf='F', 'F', 'F'
            else:                       iatom.xf, iatom.yf, iatom.zf='T', 'T', 'T'
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    writeposcars([posx],path+nameinp, opt)
#------------------------------------------------------------------------------------------
def make_ainput_vasp(moleculein, path, stage, opt):
    nameout='OUTCAR'
    if os.path.isfile(path+nameout):
        log_file =get_a_str('output_file','solids_out.txt')
        fopen = open(log_file,'a')
        print("%sInputFiles were not built because %s exist" %(path, nameout), file=fopen)
        fopen.close()
    else:
        if not os.path.isfile(path+'POSCAR'):
            make_a_poscar(moleculein, path, opt)
        if not os.path.isfile(path+'INCAR'):
            string="cp %s %sINCAR" %(incar_list[stage], path)
            os.system(string)
        if not os.path.isfile(path+'KPOINTS'):
            string="cp %s %sKPOINTS" %(kpoints_list[stage], path)
            os.system(string)
        if not os.path.isfile(path+'POTCAR'):
            string='cp POTCAR '+path
            os.system(string)
#------------------------------------------------------------------------------------------
