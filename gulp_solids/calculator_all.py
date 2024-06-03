import os.path
import numpy as np
from inout_solids.readbil      import read_block_of_bil
from inout_solids.getbilparam  import get_a_int, get_a_float, get_a_str, get_str_list
from gulp_solids.make_inputs   import make_ainput_gulp
from gulp_solids.get_geometry  import get_normaltermination_gulp, get_energy_gulp, get_all_xt_geometry_gulp
from vasp_solids.libperiodicos import molecule2poscar
#------------------------------------------------------------------------------------------
qsys=get_a_str('qsys','pbs')
if qsys=='pbs':
    from gulp_solids.make_pbs   import make_apbs_gulp
    from queuing_solids.qpbs    import send_pbs_files_to_queue
    ppj = get_a_int('nprocshared',1)
    ram = get_a_int('memory_in_gb',1)
    queue = get_a_str('queue','qintel')
    njobs = get_a_int('njobs',4)
    wallt = get_str_list('walltime','[00:30:00, 00:30:00]')
    sleep = get_a_float('timesleep',1.0)
if qsys=='local':
    from gulp_solids.make_sh     import make_ash_gulp
    from queuing_solids.qlocal   import send_sh_files_to_local
#------------------------------------------------------------------------------------------
def calculator_gulp_queue(poscarlist, folder='./', blockname='gulp', stage=0):
    headgulp=read_block_of_bil(blockname)
    conf_gulp=read_block_of_bil('gulp.conf')
    if not os.path.exists(folder): os.system('mkdir %s' %(folder))
    for iposcar in poscarlist:
        if len(iposcar.m) == 0: iposcar=molecule2poscar(iposcar)
        make_ainput_gulp(headgulp, iposcar, folder)
        if qsys=='pbs':   make_apbs_gulp(iposcar.i, ppj, ram, queue, wallt[stage], conf_gulp, folder)
        if qsys=='local': make_ash_gulp(iposcar.i, conf_gulp, folder)
    if qsys=='pbs':   send_pbs_files_to_queue(njobs, sleep)
    if qsys=='local': send_sh_files_to_local()
#------------------------------------------------------------------------------------------
def calculator_gulp_all_check(poscarlist, folder='./', blockname='gulp', stage=0):
    if not os.path.exists(folder):
        calculator_gulp_queue(poscarlist, folder, blockname, stage)
    list_ne, list_at, list_nt, listall=[],[],[],[]
    n0=len(poscarlist)
    for iposcar in poscarlist:
        iname=iposcar.i
        listall.append(iname)
        id=get_normaltermination_gulp(folder, iname+'.got')
        if (id is False): list_ne.append(iname)
        elif (id == 0):   list_at.append(iname)
        else:             list_nt.append(iname)
    nnt, nat, nne = len(list_nt), len(list_at), len(list_ne)
    pnt=int(float(nnt)*100.0/float(n0))
    pat=int(float(nat)*100.0/float(n0))
    log_file =get_a_str('output_file','solids_out.txt')
    fopen = open(log_file,'a')
    print("Calculations with N.T. = %d (%d percent)" %(nnt,pnt), file=fopen)
    print("Calculations with A.T. = %d (%d percent)" %(nat,pat), file=fopen)
    fopen.close()
    moleculeout=get_all_xt_geometry_gulp(poscarlist, folder)
    if moleculeout is False:
        fopen = open(log_file,'a')
        print("ZERO calculation finished satisfactorily", file=fopen)
        fopen.close()
        exit()
    return moleculeout
#------------------------------------------------------------------------------------------
def run_sample():
    from vasp_solids.libperiodicos import readposcars, writeposcars, expand_poscar
    from utils_solids.libmoleculas import readxyzs, writexyzs, rename_molecule, sort_by_energy
    #poscars=readposcars('initial.vasp')
    poscars=readxyzs('initial.xyz')
    poscars=rename_molecule(poscars, 'test', 3)
    ###calculator_gulp_queue(poscars, 'stage0/', 'gulp', 0)
    moleculeoutlist=calculator_gulp_all_check(poscars, 'stage0/', 'gulp', 0)
    moleculeoutlist=sort_by_energy(moleculeoutlist,0)
    writeposcars(moleculeoutlist, 'stage0.vasp', 'D')
    poscarliste=expand_poscar(moleculeoutlist, 1, 1, 1)
    writexyzs(poscarliste, 'stage0.xyz')
#run_sample()
