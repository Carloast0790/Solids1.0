import os.path
from inout_solids.getbilparam  import get_a_str, get_a_int, get_a_float, get_str_list, get_int_list, get_float_list
from vasp_solids.make_potcar   import make_potcar_file
from vasp_solids.make_inputs   import make_ainput_vasp
from vasp_solids.get_geometry  import get_normaltermination_vasp, get_all_nt_geometry_vasp, get_all_at_geometry_vasp, get_all_xt_geometry_vasp
from vasp_solids.clean         import clean_vasp
from discriminate_solids.usrp  import molin_sim_molref, kick_similar_molecules
from utils_solids.libmoleculas import readxyzs
#------------------------------------------------------------------------------------------
listnoa  =get_int_list('no_attempts_opt',[1, 1])
listper  =get_float_list('percent_of_convergence',[90.0, 90.0])
tol_opt  =get_a_float('similarity_tol', 0.95)
log_file =get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
qsys=get_a_str('qsys','pbs')
if qsys=='pbs':
    from vasp.make_pbs import make_apbs_vasp
    from queuing.qpbs  import send_pbs_files_to_queue
    from inout.readbil import read_block_of_bil
    ppj      =get_a_int('nprocshared',4)
    queue    =get_a_str('queue','qintel')
    walltime =get_str_list('walltime','[04:00:00, 08:00:00]')
    njobs    =get_a_int('njobs',4)
    timesleep=get_a_float('timesleep',1.0)
    conf_vasp=read_block_of_bil('vasp.conf')
#------------------------------------------------------------------------------------------
def calculator_vasp_all(moleculein, folder, stage=0):
    wt=walltime[stage]
    if not os.path.exists(folder): os.system('mkdir '+folder)
    make_potcar_file()
    for imol in moleculein:
        foldername=imol.i
        path=folder+foldername+'/'
        if not os.path.exists(path): os.system('mkdir '+path)
        make_ainput_vasp(imol, path, stage, 'D')
        make_apbs_vasp(ppj, queue, wt, path, foldername, conf_vasp)
    if njobs == 0:
        fopen = open(log_file,'a')
        print("Number of Jobs is zero. FINISH", file=fopen)
        fopen.close()
        exit()
    send_pbs_files_to_queue(njobs,timesleep)
    for imol in moleculein:
        foldername=imol.i
        path=folder+foldername+'/'
        clean_vasp(path, foldername)
#------------------------------------------------------------------------------------------
def calculator_vasp_kernel(moleculein, folder, stage=0):
    success,count=0,1
    n0=len(moleculein)
    noa=listnoa[stage]
    per=listper[stage]
    while success==0:
        fopen = open(log_file,'a')
        print("---------------------------------------- BEGIN CALC-ATTEMP %d" %(count), file=fopen)
        fopen.close()
        if count==1:
            calculator_vasp_all(moleculein, folder, stage)
            moleculeout=get_all_nt_geometry_vasp(folder, moleculein)
        else:
            calculator_vasp_all(moleculeux, folder, stage)
            moleculetmp=get_all_nt_geometry_vasp(folder, moleculeux)
            moleculeout.extend(moleculetmp)
            moleculeux.clear()
        n1=len(moleculeout) if moleculeout != [] else 0
        p=float(n1)*100.0/float(n0)
        if (p < float(per)):
            fopen = open(log_file,'a')
            print("Calculations with NT = %2.1f percent, less than the requested (%2.1f percent)" %(p,per), file=fopen)
            fopen.close()
            if count==noa:
                fopen = open(log_file,'a')
                print("We have reached the maximum number of attempts (%d)." %(noa), file=fopen)
                print("------------------------------------------ END CALC-ATTEMP %d" %(count), file=fopen)
                fopen.close()
                break
            else:
                erase=0
                moleculeux=get_all_at_geometry_vasp(folder, moleculein, erase)
                if moleculeux != []:
                    n2a=len(moleculeux)
                    moleculeux=molin_sim_molref(moleculeux, moleculeout, tol_opt, 0)
                    memory=readxyzs('memory.xyz')
                    moleculeux=molin_sim_molref(moleculeux, memory, tol_opt, 0)

                    n2d=len(moleculeux)
                    n0=n0-(n2a-n2d)
                    if moleculeux != []:
                        n3a=len(moleculeux)
                        moleculeux=kick_similar_molecules(moleculeux, tol_opt, 0)
                        n3d=len(moleculeux)
                        n0=n0-(n3a-n3d)
                        erase=1
                        moleculeux=get_all_at_geometry_vasp(folder, moleculeux, erase)
                else:
                    fopen = open(log_file,'a')
                    print("None AT calculation provide structural information", file=fopen)
                    print("------------------------------------------ END CALC-ATTEMP %d" %(count), file=fopen)
                    fopen.close()
                    if n1==0: exit()
                    else: break
            fopen = open(log_file,'a')
            print("------------------------------------------ END CALC-ATTEMP %d" %(count), file=fopen)
            fopen.close()
            count=count+1
        else:
            fopen = open(log_file,'a')
            print("Calculations with NT = %2.1f percent satisfying the requested (%2.1f percent)" %(p,per), file=fopen)
            print("------------------------------------------ END CALC-ATTEMP %d" %(count), file=fopen)
            fopen.close()
            success=1
#------------------------------------------------------------------------------------------
def calculator_vasp_all_check(moleculein, folder='generation000/', stage=0):
    list_ne, list_at, list_nt=[],[],[]
    n0=len(moleculein)
    per=listper[stage]
    for imol in moleculein:
        foldername=imol.i
        path=folder+foldername+'/'
        id=get_normaltermination_vasp(path)
        if (id is False): list_ne.append(foldername)
        elif (id == 0):   list_at.append(foldername)
        else:             list_nt.append(foldername)
    nnt, nat, nne = len(list_nt), len(list_at), len(list_ne)
    pnt=float(nnt)*100.0/float(n0)
    if (pnt >= float(per)):
        fopen = open(log_file,'a')
        print("Calculations with N.T. = %d" %(nnt), file=fopen)
        print("Calculations with A.T. = %d" %(nat), file=fopen)
        print("Calculations with NT = %2.1f percent satisfies the requested (%2.1f percent)" %(pnt,per), file=fopen)
        fopen.close()
    if (nnt+nat > 0 and pnt < float(per)):
        fopen = open(log_file,'a')
        print("Calculations with N.T. = %d" %(nnt), file=fopen)
        print("Calculations with A.T. = %d" %(nat), file=fopen)
        print("Calculations with NT = %2.1f percent less than the requested (%2.1f percent)" %(pnt,per), file=fopen)
        fopen.close()
    ctn=get_a_int('cont_calc',0)
    if (nne == n0) or (ctn == int(1)):
        calculator_vasp_kernel(moleculein, folder, stage)
    moleculeout=get_all_xt_geometry_vasp(folder, moleculein)
    if moleculeout is False:
        fopen = open(log_file,'a')
        print ("ZERO calculation finished satisfactorily", file=fopen)
        fopen.close()
        exit()
    return moleculeout
#------------------------------------------------------------------------------------------
