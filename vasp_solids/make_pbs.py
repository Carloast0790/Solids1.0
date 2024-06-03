import os
import os.path
#------------------------------------------------------------------------------------------
def make_apbs_vasp(ppj, queue_name, walltime, path, basename, conf):
    filepbs=basename+'.pbs'
    fileout=basename+'.out'
    if not os.path.isfile(path+fileout):
        fh=open(filepbs,'w')
        print("#!/bin/bash", file=fh)
        print("#PBS -N %s" %(basename), file=fh)
        print("#PBS -l nodes=1:ppn=%d" %(ppj), file=fh)
        print("#PBS -q %s" %(queue_name), file=fh)
        print("#PBS -l walltime=%s\n" %(walltime), file=fh)
        #-------------------------------------------------------
        for ii in conf:
            if ii=='RUNLINE\n':
                print("cd $PBS_O_WORKDIR/%s" %(path), file=fh)
                print("$exe_mpirun -np %d $exe_vasp > $PBS_O_WORKDIR/%s%s" %(ppj, path, fileout), file=fh)
            else:
                fh.write(ii)
        #-------------------------------------------------------
        fh.close()
    else:
        os.system("rm -f "+filepbs)
#------------------------------------------------------------------------------------------
