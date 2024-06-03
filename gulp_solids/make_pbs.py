import os.path
#------------------------------------------------------------------------------------------
def make_apbs_gulp(basename, ppj, ram, queue, walltime, conf_gulp, path):
    filepbs=basename+'.pbs'
    fileinp=basename+'.gin'
    fileout=basename+'.got'
    if not os.path.isfile(path+fileout):
        fh=open(filepbs,'w')
        print("#!/bin/bash", file=fh)
        print("#PBS -N %s" %(basename), file=fh)
        print("#PBS -l nodes=1:ppn=%d" %(ppj), file=fh)
        print("#PBS -l mem=%dgb" %(ram), file=fh)
        print("#PBS -q %s" %(queue), file=fh)
        print("#PBS -l walltime=%s\n" %(walltime), file=fh)
        print("cd $PBS_O_WORKDIR\n", file=fh)
        #-------------------------------------------------------
        for ii in conf_gulp: fh.write(ii)
        #-------------------------------------------------------
        print("$exe_gulp < $PBS_O_WORKDIR/%s > $PBS_O_WORKDIR/%s\n" %(path+fileinp,path+fileout), file=fh)
        fh.close()
    else:
        os.system("rm -f -v "+filepbs)
#------------------------------------------------------------------------------------------
def run_sample():
    from inout.getbilparam  import get_a_int, get_a_str, get_str_list, get_a_float
    from inout.readbil import read_block_of_bil
    ppj = get_a_int('nprocshared',1)
    ram = get_a_int('memory_in_gb',1)
    queue = get_a_str('queue','qintel')
    wallt = get_str_list('walltime','[01:00:00, 01:00:00]')
    conf_gulp=read_block_of_bil('gulp.conf')
    print("nprocshared  = %d" %(ppj))
    print("memory_in_gb = %d" %(ram))
    print("queue        = %s" %(queue))
    print("walltime     = %s" %(wallt))
    for ii in conf_gulp: print(ii)
    make_apbs_gulp('anatasa', ppj, ram, queue, wallt[0], conf_gulp, 'ejemplos/')
#run_sample()
#------------------------------------------------------------------------------------------
