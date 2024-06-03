import os.path
#------------------------------------------------------------------------------------------
def make_ash_gulp(basename, conf_gulp, path):
    filesh=basename+'.sh'
    fileinp=basename+'.gin'
    fileout=basename+'.got'
    if not os.path.isfile(path+fileout):
        fh=open(filesh,'w')
        print("#!/bin/bash", file=fh)
        for ii in conf_gulp: fh.write(ii)
        print("wait $PID ; $exe_gulp < %s > %s" %(path+fileinp,path+fileout), file=fh)
        fh.close()
    else:
        os.system("rm -f -v "+filesh)
#------------------------------------------------------------------------------------------
def run_sample():
    from inout.readbil import read_block_of_bil
    conf_gulp=read_block_of_bil('gulp.conf')
    for ii in conf_gulp: print(ii)
    make_ash_gulp('anatasa', conf_gulp, 'ejemplos/')
#run_sample()
#------------------------------------------------------------------------------------------
