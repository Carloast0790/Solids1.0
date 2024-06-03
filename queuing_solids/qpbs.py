from subprocess import Popen, PIPE
from sys import stdout
import time
import glob
import os.path
from inout_solids.getbilparam import get_a_str
log_file=get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
def queued_processes():
    workid=[]
    qs = os.popen('qstat')
    lines = qs.readlines()
    for line in lines[2:]:
        pid,job,user,time,status,queue = line.split()
        workid.append(pid)
    return workid
#------------------------------------------------------------------------------------------
def qsub(args, iii):
    cmd = ['qsub', args]
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    output, errmsg = p.communicate()
    if p.returncode != 0:
        print("command %r failed: %s" %(" ".join(cmd), errmsg.decode('ascii')))
        anspbs=False
    else:
        anspbs=output.strip()
        anspbs=anspbs.decode('ascii')
        ii=str(iii+1).zfill(5)
        mytime=time.strftime("%c")
        fopen = open(log_file,'a')
        print("#%s Date %s --> %s JobID: %s" %(ii, mytime, args, anspbs), file=fopen)
        fopen.close()
    return anspbs
#------------------------------------------------------------------------------------------
def clean_pbs():
    fopen = open(log_file,'a')
    print("Deleting pbs files ...", file=fopen)
    fopen.close()
    for file in sorted(glob.glob("*.pbs")):
        basename=file.split('.')[0]
        if os.path.isfile(file): os.remove(file)
        os.system('rm -f %s.o[0-9][0-9][0-9][0-9]*' %(basename))
        os.system('rm -f %s.e[0-9][0-9][0-9][0-9]*' %(basename))
#------------------------------------------------------------------------------------------
def send_pbs_files_to_queue(njobs, time_sleep):
    #print("Intentional latency time: 2 s ...")
    time.sleep(1.0)
    jobslist, pids, jobindex, inqueue = [], [], 0, 1
    for file in sorted(glob.glob("*.pbs")): jobslist.append(file)
    if len(jobslist) == 0:
        fopen = open(log_file,'a')
        print("There are not pbs files type or they already have been done", file=fopen)
        fopen.close()
        return 0 
    totaljobs=len(jobslist)
    timein=time.strftime("%c")
    fopen = open(log_file,'a')
    print("Total jobs found = %d" %(totaljobs), file=fopen)
    print("Enter to the pool at : Date %s" %(timein), file=fopen)
    fopen.close()
    inqueue=totaljobs
    while inqueue >= 1:
        time.sleep(float(time_sleep))
        workid=queued_processes() 
        work_queue=[]
        for elem in pids:
            if elem in workid: work_queue.append(elem)
        inqueue=len(work_queue)
        if (njobs > inqueue) and (jobindex < totaljobs): 
            #print("Intentional latency time: 1 s ...")
            time.sleep(0.5)
            inqueue=1 ## at least one in queue
            spearhead=jobslist[jobindex]
            jobid=qsub(spearhead, jobindex)
            pids.append(jobid)
            jobindex=jobindex+1
    clean_pbs()
    timeout=time.strftime("%c")
    fopen = open(log_file,'a')
    print("Out of the pool at : Date %s" %(timeout), file=fopen)
    fopen.close()
    return 0
#------------------------------------------------------------------------------------------
