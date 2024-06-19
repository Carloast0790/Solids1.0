import os
import time
import glob
from inout_solids.getbilparam import get_a_str
log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
def send_sh_files_to_local():
    jobslist=[]
    for file in sorted(glob.glob("*.sh")): jobslist.append(file)
    if len(jobslist) == 0:
        fopen = open(log_file,'a')
        print("There are not sh files type or they already have been done", file=fopen)
        fopen.close()
        return 0
    for xlocalsh in jobslist:
        timein=time.strftime("%c")
        fopen = open(log_file,'a')
        print("Enter to local submit: %s --> %s" %(xlocalsh, timein), file=fopen)
        fopen.close()
        os.system("bash "+xlocalsh)
        os.system("rm -f "+xlocalsh)
    return 0
#------------------------------------------------------------------------------------------
