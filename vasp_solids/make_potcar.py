import os
import os.path
from inout_solids.readbil import get_stchiom
from inout_solids.getbilparam  import get_a_str
log_file =get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
def get_potcar():
    elements=get_stchiom(1)
    pathpp = get_a_str('vasp_pp_path','ZERO')
    for ie in elements:
        target='POTCAR_'+ie
        if ( os.path.isfile(target) ):
            print('%s is found; take it' %(target))
        else:
            if (os.path.exists(pathpp)) or (pathpp=='ZERO'):
                origin=pathpp+'/'+ie+'/POTCAR'
                if ( os.path.isfile(origin) ):
                    task='cp '+origin+' '+target
                    fopen = open(log_file,'a')
                    print("Found: %s" %(origin), file=fopen)
                    fopen.close()
                    os.system(task)
                else:
                    print("WARNING!! %s do not exist!!" %(origin))
                    exit()
            else:
                print("WARNING!! %s do not exist!!" %(pathpp))
                print("Check the vasp_pp key in INPUT.txt")
                exit()
#------------------------------------------------------------------------------------------
def cat_potcar():
    elements=get_stchiom(1)
    filenames=[]
    for ie in elements:
        file='POTCAR_'+ie
        if os.path.exists(file):
            filenames.append(file)
        else:
            print("POTCAR_%s do not NOT found!!" %(ie))
            exit()
    fopen = open(log_file,'a')
    print("Creating: POTCAR file", file=fopen)
    fopen.close()
    with open('POTCAR', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
#------------------------------------------------------------------------------------------
def make_potcar_file():
    if ( os.path.isfile('./POTCAR') ):
        fopen = open(log_file,'a')
        print("POTCAR is found; take it", file=fopen)
        fopen.close()
    else:
        get_potcar()
        cat_potcar()
        elements=get_stchiom(1)
        for ie in elements: os.remove('POTCAR_'+ie)            
#------------------------------------------------------------------------------------------
