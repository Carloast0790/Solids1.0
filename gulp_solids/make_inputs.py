import os.path
from utils_solids.libmoleculas import molecular_stoichiometry
from inout_solids.getbilparam import get_a_str
from vasp_solids.libperiodicos import cartesian2direct
log_file =get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
def make_ainput_gulp(headgulp, singleposcar, folder):
    iposname=singleposcar.i
    nameinp=iposname+'.gin'
    nameout=iposname+'.got'
    if os.path.isfile(folder+nameout):
        fopen = open(log_file,'a')
        print ("%s ... was not built because %s exists" %(folder+nameinp, folder+nameout), file=fopen)
        fopen.close()
    else:
        fopen = open(log_file,'a')
        print("Making input file = %s" %(folder+nameinp), file=fopen)
        fopen.close()
        fh=open(folder+nameinp,"w")
        for iline in headgulp:
            if iline=='LATTICEVECTORS\n':
                matrix=singleposcar.m
                print("%12.9f %16.9f %16.9f" %(matrix[0,0],matrix[0,1],matrix[0,2]), file=fh)
                print("%12.9f %16.9f %16.9f" %(matrix[1,0],matrix[1,1],matrix[1,2]), file=fh)
                print("%12.9f %16.9f %16.9f" %(matrix[2,0],matrix[2,1],matrix[2,2]), file=fh)
            elif iline=='COORDINATES\n':
                matrix=singleposcar.m
                for iatom in singleposcar.atoms:
                    s,xc,yc,zc=iatom.s,iatom.xc,iatom.yc,iatom.zc
                    xd,yd,zd=cartesian2direct(xc,yc,zc, matrix)
                    print("%-2s %16.9f %16.9f %16.9f" %(s,xd,yd,zd), file=fh)
            else:
                fh.write(iline)
        fh.close()
#------------------------------------------------------------------------------------------
def run_sample():
    from inout_solids.readbil import read_block_of_bil
    from vasp_solids.libperiodicos import readposcars
    headgulp=read_block_of_bil('gulp')
    singleposcar=readposcars('ejemplos/anatasa.vasp')[0]
    singleposcar.i='anatasa'
    folder='ejemplos/'
    make_ainput_gulp(headgulp, singleposcar, folder)
#run_sample()
#------------------------------------------------------------------------------------------
