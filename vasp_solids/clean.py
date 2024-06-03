import os
import sys
sys.path.append('./')
#------------------------------------------------------------------------------------------
def clean_one_calc_vasp(path):
    for file in ("CHG","CHGCAR","DOSCAR","EIGENVAL","IBZKPT","PCDAT","PROCAR","REPORT","WAVECAR","XDATCAR","KPOINTS","POTCAR","OSZICAR","vasprun.xml"):
        pathfile=path+file
        if os.path.isfile(pathfile):
            #print "Remove:",pathfile
            os.remove(pathfile)
#------------------------------------------------------------------------------------------
def clean_vasp(path, foldername):
    clean_one_calc_vasp(path)
    b1="rm -f %s.e[0-9][0-9][0-9][0-9]*" %(foldername)
    b2="rm -f %s.o[0-9][0-9][0-9][0-9]*" %(foldername)
    os.system(b1)
    os.system(b2)
#------------------------------------------------------------------------------------------
