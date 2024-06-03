from utils_solids.libmoleculas import sort_by_energy
from inout_solids.getbilparam  import get_a_str
log_file=get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
""" THIS ROUTINE RETURNS THE LIST OF MOLECULES THAT SATISFY THE
    CUT ENERGY CRITERION. 
    HOW TO USE:
    >>> from inout.xyz import readxyz, writexyz
    >>> from discriminate.energy import cutter_energy
    >>> molecule=readxyz('all_reoptim.xyz')
    >>> mol=cutter_energy(molecule,0.5)
    >>> writexyz(mol,'output.xyz')
"""
def cutter_energy(moleculein, enemax, silence=0):
    moleculeout, count=[],0
    if moleculein==[]: return moleculeout
    moleculesort=sort_by_energy(moleculein,1)
    emin0=moleculesort[0].e
    if silence==0:
        fopen = open(log_file,'a')
        print("\nmax_energy_allowed = %3.2f kcal/mol" %(enemax), file=fopen)
        fopen.close()
    for imol in moleculein:
        de=imol.e - emin0
        if ( de < float(enemax) ):
            moleculeout.extend([imol])
        else:
            count=count+1
            jj=str(count).zfill(5)
            if silence==0:
                fopen = open(log_file,'a')
                print("%s %15s ... DISCRIMINATED: DeltaE = %3.2f" %(jj, imol.i, de), file=fopen)
                fopen.close()
    if count==0 and silence==0:
        fopen = open(log_file,'a')
        print("ZERO elements discriminated by Energy", file=fopen)
        fopen.close()
    elif moleculeout==[] and silence==0:
        fopen = open(log_file,'a')
        print("All the elements were discriminated by Energy", file=fopen)
        print("Please choose another value for enemax", file=fopen)
        fopen.close()
        exit()
    return moleculeout
#------------------------------------------------------------------------------------------
