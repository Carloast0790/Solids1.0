from utils_solids.libmoleculas import sort_by_energy
from inout_solids.getbilparam  import get_a_str

log_file = get_a_str('output_file','solids_out.txt')
#------------------------------------------------------------------------------------------
def cutter_energy(xtalist_in, enemax, silence=0):
    """Returns a list of structures whose energies lies under the cut energy criteria

    in:
    xtalist_in (List); List of Molecule objects that will undergoe the cut

    out:
    xtalis_out (List); Filtered list with only structures under energetic cutoff
    """
    xtalist_out, count = [],0
    if xtalist_in == []: 
        return xtalist_out
    moleculesort = sort_by_energy(xtalist_in,1)
    emin0 = moleculesort[0].e
    fopen = open(log_file,'a')
    print('\n-------------------------------------------------------------------',file=fopen)
    print('------------------ Structure Removal by Energy --------------------',file=fopen)
    if silence == 0:
        print("\nMaximum energy gap allowed = %3.2f eV \n" %(enemax), file=fopen)
    for imol in xtalist_in:
        de = imol.e - emin0
        if de < enemax:
            xtalist_out.extend([imol])
        else:
            count = count+1
            jj = str(count).zfill(5)
            if silence == 0:
                print("%15s Removed: DeltaE = %3.2f" %(imol.i, de), file=fopen)
    if count == 0 and silence == 0:
        print("ZERO elements Removed by Energy", file=fopen)
    elif xtalist_out == [] and silence == 0:
        print("All the elements were Removed by Energy", file=fopen)
        print("Please choose another value for energy_range", file=fopen)
        fopen.close()
        exit()
    fopen.close()
    return xtalist_out
#------------------------------------------------------------------------------------------
