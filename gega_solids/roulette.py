import random
import numpy as np
from gega_solids.aptitude import get_aptitude
#from inout.getbilparam import get_a_str
#log_file=get_a_str('output_file','glomos_out.txt')
#------------------------------------------------------------------------------------------
def get_roulette_wheel_selection_beta(moleculein, nmatings):
    fitness=list(get_aptitude(moleculein))
    nlist=list(range(len(moleculein)))
    roulette = []
    for iix in range(nmatings):
        tryv=0
        while tryv == 0:
            ni=random.choice(nlist)
            rr01 = random.random()
            if(fitness[ni] > rr01):
                moleculein[ni].c=[fitness[ni]]
                roulette.extend([moleculein[ni]])
                tryv=1
    return roulette
#------------------------------------------------------------------------------------------
def get_roulette_wheel_selection(moleculein, nmatings):
    fitness=list(get_aptitude(moleculein))
    sum_of_fitness=sum(fitness)
    nlist=list(range(len(moleculein)))
    previous_probability = 0.0
    pp=[previous_probability]
    for ii in nlist:
        previous_probability=previous_probability + (fitness[ii] / sum_of_fitness)
        pp.append(previous_probability)
    #logfile = open(log_file,'a')
    #print("fitness list:", file=logfile)
    #print(np.around(fitness,decimals=2), file=logfile)
    #print("probab range:", file=logfile)
    #print(np.around(pp,decimals=2), file=logfile)
    #logfile.close()
    roulette = []
    for iix in range(nmatings):
        random_number = random.random()
        for ii in nlist:
            li, ls=pp[ii], pp[ii+1]
            if (random_number >= li) and (random_number < ls):
                iselect, lis, lss = ii, li, ls
                break
        #logfile = open(log_file,'a')
        #print("random_number=%3.2f --> index=%2d (%s)(f=%3.2f)" %(random_number, iselect+1, moleculein[iselect].i, fitness[iselect]), file=logfile)
        #logfile.close()
        moleculein[iselect].c=[fitness[iselect], iselect, random_number, lis, lss]
        roulette.extend([moleculein[iselect]])
    return roulette
#------------------------------------------------------------------------------------------
mole="""6
-549142.86818553     stage200001
Ag     -0.984025000      1.268430000      0.000441000
Ag     -0.606747000     -1.486261000      0.000346000
Ag      1.223362000      2.995596000     -0.000035000
Ag      1.590490000      0.217712000     -0.000373000
Ag      1.983485000     -2.556773000      0.000012000
Ag     -3.206566000     -0.438705000     -0.000392000
6
-549141.46446082     stage200002
Ag      0.000194000     -0.495234000     -1.449945000
Ag     -2.300654000     -1.125987000      0.000155000
Ag     -0.000144000     -0.493479000      1.450760000
Ag      2.300729000     -1.125936000      0.000144000
Ag      1.366330000      1.620389000     -0.000368000
Ag     -1.366455000      1.620246000     -0.000746000
6
-549140.48342833     stage200003
Ag      0.000000000     -0.552909000      0.000316000
Ag     -0.000001000      2.391552000     -0.000494000
Ag     -2.351733000      0.957631000      0.000301000
Ag      2.351732000      0.957632000      0.000296000
Ag     -2.439119000     -1.876954000     -0.000212000
Ag      2.439121000     -1.876952000     -0.000207000
6
-549139.97362840     stage200004
Ag     -0.656377000     -1.506650000     -0.941700000
Ag     -0.664533000      1.509773000     -0.934828000
Ag     -2.566349000     -0.005774000      0.471227000
Ag      0.039251000     -0.000884000      1.450489000
Ag      1.926466000     -1.465371000     -0.021452000
Ag      1.921542000      1.468907000     -0.023735000
"""
def run_sample():
    exfile = open("mol.xyz", "w")
    exfile.write(mole)
    exfile.close()
    from utils.libmoleculas import readxyzs, writexyzs
    mol=readxyzs('mol.xyz')
    nchilds=10
    mama = get_roulette_wheel_selection(mol, nchilds)
    for imol in mama: print(imol.c)
    writexyzs(mama,'mama.xyz')
#run_sample()
#------------------------------------------------------------------------------------------
##https://en.wikipedia.org/wiki/Fitness_proportionate_selection
##imprimir fitness, probabilidad, ruleta
