import sys
sys.path.insert(0, '/home/carlos/installdir/solids/glomos_solidsGH/')

import os.path
import random
from discriminate.energy import cutter_energy
from inout.getbilparam import get_a_int, get_a_str, get_a_float
from inout.readbil import read_var_composition, clustername
from utils.libmoleculas import readxyzs, writexyzs, rename_molecule, sort_by_energy
from vasp.libperiodicos import readposcars, writeposcars, expand_poscar
from gulp.calculator_all import calculator_gulp_all_check
#------------------------------------------------------------------------------------------------
from gega.aptitude import get_aptitude
from discriminator import discriminate_calculated, discriminate_calculated_vs_pool
from crossoversolids import crossover, popgen_childs
from mutation_and_heredity import popgen_mutants, make_mutants
#------------------------------------------------------------------------------------------------
def get_xcomp(composition):
    x = len(composition[0])
    species = []
    atms_per_specie = []
    for ii in range(x):
        s = composition[0][ii][0]
        species.append(s)
        n = composition[0][ii][1]
        atms_per_specie.append(n)
    return species, atms_per_specie
#------------------------------------------------------------------------------------------------
composition = read_var_composition('composition')
atms_specie,atms_per_specie = get_xcomp(composition)
total_structures = get_a_int('total_structures', 4)
formula_units = get_a_int('formula_units',1)
dimension = get_a_int('dimension',3)
volume_factor = get_a_float('volume_factor', 1.1)
log_file = get_a_str('output_file','glomos_out.txt')
#------------------------------------------------------------------------------------------------
max_number_inputs = get_a_int('max_number_inputs',20)
emax = get_a_float('energy_range', 100.0)
nmaxgen = get_a_int('max_number_gens',50)
nmaxrep=get_a_int('crit_stop_nrep',10)
#------------------------------------------------------------------------------------------------
pid = os.getpid()
fopen = open(log_file,'w')
print("-------------------------------------------------------------------", file=fopen)
print("Current PID          = %s" %(pid),file=fopen)
cf=clustername(composition[0])
print("Chemical Formula     = %s" %(cf), file=fopen)
print("Formula Units        = %s" %(formula_units), file=fopen)
fopen.close()
#------------------------------------------------------------------------------------------------
def build_population_0():
    initialfile = 'initial000.vasp'
    if not os.path.isfile(initialfile):
        from randxtal import random_crystal_gen
        fopen = open(log_file,'a')
        print("Making initial file = %s" %(initialfile),file=fopen)
        print("-------------------------------------------------------------------",file=fopen)
        print("-----------------------POPULATION  GENERATOR-----------------------",file=fopen)
        fopen.close()
        poscarlist=random_crystal_gen(total_structures,atms_specie,atms_per_specie,formula_units,dimension,volume_factor)
        poscarlist = rename_molecule(poscarlist, 'random_000_', 4)
        writeposcars(poscarlist, initialfile, 'D')
    else:
        poscarlist=readposcars(initialfile)
        fopen = open(log_file,'a')
        print("%s exist .... we take it" %(initialfile), file=fopen)
        fopen.close()
    return poscarlist

#------------------------------------------------------------------------------------------------
def run_calculator(poscarlistin, folder, generation=0):
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------", file=fopen)
    print("---------------- LOCAL OPTIMIZATION: GENERATION %d   ----------------" %(generation), file=fopen)
    fopen.close()
    blockname = 'gulp'
    poscarlistout = calculator_gulp_all_check(poscarlistin, folder, blockname, generation)
    return poscarlistout

#------------------------------------------------------------------------------------------------
def display_info(poscarlist, flagsum=0, generation=0):
    fopen = open(log_file,'a')
    if flagsum==0:
        print("\n--------------------------GLOBAL  SUMMARY--------------------------", file=fopen)
        print("energy_range = %f" %(emax), file=fopen)
        print("max_number_inputs = %d" %(max_number_inputs), file=fopen)
        print("Number    File-Name        Energy (eV)    Delta-E   NT", file=fopen)
    elif flagsum == 1:
        gen=str(generation).zfill(2)
        print("\n--------------------- GENERATION %s SUMMARY ----------------------" %(gen), file=fopen)
        print("Number    File-Name        Energy (eV)    Delta-E  Fitness", file=fopen)
    fopen.close()
    poscarzz = sort_by_energy(poscarlist,1)
    emin = poscarzz[0].e
    fopen = open(log_file,'a')
    for ii, iposcar in enumerate(poscarzz):
        deltae = iposcar.e - emin
        nt = iposcar.c[0]
        jj=str(ii+1).zfill(5)
        if flagsum==0:
            print("#%s %16s  %12.8f %10.6f  %2d" %(jj, iposcar.i, iposcar.e,deltae, nt), file=fopen)
        else:
            fitness = list(get_aptitude(poscarlist))
            print("#%s %16s  %12.8f %10.6f  %5.2f" %(jj, iposcar.i, iposcar.e,deltae, fitness[ii]), file=fopen)
    fopen.close()

#------------------------------------------------------------------------------------------------
def build_population_n(poscarlist, generation=1):
    initialfile = 'initial' + str(generation).zfill(3) + '.vasp'
    if os.path.isfile(initialfile):
        fopen = open(log_file,'a')
        print("\n------------------------- NEW  GENERATION -------------------------", file=fopen)
        print("%s exist .... we take it" %(initialfile), file=fopen)
        fopen.close()
        xtal_out=readposcars(initialfile)
    else:
        xtal_out = []
        cross = popgen_childs(poscarlist, generation)
        xtal_out.extend(cross)
        mut = popgen_mutants(poscarlist,generation)
        xtal_out.extend(mut)
        writeposcars(xtal_out, initialfile, 'D')
    return xtal_out

#------------------------------------------------------------------------------------------------
poscar00 = build_population_0()
generation =  0
folder = 'generation000/'
poscar01 = run_calculator(poscar00, folder, generation)
poscar01 = sort_by_energy(poscar01,1)
poscarxx = poscar01[0:max_number_inputs]
display_info(poscarxx, 0)
poscaryy = sort_by_energy(poscarxx, 0)
writeposcars(poscaryy, 'summary.vasp', 'D')
emin = poscarxx[0].e
cont = 0
for generation in range(1,nmaxgen + 1):
    folder = 'generation' + str(generation).zfill(3) + '/'
    poscar00 = build_population_n(poscarxx, generation)
    poscar01 = run_calculator(poscar00, folder, generation)
    poscar01 = sort_by_energy(poscar01,1)
    #discrimination among local structures
    poscar01 = discriminate_calculated(poscar01)
    #discrimination btwn remaining local vs pool
    poscar02 = discriminate_calculated_vs_pool(poscar01,poscarxx)
    if len(poscar02) == 0:
        fopen = open(log_file,'a')
        print("-------------------------------------------------------------------", file=fopen)
        print("-------------------------------------------------------------------", file=fopen)
        print('The process has stopped, we were unable to build new structures', file=fopen)
        fopen.close()
        exit()        
    display_info(poscar02, 1,generation)
    poscarxx.extend(poscar02)
    poscarxx = sort_by_energy(poscarxx, 1)
    poscarxx = poscarxx[0:max_number_inputs]
    poscarxx = cutter_energy(poscarxx, emax, 1)
    display_info(poscarxx, 0)
    poscaryy=sort_by_energy(poscarxx, 0)
    writeposcars(poscaryy, 'summary.vasp', 'D')
    emini = poscarxx[0].e
    if emini < emin:
        emin = emini
        cont = 0
    else:
        cont = cont + 1
        if cont >= nmaxrep:
            fopen = open(log_file,'a')
            print("-------------------------------------------------------------------", file=fopen)
            print("-------------------------------------------------------------------", file=fopen)
            print("STOP criterion satisfied (crit_stop_nrep = %s). Solids has finished" %(nmaxrep), file=fopen)
            fopen.close()
            exit()

fopen = open(log_file,'a')
print ("-------------------------------------------------------------", file=fopen)
print ("SOLIDS HAS FINISHED SUCCESSFULLY", file=fopen)
print ("-------------------------------------------------------------", file=fopen)
fopen.close()
