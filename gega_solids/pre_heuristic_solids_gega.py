import sys
sys.path.insert(0, '/home/carlos/installdir/solids/GLOMOSolids1.0/')

import os.path
import random
import numpy as np
from discriminate.energy import cutter_energy
from inout.getbilparam import get_a_int, get_a_str, get_a_float
from inout.readbil import read_var_composition, clustername
from utils.libmoleculas import rename_molecule, sort_by_energy, sort_by_stoichiometry
from vasp.libperiodicos import readposcars, writeposcars, expand_poscar
#------------------------------------------------------------------------------------------------
from gega.aptitude import get_aptitude
from discriminator import discriminate_calculated, discriminate_calculated_vs_pool
from crossoversolids import crossover, popgen_childs
from mutation_and_heredity import popgen_mutants, make_mutants
from miscellaneous import get_xcomp, uc_restriction, rescale_str,get_tolerances
#-------------------------------------------------------------------------------------------------
vol_restriction = uc_restriction()
flag = get_a_str('calculator','vasp')
composition = read_var_composition('composition')
atms_specie,atms_per_specie = get_xcomp(composition)
total_structures = get_a_int('initial_structures', 4)
formula_units = get_a_int('formula_units',2)
dimension = get_a_int('dimension',3)
volume_factor = get_a_float('volume_factor', 1.0)
log_file = get_a_str('output_file','glomos_out.txt')

l_tol,p_tol = get_tolerances(atms_specie)
#------------------------------------------------------------------------------------------------
max_number_inputs = get_a_int('max_number_inputs',20)
emax = get_a_float('energy_range', 3.0)
nmaxgen = get_a_int('max_number_gens',10)
nmaxrep = get_a_int('crit_stop_nrep',10)
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
        # poscarlist = random_crystal_gen(total_structures,atms_specie,atms_per_specie,formula_units,dimension,volume_factor,vol_restriction)
        poscarlist = random_crystal_gen(total_structures,atms_specie,atms_per_specie,p_tol,formula_units,dimension,volume_factor,vol_restriction)
        poscarlist = rename_molecule(poscarlist, 'random_000_', 4)
        writeposcars(poscarlist, initialfile, 'D')
    else:
        poscarlist = readposcars(initialfile)
        poscarlist = rename_molecule(poscarlist, 'restart_000_', 4)
        fopen = open(log_file,'a')
        print("%s exists... we take it" %(initialfile), file=fopen)
        fopen.close()
    return poscarlist

#------------------------------------------------------------------------------------------------
def run_calculator(poscarlistin, folder, generation=0):
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------", file=fopen)
    print("---------------- LOCAL OPTIMIZATION: GENERATION %d   ----------------" %(generation), file=fopen)
    fopen.close()
    if flag == 'vasp':
        from vasp.calculator_all import calculator_vasp_all_check
        poscarlistin = sort_by_stoichiometry(poscarlistin)
        poscarlistout = calculator_vasp_all_check(poscarlistin,folder,0)
    if flag == 'gulp':
        from gulp.calculator_all import calculator_gulp_all_check
        poscarlistin = sort_by_stoichiometry(poscarlistin)
        poscarlistout = calculator_gulp_all_check(poscarlistin,folder,'gulp',generation)
    return poscarlistout

#------------------------------------------------------------------------------------------------
def display_info(poscarlist, flagsum=0, generation=0):
    fopen = open(log_file,'a')
    if flagsum == 0:
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
        jj = str(ii+1).zfill(5)
        if flagsum == 0:
            print("#%s %16s  %12.8f %10.6f  %2d" %(jj, iposcar.i, iposcar.e,deltae, nt), file=fopen)
        else:
            fitness = list(get_aptitude(poscarlist))
            print("#%s %16s  %12.8f %10.6f  %5.2f" %(jj, iposcar.i, iposcar.e,deltae, fitness[ii]), file=fopen)
    fopen.close()

#------------------------------------------------------------------------------------------------
def build_population_n(poscarlist,ref_d,generation=1):
    initialfile = 'initial' + str(generation).zfill(3) + '.vasp'
    if os.path.isfile(initialfile):
        fopen = open(log_file,'a')
        print("\n------------------------- NEW  GENERATION -------------------------", file=fopen)
        print("%s exists... we take it" %(initialfile), file=fopen)
        fopen.close()
        xtal_out = readposcars(initialfile)
    else:
        xtal_out = []
        cross = popgen_childs(poscarlist, ref_d, generation)
        xtal_out.extend(cross)
        mut = popgen_mutants(poscarlist,generation)
        xtal_out.extend(mut)
        if vol_restriction:
            for x in xtal_out:
                x = rescale_str(x,vol_restriction)
        else:
            get_vol = lambda v0,v1,v2: abs(np.dot(np.cross(v0,v1),v2))
            best_mtx = poscarlist[0].m
            best_vol = get_vol(best_mtx[0],best_mtx[1],best_mtx[2])
            for x in xtal_out:
                x = rescale_str(x,best_vol)
        writeposcars(xtal_out, initialfile, 'D')
    return xtal_out

#------------------------------------------------------------------------------------------------
poscar00 = build_population_0()
poscar01 = run_calculator(poscar00, 'generation000/', 0)
poscar01 = discriminate_calculated(poscar01,vol_restriction)
poscar01 = sort_by_energy(poscar01,1)
poscarxx = poscar01[0:max_number_inputs]
# list_ref_distances = get_min_interatomic_distance(poscarxx[0],min_rad,scale_fact)
display_info(poscarxx, 0)
poscaryy = sort_by_energy(poscarxx, 0)
writeposcars(poscaryy,'relaxed_gen_'+str(000)+'.vasp','D')
writeposcars(poscaryy, 'summary.vasp', 'D')
emin = poscarxx[0].e
cont = 0
for generation in range(1,nmaxgen + 1):
    folder = 'generation' + str(generation).zfill(3) + '/'
    poscar00 = build_population_n(poscaryy,l_tol,generation)
    poscar01 = run_calculator(poscar00, folder, generation)
    poscar01 = cutter_energy(poscar01,emax,1)
    #discrimination among local structures
    poscar01 = discriminate_calculated(poscar01,vol_restriction)
    #discrimination btwn remaining local vs pool
    poscar02 = discriminate_calculated_vs_pool(poscar01,poscarxx,vol_restriction)
    if len(poscar02) == 0:
        fopen = open(log_file,'a')
        print("-------------------------------------------------------------------", file=fopen)
        print("-------------------------------------------------------------------", file=fopen)
        print('The Process has Stopped, We Were Unable to Build New Structures', file=fopen)
        fopen.close()
        exit()        
    poscar02 = sort_by_energy(poscar02,1)
    display_info(poscar02,1,generation)
    poscarxx.extend(poscar02)
    poscarxx = sort_by_energy(poscarxx, 1)
    poscarxx = poscarxx[0:max_number_inputs]
    display_info(poscarxx, 0)
    poscaryy = sort_by_energy(poscarxx, 0)
    writeposcars(poscaryy,'relaxed_gen_'+str(generation)+'.vasp','D')
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