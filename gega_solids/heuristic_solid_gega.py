import os.path
import numpy as np
from discriminate_solids.energy import cutter_energy
from discriminate_solids.fp_discriminator import discriminate_calculated, discriminate_calculated_vs_pool
from inout_solids.getbilparam import get_a_int, get_a_str, get_a_float
from inout_solids.readbil import read_var_composition, clustername
from inout_solids.messag import write_welcome
from utils_solids.libmoleculas import rename_molecule, sort_by_energy, sort_by_stoichiometry
from utils_solids.miscellaneous import get_xcomp, uc_restriction, rescale_str,get_tolerances
from utils_solids.fresh_random_gen import popgen_fresh_random_gen
from gega_solids.aptitude import get_aptitude
from gega_solids.crossoversolids import crossover, popgen_childs
from gega_solids.mutation_and_heredity import popgen_mutants, make_mutants
from vasp_solids.libperiodicos import readposcars, writeposcars
#------------------------------------------------------------------------------------------------
# Variables
vol_restriction = uc_restriction()
flag = get_a_str('calculator','vasp')
composition = read_var_composition('composition')
atms_specie,atms_per_specie = get_xcomp(composition)
formula_units = get_a_int('formula_units',2)
#------------------
z = len(atms_per_specie)
for i in range(z):
    atms_per_specie[i] = atms_per_specie[i] * formula_units
#-----------------
total_structures   = get_a_int('initial_structures', 4)
max_number_inputs  = get_a_int('max_number_inputs',20)
nummber_of_randoms = get_a_int('number_of_randoms',10)
volume_factor = get_a_float('volume_factor', 1.0)
l_tol,p_tol   = get_tolerances(atms_specie)
dimension = get_a_int('dimension',3)
log_file  = get_a_str('output_file','solids_out.txt')
nmaxgen   = get_a_int('max_number_gens',10)
nmaxrep   = get_a_int('crit_stop_nrep',10)
emax = get_a_float('energy_range', 3.0)
#------------------------------------------------------------------------------------------------
pid = os.getpid()
fopen = open(log_file,'w')
write_welcome(fopen)
print("-------------------------------------------------------------------", file=fopen)
print("Current PID          = %s" %(pid),file=fopen)
cf=clustername(composition[0])
print("Chemical Formula     = %s" %(cf), file=fopen)
print("Formula Units        = %s" %(formula_units), file=fopen)
fopen.close()
#------------------------------------------------------------------------------------------------
def build_population_0():
    '''Creates the initial population of crystalline structures.

    out: 
    poscarlist (list); List with Molecule objects  
    '''
    initialfile = 'initial000.vasp'
    if not os.path.isfile(initialfile):
        from utils_solids.randxtal import random_crystal_gen
        fopen = open(log_file,'a')
        print("Making initial file  = %s" %(initialfile),file=fopen)
        print("-------------------------------------------------------------------",file=fopen)
        print("-----------------------POPULATION  GENERATOR-----------------------",file=fopen)
        fopen.close()
        poscarlist = random_crystal_gen(total_structures,atms_specie,atms_per_specie,p_tol,formula_units,dimension,volume_factor,vol_restriction)
        poscarlist = rename_molecule(poscarlist, 'random_000_', 3)
        writeposcars(poscarlist, initialfile, 'D')
    else:
        poscarlist = readposcars(initialfile)
        poscarlist = rename_molecule(poscarlist, 'restart_000_', 3)
        fopen = open(log_file,'a')
        print("%s exists... we take it" %(initialfile), file=fopen)
        fopen.close()
    return poscarlist

#------------------------------------------------------------------------------------------------
def run_calculator(poscarlistin, folder, generation=0):
    '''Send every crystal structure to the desired relaxation software

    in: 
    poscarlistin (list); Every to-be-relaxed structure
    folder (str); The name of the folder where POSCAR files will be stored
    generation (int); The index to be used for the  
    '''
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------", file=fopen)
    print("---------------- LOCAL OPTIMIZATION: GENERATION %d   ----------------" %(generation), file=fopen)
    fopen.close()
    if flag == 'vasp':
        from vasp_solids.calculator_all import calculator_vasp_all_check
        poscarlistin = sort_by_stoichiometry(poscarlistin)
        poscarlistout = calculator_vasp_all_check(poscarlistin,folder,0)
    if flag == 'gulp':
        from gulp_solids.calculator_all import calculator_gulp_all_check
        poscarlistin = sort_by_stoichiometry(poscarlistin)
        poscarlistout = calculator_gulp_all_check(poscarlistin,folder,'gulp',generation)
    poscarlistout = sort_by_energy(poscarlistout, 1)
    return poscarlistout

#------------------------------------------------------------------------------------------------
def display_info(poscarlist, flagsum=0, generation=0):
    disp_emax = round(emax,2)
    fopen = open(log_file,'a')
    if flagsum == 0:
        print("\n--------------------------GLOBAL  SUMMARY--------------------------", file=fopen)
        print("energy_range = %f eV" %(disp_emax), file=fopen)
        print("max_number_inputs = %d" %(max_number_inputs), file=fopen)
        print("Number    File-Name        Energy (eV)    Delta-E   NT", file=fopen)
    elif flagsum == 1:
        gen=str(generation).zfill(3)
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
            print("#%s %16s  %12.8f %10.6f  %2d" %(jj, iposcar.i, iposcar.e, deltae, nt), file=fopen)
        else:
            fitness = list(get_aptitude(poscarlist))
            print("#%s %16s  %12.8f %10.6f  %5.2f" %(jj, iposcar.i, iposcar.e, deltae, fitness[ii]), file=fopen)
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
        new_strs = popgen_fresh_random_gen(nummber_of_randoms,atms_specie,atms_per_specie,generation,formula_units,dimension)
        xtal_out.extend(new_strs)
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
s_initial = build_population_0()
s_relax = run_calculator(s_initial, 'generation000/', 0)
s_clean = discriminate_calculated(s_relax,vol_restriction)
s_ready = s_clean[:max_number_inputs]
display_info(s_ready, 0)
writeposcars(s_ready, 'summary.vasp', 'D')
emin = s_ready[0].e
cont = 0
for generation in range(1,nmaxgen + 1):
    folder = 'generation' + str(generation).zfill(3) + '/'
    s_initial_n = build_population_n(s_ready,l_tol,generation)
    s_relax_n = run_calculator(s_initial_n, folder, generation)
    s_clean_gen_n = discriminate_calculated(s_relax_n,vol_restriction)
    s_cleanpool_n = discriminate_calculated_vs_pool(s_clean_gen_n,s_ready,vol_restriction)
    if len(s_cleanpool_n) == 0:
        fopen = open(log_file,'a')
        print("-------------------------------------------------------------------", file=fopen)
        print("-------------------------------------------------------------------", file=fopen)
        print('The Process has Stopped, Solids was Unable to Build New Structures', file=fopen)
        fopen.close()
        exit()
    s_ready_n = sort_by_energy(s_cleanpool_n,1)
    s_ready_n = cutter_energy(s_ready_n,emax)
    s_ready_n = s_ready_n[0:max_number_inputs]
    display_info(s_ready_n,1,generation)
    s_ready.extend(s_ready_n)
    s_ready = sort_by_energy(s_ready, 1)
    display_info(s_ready, 0)
    writeposcars(s_ready, 'summary.vasp', 'D')
    emin_n = s_ready[0].e
    if emin_n < emin:
        emin = emin_n
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
