import os.path
from inout_solids.getbilparam import get_a_int, get_a_str, get_a_float
from inout_solids.readbil import read_var_composition,  clustername
from inout_solids.messag import write_welcome
from utils_solids.libmoleculas import rename_molecule, sort_by_energy, sort_by_stoichiometry
from discriminate_solids.energy import cutter_energy
from vasp_solids.libperiodicos import readposcars, writeposcars
from gulp_solids.calculator_all import calculator_gulp_all_check
from discriminate_solids.fp_discriminator import discriminate_calculated
from utils_solids.miscellaneous import get_xcomp, uc_restriction, get_tolerances
from discriminate_solids.removal_by_descriptors import descriptor_comparison_calculated
#------------------------------------------------------------------------------------------------
# Variables
vol_restriction = uc_restriction() 
flag = get_a_str('calculator','vasp')
composition = read_var_composition('composition')
emax =  get_a_float('energy_range', 5.0)
atms_specie,atms_per_specie = get_xcomp(composition)
formula_units = get_a_int('formula_units',4)
#------------
z = len(atms_per_specie)
for i in range(z):
    atms_per_specie[i] = atms_per_specie[i] * formula_units
#------------
restart = get_a_str('restart','FALSE')
total_structures = get_a_int('initial_structures', 30)
dimension = get_a_int('dimension',3)
volume_factor = get_a_float('volume_factor', 1.0)
nofstages = get_a_int('number_of_stages', 1)
log_file = get_a_str('output_file','solids_out.txt')
l_tol,p_tol = get_tolerances(atms_specie)
simil_tol = get_a_float('similarity_tolerance',0.8)
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
def build_population_0(run):
    '''It builds the initial set of structures

    out:
    xtalist_out (list); Structures randomly generated using symmetry 
    '''
    initialfile = 'initial'+'_'+str(run)+'.vasp'
    if not os.path.isfile(initialfile):
        from utils_solids.randxtal import random_crystal_gen
        fopen = open(log_file,'a')
        print("Making initial file  = %s" %(initialfile),file=fopen)
        print("-------------------------------------------------------------------",file=fopen)
        print("-----------------------POPULATION  GENERATOR-----------------------",file=fopen)
        fopen.close()
        xtalist_out = random_crystal_gen(total_structures,atms_specie,atms_per_specie,p_tol,formula_units,dimension,volume_factor,vol_restriction)
        xtalist_out = rename_molecule(xtalist_out, 'random', 4)
        writeposcars(xtalist_out, initialfile, 'D')
    else:
        xtalist_out = readposcars(initialfile)
        if restart == 'FALSE':
            xtalist_out = rename_molecule(xtalist_out, 'restart', 3)
        fopen = open(log_file,'a')
        print("%s exists... we take it" %(initialfile), file=fopen)
        fopen.close()
    return xtalist_out

#------------------------------------------------------------------------------------------------
def run_calculator(poscarlistin, folder, stage=0):
    '''Receives a list of Molecule structures and send them to a relaxation code

    in:
    poscarlistin (list); Molecule list of all to-be-relaxed structures
    folder (str); Name of the folder where structures will be stored
    stage (int); Index number for the corresponding stage of the calculation

    out:
    poscarlistout (list); Molecule list of relaxed structures
    '''
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------",file=fopen)
    print("------------------- LOCAL OPTIMIZATION: STAGE %d -------------------" %(stage+1),file=fopen)
    fopen.close()
    if flag == 'vasp':
        from vasp_solids.calculator_all import calculator_vasp_all_check
        poscarlistin = sort_by_stoichiometry(poscarlistin)
        poscarlistout = calculator_vasp_all_check(poscarlistin,folder,0)
    if flag == 'gulp':
       from gulp_solids.calculator_all import calculator_gulp_all_check
       blockname = 'gulp'+str(stage+1)
       poscarlistout = calculator_gulp_all_check(poscarlistin, folder, blockname, stage)
       poscarlistout = sort_by_energy(poscarlistout,1)
    if flag=='mopac':
       blockg='mopac'+str(stage+1)
       from mopac.calculator_all import calculator_mopac_all_check
       poscarlistout=calculator_mopac_all_check(poscarlistin, folder, blockg, stage)
    poscarlistout = sort_by_energy(poscarlistout, 1)
    return poscarlistout

#------------------------------------------------------------------------------------------------
def display_mol_info(moleculein, stage=0, opt=0):
    fopen = open(log_file,'a')
    if opt == 0:
        print("-------------------------------------------------------------------",file=fopen)
        print("------------------------- SUMMARY STAGE %d -------------------------" %(stage+1), file=fopen)
        print("Number  File-Name   Energy (ev)    Delta-E  NT", file=fopen)
    fopen.close()
    molzz = sort_by_energy(moleculein,1)
    emin = molzz[0].e
    for ii, imol in enumerate(molzz):
        deltae  =  imol.e - emin
        nt = imol.c[0]
        fopen = open(log_file,'a')
        jj=str(ii+1).zfill(5)
        print("%5s %11s %14.8f %10.6f %d" %(jj, imol.i, imol.e, deltae, nt), file=fopen)
        fopen.close()

#------------------------------------------------------------------------------------------------
run = 1
for i in range(4):
    fopen = open(log_file,'a')
    print ("\n-------------------------------------------------------------", file=fopen)
    print ("START of RUN", file=fopen)
    print ("-------------------------------------------------------------", file=fopen)
    fopen.close()
    #poscar00 = build_population_0()
    poscar00 = build_population_0(run)    
    for stage in range(nofstages):
        basenm = 'stage'+str(stage+1)+'_'+str(run)
        folder = basenm+'/'
        poscar01 = rename_molecule(poscar00, basenm, 3)
        poscar00 = run_calculator(poscar01, folder, stage)
        poscar00 = cutter_energy(poscar00,emax,0)
        poscar00 = descriptor_comparison_calculated(poscar00,simil_tol)
        display_mol_info(poscar00,stage)
        writeposcars(poscar00, basenm + '.vasp', 'D')
    run = run + 1 
fopen = open(log_file,'a')
print ("-------------------------------------------------------------", file=fopen)
print ("SOLIDS HAS FINISHED SUCCESSFULLY", file=fopen)
print ("-------------------------------------------------------------", file=fopen)
fopen.close()
