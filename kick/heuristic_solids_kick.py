import sys
sys.path.insert(0, '/home/carlos/installdir/solids/GLOMOSolids1.0/')

import os.path
from inout.getbilparam import get_a_int, get_a_str, get_a_float
from inout.readbil import read_var_composition,  clustername
from utils.libmoleculas import readxyzs, writexyzs, rename_molecule, sort_by_energy, sort_by_stoichiometry
from discriminate.energy import cutter_energy
from vasp.libperiodicos import readposcars, writeposcars, expand_poscar
from gulp.calculator_all import calculator_gulp_all_check
from discriminator import discriminate_calculated
from miscellaneous import get_xcomp, uc_restriction, get_tolerances
#------------------------------------------------------------------------------------------------
vol_restriction = uc_restriction() 
flag = get_a_str('calculator','vasp')
composition = read_var_composition('composition')
emax =  get_a_float('energy_range', 99.0)
atms_specie,atms_per_specie = get_xcomp(composition)
total_structures = get_a_int('initial_structures', 10)
formula_units = get_a_int('formula_units',4)
dimension = get_a_int('dimension',3)
volume_factor = get_a_float('volume_factor', 1.0)
nofstages = get_a_int('number_of_stages', 1)
log_file = get_a_str('output_file','glomos_out.txt')
l_tol,p_tol = get_tolerances(atms_specie)
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
    initialfile = 'initial.vasp'
    if not os.path.isfile(initialfile):
        from randxtal import random_crystal_gen
        fopen = open(log_file,'a')
        print("Making initial file  = %s" %(initialfile),file=fopen)
        print("-------------------------------------------------------------------",file=fopen)
        print("-----------------------POPULATION  GENERATOR-----------------------",file=fopen)
        fopen.close()
        poscarlist = random_crystal_gen(total_structures,atms_specie,atms_per_specie,p_tol,formula_units,dimension,volume_factor,vol_restriction)
        poscarlist = rename_molecule(poscarlist, 'random', 4)
        writeposcars(poscarlist, initialfile, 'D')
    else:
        poscarlist = readposcars(initialfile)
        poscarlist = rename_molecule(poscarlist, 'restart', 4)
        fopen = open(log_file,'a')
        print("%s exists... we take it" %(initialfile), file=fopen)
        fopen.close()
    return poscarlist

#------------------------------------------------------------------------------------------------
def run_calculator(poscarlistin, folder, stage=0):
    fopen = open(log_file,'a')
    print("-------------------------------------------------------------------",file=fopen)
    print("------------------- LOCAL OPTIMIZATION: STAGE %d -------------------" %(stage+1),file=fopen)
    fopen.close()
    if flag == 'vasp':
        from vasp.calculator_all import calculator_vasp_all_check
        poscarlistin = sort_by_stoichiometry(poscarlistin)
        poscarlistout = calculator_vasp_all_check(poscarlistin,folder,0)
    if flag == 'gulp':
       from gulp.calculator_all import calculator_gulp_all_check
       blockname = 'gulp'+str(stage+1)
       poscarlistout = calculator_gulp_all_check(poscarlistin, folder, blockname, stage)
       poscarlistout = sort_by_energy(poscarlistout,1)
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
poscar00 = build_population_0()
for stage in range(nofstages):
    basenm = 'stage'+str(stage+1)
    folder = basenm+'/'
    poscar01 = rename_molecule(poscar00, basenm, 3)
    poscar00 = run_calculator(poscar01, folder, stage)
    poscar00 = cutter_energy(poscar00,emax,1)
    poscar00 = discriminate_calculated(poscar00,vol_restriction)
    display_mol_info(poscar00,stage)
    writeposcars(poscar00, basenm + '.vasp', 'D')
fopen = open(log_file,'a')
print ("-------------------------------------------------------------", file=fopen)
print ("SOLIDS HAS FINISHED SUCCESSFULLY", file=fopen)
print ("-------------------------------------------------------------", file=fopen)
fopen.close()
