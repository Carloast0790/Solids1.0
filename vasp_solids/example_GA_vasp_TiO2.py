from inout_solids.messag import *

long_string = """
========= Optimization scheme =========

option      GA

==== Crystal structure information ====

---COMPOSITION---
Ti   1
O    2
---COMPOSITION---

formula_units   4

=========== Constrictions =============

tol_atomic_overlap 0.95

========= Initial population ==========

initial_structures  20

======== Algorithm parameters =========

max_number_inputs   10
number_of_matings   20
number_of_xchange   6
number_of_strains   6
number_of_randoms   3
max_number_gens     30
crit_stop_nrep      10

========== Discrimination =============

similarity_tolerance    0.95
energy_range            2.0

====== Calculation parameters =========

no_attempts_opt         2
percent_of_convergence  85.0

## Queue system and Computer resources
queue                   {queue}
njobs                   10
nprocshared             10
memory_in_gb            8
walltime                08:00:00
timesleep               15.0

## Periodic/VASP parameters
calculator              vasp
vasp_pp_path            Directory/PBE/potpawPBE54
incar_files             incar_1
kpoints_files           kpoints

---VASP.CONF---
env > entorno-$PBS_JOBID.txt
export PATH=/opt/intel/impi/4.1.1.036/intel64/bin/:$PATH
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64/:/opt/intel/lib/intel64/:$PATH
exe_mpirun=/opt/intel/impi/4.1.1.036/intel64/bin/mpirun
exe_vasp=/LUSTRE/software/intel/vasp/vasp_std-intel

RUNLINE

rm -rf $PBS_JOBNAME-$PBS_JOBID
rm entorno-$PBS_JOBID.txt
---VASP.CONF---
"""

exfile = open('INPUT.txt', "w")
exfile.write(welcome_solids)
exfile.write(long_string)
exfile.close()
