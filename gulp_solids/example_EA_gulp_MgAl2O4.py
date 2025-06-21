from inout_solids.messag import *

qname = 'queue_name'

long_string = """
========= Optimization Scheme =========

option      EA

==== Crystal Structure Information ====

---COMPOSITION---
Mg   1
Al   2
O    4
---COMPOSITION---

formula_units   4

=========== Constrictions =============

tol_atomic_overlap 0.95

========= Initial Population ==========

initial_structures  30

======== Algorithm Parameters =========

max_number_inputs   5
number_of_matings   20
number_of_xchange   6
number_of_strains   6
number_of_randoms   3
max_number_gens     30
crit_stop_nrep      10

========== Discrimination =============

similarity_tolerance    0.95
energy_range            2.0

====== Calculation Parameters =========

calculator        gulp
qsys              local

---GULP---
opti conjugate nosymmetry conp
switch_minimiser bfgs gnorm 0.01
pressure 100 GPa
vectors
LATTICEVECTORS
frac
COORDINATES
space
1
species
Mg  2.0
Al  3.0
O  -2.0
lennard 12 6
Mg O   1.50 0.00 0.00 6.0
Al O   1.50 0.00 0.00 6.0
O O    1.50 0.00 0.00 6.0
Mg Mg  1.50 0.00 0.00 6.0
Mg Al  1.50 0.00 0.00 6.0
Al Al  1.50 0.00 0.00 6.0
buck
Mg O 1428.5 0.2945 0.0 0.0 7.0
Al O 1114.9 0.3118 0.0 0.0 7.0
O O  2023.8 0.2674 0.0 0.0 7.0
maxcyc 850
switch rfo 0.010
---GULP---

---GULP.CONF---
exe_gulp=GULP_Dir/gulp
---GULP.CONF---
""".format(queue=qname)

exfile = open('INPUT.txt', "w")
exfile.write(welcome_solids)
exfile.write(long_string)
exfile.close()
