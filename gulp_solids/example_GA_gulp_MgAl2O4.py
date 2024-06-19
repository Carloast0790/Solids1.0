from inout_solids.messag import *

qname = 'queue_name'

long_string = """
========= Optimization Scheme =========

option      GA

==== Crystal Structure Information ====

---COMPOSITION---
Mg   1
Al   2
O    4
---COMPOSITION---

formula_units   4

=========== Constrictions =============

interatom_scale_value 0.9

========= Initial Population ==========

initial_structures  40

======== Algorithm Parameters =========

max_number_inputs   10
number_of_matings   60
number_of_mutants   20
number_of_randoms   20
max_number_gens     20
crit_stop_nrep      10

========== Discrimination =============

min_energy_difference   0.01
min_volume_difference   0.01
energy_range            10.0

====== Calculation Parameters =========

calculator        gulp
qsys              local

---GULP---
opti conjugate nosymmetry conv
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
