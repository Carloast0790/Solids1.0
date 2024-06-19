from inout_solids.messag import *

long_string = """
========= Optimization scheme =========

option      GA

==== Crystal structure information ====

---COMPOSITION---
Sr   1
Ti   1
O    3
---COMPOSITION---

formula_units   10

=========== Constrictions =============

interatom_scale_value 0.9

========= Initial population ==========

initial_structures  100

======== Algorithm parameters =========

max_number_inputs   40
number_of_matings   60
number_of_mutants   20
max_number_gens     20
crit_stop_nrep      10

========== Discrimination =============

min_energy_difference   0.01
min_volume_difference   0.01
energy_range            10.0

====== Calculation parameters =========

calculator        gulp
qsys              local

---GULP---
opti conjugate nosymmetry conp
switch_minimiser bfgs gnorm 0.01
vectors
LATTICEVECTORS
frac
COORDINATES
space
1
species
Sr 2.00
Ti 4.00
O -2.00
lennard 12 6
Sr Sr 1.0 0.0 0. 6.0
Sr Ti 1.0 0.0 0. 6.0
Sr O  2.0 0.0 0. 6.0
Ti Ti 1.0 0.0 0. 6.0
Ti O  2.0 0.0 0. 6.0
O  O  2.0 0.0 0. 6.0
buck
Sr Sr 9949.1  0.2446 0.0 0. 8.0
Sr Ti 12708.1 0.2191 0.0 0. 8.0
Sr O  1805.2  0.3250 0.0 0. 8.0
Ti Ti 16963.1 0.1847 0.0 0. 8.0
Ti O  845.0   0.3770 0.0 0. 8.0
O  O  22746.3 0.1490 0.0 0. 8.0
maxcyc 3850
switch bfgs gnorm 0.010
---GULP---

---GULP.CONF---
exe_gulp=GULP_Dir/gulp
---GULP.CONF---
"""

exfile = open('INPUT.txt', "w")
exfile.write(welcome_solids)
exfile.write(long_string)
exfile.close()
