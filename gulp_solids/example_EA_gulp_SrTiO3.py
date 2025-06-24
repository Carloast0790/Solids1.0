from inout_solids.messag import *

long_string = """
========= Optimization scheme =========

option      EA

==== Crystal structure information ====

---COMPOSITION---
Sr   1
Ti   1
O    3
---COMPOSITION---

formula_units   2

=========== Constrictions =============

tol_atomic_overlap 0.95

========= Initial population ==========

initial_structures  30

======== Algorithm Parameters =========

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
---GULP---

---GULP.CONF---
exe_gulp=GULP_Dir/gulp
---GULP.CONF---
"""

exfile = open('INPUT.txt', "w")
exfile.write(welcome_solids)
exfile.write(long_string)
exfile.close()
