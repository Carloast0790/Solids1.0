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

tol_atomic_overlap 0.9

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
opti conj conp
switch_minimiser bfgs gnorm 0.5
vectors
LATTICEVECTORS
frac
COORDINATES
species
Ti  2.196
O  -1.098
buck
Ti Ti 31120.1 0.1540 5.25  15
O  O  11782.7 0.2340 30.22 15
Ti O  16957.5 0.1940 12.59 15
lennard 12 6
Ti Ti   1   0 15
O  O    1   0 15
Ti O    1   0 15
---GULP---

---GULP.CONF---
exe_gulp=Dir/gulp
---GULP.CONF---
"""

exfile = open('INPUT.txt', "w")
exfile.write(welcome_solids)
exfile.write(long_string)
exfile.close()
