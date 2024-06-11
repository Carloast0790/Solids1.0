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
