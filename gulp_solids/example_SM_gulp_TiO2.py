from inout_solids.messag import *

long_string = """
========= Optimization scheme =========

option      SM

==== Crystal structure information ====

---COMPOSITION---
Ti   1
O    2
---COMPOSITION---

formula_units   4

=========== Constrictions =============

interatom_scale_value 0.8

========= Initial population ==========

initial_structures  50

======== Algorithm parameters =========

number_of_stages    2

========== Discrimination =============

min_energy_difference   0.01
min_volume_difference   0.001
energy_range            5.0

====== Calculation parameters =========

calculator        gulp
qsys              local

---GULP1---
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
---GULP1---

---GULP2---
opti conj conp
switch_minimiser bfgs gnorm 0.01
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
---GULP2---

---GULP.CONF---
exe_gulp=GULP_Dir/gulp
---GULP.CONF---

"""

exfile = open('INPUT.txt', "w")
exfile.write(welcome_solids)
exfile.write(long_string)
exfile.close()
