from inout_solids.messag import *

long_string = """
========= Optimization scheme =========

option      SM

==== Crystal structure information ====

---COMPOSITION---
Si   2
---COMPOSITION---

formula_units   4

=========== Constrictions =============

tol_atomic_overlap 0.9

========= Initial population ==========

initial_structures  50

======== Algorithm parameters =========

number_of_stages    2

========== Discrimination =============

similarity_tolerance    0.90
energy_range            5.0

====== Calculation parameters =========

calculator        gulp
qsys              local

---GULP1---
opti conj conp
switch_minimiser bfgs gnorm 0.01
vectors
LATTICEVECTORS
frac
COORDINATES
lib GULP_Dir/Libraries/meam_2nn.lib
maxcyc 950
---GULP1---

---GULP2---
opti conj conp
switch_minimiser bfgs gnorm 0.001
vectors
LATTICEVECTORS
frac
COORDINATES
lib GULP_Dir/Libraries/meam_2nn.lib
maxcyc 950
---GULP2---

---GULP.CONF---
exe_gulp=GULP_Dir/gulp
---GULP.CONF---
"""

exfile = open('INPUT.txt', "w")
exfile.write(welcome_solids)
exfile.write(long_string)
exfile.close()
