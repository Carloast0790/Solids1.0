from inout_solids.messag import *

long_string = """
========= Optimization scheme =========

option      GA

==== Crystal structure information ====

---COMPOSITION---
Mg   1
Si   1
O    3
---COMPOSITION---

formula_units   4

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
opti conjugate nosymmetry conv
switch_minimiser bfgs gnorm 0.01
vectors
LATTICEVECTORS
frac
COORDINATES
space
1
species
Mg 1.8
Si 2.4
O -1.4
lennard 12 6
Mg O  2.5 0.0 0.0 6.0
Mg Si 1.5 0.0 0.0 6.0
Si O  1.5 0.0 0.0 6.0
Mg Mg 1.5 0.0 0.0 6.0
Si O  1.5 0.0 0.0 6.0
O  O  2.5 0.0 0.0 6.0
buck
Mg O   806.915 0.291 2.346 0.0 10.0
Si O  1122.392 0.256 0.000 0.0 10.0
O O    792.329 0.362 31.58 0.0 10.0
Mg Mg  900.343 0.220 0.174 0.0 10.0
Mg Si 1536.282 0.185 0.000 0.0 10.0
Si Si 3516.558 0.150 0.000 0.0 10.0
maxcyc
300
switch rfo cycle 350
---GULP---

---GULP.CONF---
exe_gulp=GULP_Dir/gulp
---GULP.CONF---

"""

exfile = open('INPUT.txt', "w")
exfile.write(welcome_solids)
exfile.write(long_string)
exfile.close()