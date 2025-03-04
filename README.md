# Solids1.0
Solids is a tool written in Python used for crystal structure prediction. It relies on two different algorithms to perform the energy ladscape explorations: A Stochastich Method (SM) designed for a quick exploration, and a Modified Genetic Algorithm (GA) for more complex and dedicated explorations. 
The first scheme builds random structures, relaxes them with one of the interfaced codes and presents the results in a POSCAR-type file. 
The second scheme expands the first one by applying operators of mutation and hereditary to the initial set of structures, allowing a more thorough exploration of the energy landscape.

Solids is interfaced with VASP and GULP to carry structural relaxations. 

# Dependencies 

Pyxtal
DScribe 

# Updates
- Correction of initial random generator
- Addition of GULP examples
- Addition of fresh random structures in GA
- Correction to discrimination and output formatting
