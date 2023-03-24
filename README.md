# GOMOSolids1.0
GLOMOS Solids is a tool written in Python used for crystal structure prediction.

## Updates:
The interface with VASP and GULP seems to be working correctly. Right now, the code is undergoing test for both interfaces, but seems nice.

There was a major change in the crossover scheme. We're keeping half of the main structure and filling the gaps in the chemical comp
by using the positions of the remaining half and a random selection of each atom's species. Also, there is an ongoing change in the 
minimum distance checking. We're about to start adding this change in distance checking to the crossover scheme.  

The min interatomic distances are now veryfied by two possibilities: The first is just a scalling factor, that will scale the addition
of the two covalent radius. The second one, will be provided by the user by pairs in the form 'S1 S2 min_dist'. This interatomic distances
have allowed the algorithm to find diamond and graphite structures of C8.

