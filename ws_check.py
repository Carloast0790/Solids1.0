import numpy as np
from copy import deepcopy
from pyxtal.symmetry import Group
#------------------------------------------------------------------------------------------------
def check_compatible(group, numIons):
     """
     Checks if the number of atoms is compatible with the Wyckoff
     positions. Considers the number of degrees of freedom for each Wyckoff
     position, and makes sure at least one valid combination of WP's exists.
     NOTE Comprhys: Is degrees of freedom used symnomously with multiplicity?
     perhaps standardising to multiplicity would be clearer?
     """
     # Store whether or not at least one degree of freedom exists
     has_freedom = False
     # Store the wp's already used that don't have any freedom
     used_indices = []
     # Loop over species
     for numIon in numIons:
         # Get lists of multiplicity, maxn and freedom
         l_mult0 = []
         l_maxn0 = []
         l_free0 = []
         indices0 = []
         for i_wp, wp in enumerate(group):
             indices0.append(i_wp)
             l_mult0.append(len(wp))
             l_maxn0.append(numIon // len(wp))
             if np.allclose(wp[0].rotation_matrix, np.zeros([3, 3])):
                 l_free0.append(False)
             else:
                 l_free0.append(True)
         # Remove redundant multiplicities:
         l_mult = []
         l_maxn = []
         l_free = []
         indices = []
         for mult, maxn, free, i_wp in zip(l_mult0, l_maxn0, l_free0, indices0):
             if free is True:
                 if mult not in l_mult:
                     l_mult.append(mult)
                     l_maxn.append(maxn)
                     l_free.append(True)
                     indices.append(i_wp)
             elif free is False and i_wp not in used_indices:
                 l_mult.append(mult)
                 indices.append(i_wp)
                 if mult <= numIon:
                     l_maxn.append(1)
                 elif mult > numIon:
                     l_maxn.append(0)
                 l_free.append(False)
         # Loop over possible combinations
         p = 0  # Create pointer variable to move through lists

         # Store the number of each WP, used across possible WP combinations
         n0 = [0] * len(l_mult)
         n = deepcopy(n0)
         for i, mult in enumerate(l_mult):
             if l_maxn[i] != 0:
                 p = i
                 n[i] = l_maxn[i]
                 break
         p2 = p
         if n == n0:
             return False
         while True:
             num = np.dot(n, l_mult)
             dobackwards = False
             # The combination works: move to next species
             if num == numIon:
                 # Check if at least one degree of freedom exists
                 for val, free, i_wp in zip(n, l_free, indices):
                     if val > 0:
                         if free is True:
                             has_freedom = True
                         elif free is False:
                             used_indices.append(i_wp)
                 break
             # All combinations failed: return False
             if n == n0 and p >= len(l_mult) - 1:
                 return False
             # Too few atoms
             if num < numIon:
                 # Forwards routine
                 # Move p to the right and max out
                 if p < len(l_mult) - 1:
                     p += 1
                     n[p] = min((numIon - num) // l_mult[p], l_maxn[p])
                 elif p == len(l_mult) - 1:
                     # p is already at last position: trigger backwards routine
                     dobackwards = True
             # Too many atoms
             if num > numIon or dobackwards is True:
                 # Backwards routine
                 # Set n[p] to 0, move p backwards to non-zero, and decrease by 1
                 n[p] = 0
                 while p > 0 and p > p2:
                     p -= 1
                     if n[p] != 0:
                         n[p] -= 1
                         if n[p] == 0 and p == p2:
                             p2 = p + 1
                         break
     if has_freedom:
         # All species passed: return True
         return True
     else:
         # All species passed, but no degrees of freedom: return 0
         return 0

#------------------------------------------------------------------------------------------------
def possible_sym(total_atms):
    """
    This function returns a list of the possible symmetry groups for a given stoichiometry

    in: total_atms (int); the total amount of atoms in the stoicheometry

    out: pos_sym (list); a list containing the possible simetry for the crystal in its number form
    """
    pos_sym = []
    for ii in range(1,230):
        sg = Group (ii)
        sg_symbol = str(sg.symbol)
        x = check_compatible(sg,total_atms)
        if x == True:
           pos_sym.append(sg_symbol)
    return pos_sym

#------------------------------------------------------------------------------------------------
