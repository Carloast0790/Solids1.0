from pyxtal import pyxtal
from pyxtal.interface.vasp import optimize
from ase.db import connect
from random import randint
from time import time
from vasp_solids.libperiodicos import readposcars, writeposcars
from utils_solids.miscellaneous import pyxtal2xyz

import warnings

# warnings.filterwarnings("ignore")

N = 10
elements = {"C": [6,8]}
for i in range(N):
    while True:
        sg = randint(2, 230)
        print('trying symmetry',sg)
        species = []
        numIons = []
        for ele in elements.keys():
            species.append(ele)
            if len(elements[ele]) == 2:
                num = randint(elements[ele][0], elements[ele][1])
                numIons.append(num)
            else:
                numIons.append(elements[ele])
        print(species,numIons)
        crystal = pyxtal()
        crystal.from_random(2, sg, ['C'], [6], force_pass=True)

        if crystal.valid:
            print('built structure')
            solidsxtal = pyxtal2xyz(crystal)
            writeposcars([solidsxtal],'succ.vasp','D')
            break
        else:
            print('failure\n')