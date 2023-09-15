import random
import numpy as np
from gega.aptitude import get_aptitude

def get_roulette_wheel_selection(moleculein, nmatings):
    fitness = list(get_aptitude(moleculein))
    sum_of_fitness = sum(fitness)
    nlist = list(range(len(moleculein)))
    previous_probability = 0.0
    pp = [previous_probability]
    for ii in nlist:
        previous_probability = previous_probability + (fitness[ii] / sum_of_fitness)
        pp.append(previous_probability)
    roulette = []
    for iix in range(nmatings):
        random_number = random.random()
        for ii in nlist:
            li, ls = pp[ii], pp[ii+1]
            if random_number >= li and random_number < ls:
                iselect = ii
                break
        roulette.extend([moleculein[iselect]])
    return roulette
