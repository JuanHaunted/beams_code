from multiprocessing import Barrier
import numpy as np
import sympy as sp 
from functions import *

init_values = initial_values()
bar = init_values[2]
si_loc = init_values[0]
sf_loc = init_values[1]
discret = bar/(bar*10)
force_mat = process_forces(bar, discret)
reactions = calc_reactions_sup(force_mat, si_loc, sf_loc, discret)
print(reactions)


