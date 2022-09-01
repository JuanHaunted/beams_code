import numpy as np
import sympy as sp 
from functions import *

init_values = initial_values()

supports = np.array(init_values[0],init_values[1])
bar = init_values[2]
discret = bar/(bar*10)

force_mat = process_forces(bar, discret)
