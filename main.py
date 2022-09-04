import scipy as sp
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import cumulative_trapezoid
from multiprocessing import Barrier
import numpy as np
from functions import *

init_values = initial_values()
bar = init_values[2]
si_loc = init_values[0]
sf_loc = init_values[1]
discret = bar/(bar*10)
force_mat = process_forces(bar, discret)
force_with_reactions = calc_reactions_sup(force_mat, si_loc, sf_loc, discret)
sheer_forces = calc_sheer_forces(force_with_reactions, bar, discret)
integral_list = trapezoids_sum(discret, sheer_forces[1,:])
integral_list2 = riemann_sum(discret, sheer_forces[1,:])
print(sheer_forces)


#Graphing
plt.plot(sheer_forces[0,:], sheer_forces[1,:])
plt.xlabel('Variable independiente')
plt.ylabel('Variable dependiente')
plt.show()
plt.plot(sheer_forces[0,:], integral_list)
plt.xlabel('Variable independiente')
plt.ylabel('Variable dependiente acumulada')
plt.show()
plt.plot(sheer_forces[0,:], integral_list2)
plt.xlabel('Variable independiente')
plt.ylabel('Variable dependiente acumulada')
plt.show()

