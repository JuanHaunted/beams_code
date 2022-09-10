import scipy as sp
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import cumulative_trapezoid
from multiprocessing import Barrier
import numpy as np
from functions import *

init_values = initial_values()
if init_values[-1] == 0:
    bar = init_values[2]
    si_loc = init_values[0]
    sf_loc = init_values[1]
else:
    bar = init_values[1]
    si_loc = init_values[0]
    sf_loc = -1

aditional_values = aditional_info()
discret = bar/(bar*100)
p_moments = pure_moments(bar, discret)
force_mat = process_forces(bar, discret)[0]
calc_data = calc_reactions_sup(force_mat, si_loc, sf_loc, p_moments, discret)
force_with_reactions = calc_data[0]
moments = calc_data[1]
sheer_forces = calc_sheer_forces(force_with_reactions, bar, discret)
integral_list = riemann_sum(discret, sheer_forces[1,:])
if np.sum(moments[1]) != 0:
    integral_list = add_pure_moments(moments, integral_list, bar, discret)
if aditional_values[5] == "y":
    sigma_max = max_sigma(integral_list, aditional_values[0])
    tao_max = max_tao(sheer_forces[1], aditional_values[1], aditional_values[2], aditional_values[3], aditional_values[6], aditional_values[4])
elif aditional_values[5] == "n":
    sigma_max = 0
    tao_max = 0
    
#Graphing
plt.plot(np.concatenate([[0], sheer_forces[0,:]]), np.concatenate([[0], sheer_forces[1,:]]))
plt.xlabel('Puntos de la viga')
plt.ylabel('Cortantes')
plt.grid()
plt.show()
plt.plot(sheer_forces[0,:], integral_list)
plt.xlabel('Puntos de la viga')
plt.ylabel('Momentos')
plt.grid()
plt.show()

#Pruebas tao y sigma
print(sigma_max)
print(tao_max)

