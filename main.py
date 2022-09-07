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


discret = bar/(bar*100)
p_moments = pure_moments(discret)
force_mat = process_forces(bar, discret)
calc_data = calc_reactions_sup(force_mat, si_loc, sf_loc, p_moments, discret)
force_with_reactions = calc_data[0]
moments = calc_data[1]
print(force_with_reactions)
print(moments)
sheer_forces = calc_sheer_forces(force_with_reactions, bar, discret)
integral_list = riemann_sum(discret, sheer_forces[1,:])
print(integral_list)

if len(moments) > 0:
    for i in range(len(moments[0])):
        idx = get_idx(moments[0, i], discret)
        integral_list[idx] += moments[1, i]
    
10
#Añadir los momentos puros al integral list (Momentos flexionates)
print(sheer_forces)


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

