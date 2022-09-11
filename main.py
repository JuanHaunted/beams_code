import scipy as sp
#import simpy as sp
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
func = process_forces(bar, discret)
force_mat = func[0]
calc_data = calc_reactions_sup(force_mat, si_loc, sf_loc, p_moments, discret)
force_with_reactions = calc_data[0]
moments = calc_data[1]
sheer_forces = calc_sheer_forces(force_with_reactions, bar, discret)
integral_list = riemann_sum(discret, sheer_forces[1,:])
if np.sum(moments[1]) != 0:
    integral_list = add_pure_moments(moments, integral_list, bar, discret)
if aditional_values[5] == "y":
    sigma_max = max_sigma(integral_list, float(aditional_values[0]))
    tao_max = max_tao(sheer_forces[1], float(aditional_values[1]), float(aditional_values[2]), float(aditional_values[3]), aditional_values[6], float(aditional_values[4]))
elif aditional_values[5] == "n":
    sigma_max = 0
    tao_max = 0

#Graphing, reemplazar por el viejito
bar_typ = init_values[3]
locs_sup = [si_loc, -3]
loc_der_sup, loc_iz_sup = find_point(locs_sup)
locs_sup2 = [sf_loc, -3]
loc_der_sup2, loc_iz_sup2 = find_point(locs_sup2)
force_m = func[1]
distri = func[2]

gen_graph(bar, bar_typ, locs_sup, loc_der_sup, loc_iz_sup, locs_sup2,loc_der_sup2, loc_iz_sup2,
          force_m, distri, p_moments)

cool_graphs(sheer_forces,moments,discret,bar)

#Pruebas tao y sigma
print(sigma_max)
print(tao_max)
