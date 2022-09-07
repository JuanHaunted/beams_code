import scipy as sp
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import cumulative_trapezoid
from multiprocessing import Barrier
import numpy as np
from functions import *
from math import pi

#Only for rectangular profiles, will recerive circular too
def inertial_moment():
    r_thickness = float(input("Ingrese el espesor del perfil rectangular: "))
    r_height = float(input("Ingrese la altura del perfil rectangular: "))    
    i_ = (r_thickness*(r_height**3))/12

    return i_

def deflexion_slope(moments, E, I ):
    slopes = np.zeros(len(moments))
    moments_int = riemann_sum(discret, moments)
    for i in range(0, len(moments)):
        slopes[i] = (moments_int[i]/(E*I))
        i+=1
    return slopes

def deflexion(slopes, E, I):
    defx = np.zeros(len(slopes))
    slopes_int = riemann_sum(discret, slopes)
    for i in range(0, len(slopes)):
        defx[i] = (slopes_int[i]/(E*I))
        i+=1
    return defx

def max_sigma(moments, section):
    moment = np.zeros[2]
    moment[0] = abs(max(moments))
    moment[1] = abs(min(moments))
    max_moment = max(moment)
    sigma_max = (max_moment)/(section)
    return sigma_max
    

def max_tao(sheers, Q, I, t):
    sheer = np.zeros[2]
    sheer[0] = abs(max(sheers))
    sheer[1] = abs(min(sheers))
    max_sheer = max(sheer)
    tao_max = (max_sheer*Q)/(I*t)
    return tao_max

#Data processing
init_values = initial_values()
#Specifies bar conditions
bar_type = input("Ingrese el perfil de la barra (Circular (c), o Rectangular (r)): ")
my_S = float(input("Ingrese el modulo de sección de la viga: "))
my_E = float(input("Ingrese el modulo de elasticidad de la viga: "))
my_I = inertial_moment()
bar = init_values[2]
#Beam processing
si_loc = init_values[0]
sf_loc = init_values[1]
discret = bar/(bar*10)
force_mat = process_forces(bar, discret)
force_with_reactions = calc_reactions_sup(force_mat, si_loc, sf_loc, discret)
sheer_forces = calc_sheer_forces(force_with_reactions, bar, discret)
print(sheer_forces)
integral_list2 = riemann_sum(discret, sheer_forces[1,:])
#Aditional bar info, graphs of deflection
slopes_deflexion = deflexion_slope(integral_list2, my_E, my_I )
pure_deflexion = deflexion(slopes_deflexion, my_E, my_I )


#Graphing
plt.plot(sheer_forces[0,:], sheer_forces[1,:])
plt.xlabel('Puntos de la viga')
plt.ylabel('Cortantes')
plt.grid()
plt.show()
plt.plot(sheer_forces[0,:], integral_list2)
plt.xlabel('Puntos de la viga')
plt.ylabel('Momentos')
plt.grid()
plt.show()
plt.plot(sheer_forces[0,:], slopes_deflexion)
plt.xlabel('Puntos de la viga')
plt.ylabel('Momentos')
plt.grid()
plt.show()
plt.plot(sheer_forces[0,:], pure_deflexion)
plt.xlabel('Puntos de la viga')
plt.ylabel('Momentos')
plt.grid()
plt.show()








