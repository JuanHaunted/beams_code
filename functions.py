from operator import mod
from math import cos, sin
from tkinter import X
import numpy as np
import pandas as pd

#Se ingresan los inputs de tipo,posición de soportes,discretización ylongitud de barra

#Aprpxima números al diferencial m
def aprox_diff(num, diff_x):
    n_diff = num//diff_x
    residual = num % diff_x
    diff_mean = diff_x/2
    if residual >= diff_mean:
        return (n_diff*diff_x) + diff_x
    else:
        return n_diff * diff_x

def str_to_function(math_expr):
    if "x" in math_expr:
        return lambda x: eval(math_expr)
    else:
        return lambda x: (x*0) + (eval(math_expr) * 0.1)
    
def get_idx(xi, diff):
    return int(round(xi/diff))

def gen_beam_mat(b_len, x_diff):
    row_diff = np.arange(0 ,b_len + x_diff, x_diff)
    empty_row = np.zeros(len(row_diff))
    return np.stack([row_diff, empty_row])


def shift_mat(mat, s_value, diff):
    shift = get_idx(s_value, diff)
    new_ini_p = (-1) * diff * shift
    new_fnl_p = mat[-1] + new_ini_p
    shifted_mat = np.arange(new_ini_p, new_fnl_p + diff, diff)
    return shifted_mat

def moment_sum(mat, react_indx, mat_pure_moments = np.zeros(0)): #Editar para que funcione con empotradas
    #Rf = (sumatoria(momentos)+sumatoria(momentos_puros))/dist_rf
    dist_rf = mat[0, react_indx]
    neg_mat_pure_moments = mat_pure_moments * (-1)
    pm_sum = np.sum(neg_mat_pure_moments)
    m_sum = 0
    for i in range(0, len(mat[0])):
        m_sum += (mat[0, i]) * (mat[1, i]) 

    if dist_rf == 0:
        return m_sum + pm_sum
    else:
        return (m_sum + pm_sum)/dist_rf

def force_sum_y(force_arr):
    force_sum = 0
    for i in range(0, len(force_arr)):
        force_sum += force_arr[i]

    return force_sum

#Hacer funcion para empotradas!!!!
def calc_reactions_sup(beam_mat, sup_i, sup_f, pure_moments, diff):
    if sup_f == -1:
        sup_i_idx = 0
        print(beam_mat)
        react_y = force_sum_y(beam_mat[1])
        react_ang = moment_sum(beam_mat, sup_i_idx, pure_moments)
        pure_moments[1, 0] += react_ang  
        beam_mat[1, 0] -= react_y

        return [beam_mat, pure_moments]

    else:
        sup_i_idx = get_idx(sup_i, diff)
        sup_f_idx = get_idx(sup_f, diff)
        dist_mat = shift_mat(beam_mat[0], sup_i, diff)
        sh_beam_mat = np.stack((dist_mat, beam_mat[1]))
        print(sh_beam_mat)
        react_sf = moment_sum(sh_beam_mat, sup_f_idx)
        beam_mat[1, sup_f_idx] -= react_sf #corregir
        react_si = force_sum_y(beam_mat[1]) 
        beam_mat[1, sup_i_idx] -= react_si
        print(np.array([react_si, react_sf]))

        return [beam_mat, pure_moments]

def initial_values():
    bar_len = float(input("Ingrese la longitud de la barra: "))
    diff = bar_len/(bar_len*10)
    bar_type = int(input("Ingrese el tipo de barra, 0=con dos soportes o 1=barra empotrada: "))
    if bar_type == 0:

        sup_1 = aprox_diff(float(input("Donde quiere localizar su soporte 1: ")), diff)
        while sup_1<0 or sup_1>bar_len:
            print("Su localización de soporte está por fuerza de la barra")
            sup_1 = aprox_diff(float(input("Donde quiere localizar su soporte 1: ")), diff)

        sup_2 = aprox_diff(float(input("Donde quiere localizar su soporte 2: ")), diff)
        while sup_2<0 or sup_2>bar_len:
            print("Su localización de soporte está por fuerza de la barra")
            sup_2 = aprox_diff(float(input("Donde quiere localizar su soporte 2: ")), diff)

        return np.array([sup_1,sup_2,bar_len, bar_type]) 

    if bar_type == 1:
        print("El análisis se realizará automáticamente con la barra empotrada a la derecha")
        sup_1 = 0

        return np.array([sup_1, bar_len, bar_type])

def pure_moments(diff_x):
    flag = input("¿Desea ingresar momentos puros?(y-yes, n-no): ")
    if flag == 'n':
        return np.zeros((2, 1))
    n = int(input("¿Cuántos momentos puros desea ingresar?: "))
    moment_mat = np.zeros((2, n))
    while n:
        n -= 1
        loc = aprox_diff(float(input("Ingresa la localización del momento: ")), diff_x)
        value = float(input("Ingrese la magnitud del momento: "))
        moment_mat[0, n] += loc
        moment_mat[1, n] += value
    return moment_mat

def process_forces(bar_len, diff_x):
    f_mat = gen_beam_mat(bar_len, diff_x)

    
    punt_num = int(input("Ingrese la cantidad de cargas puntuales: "))
    while punt_num:
        punt_num -= 1   
        force_value = float(input("Ingrese los Newtons de dicha fuerza puntual: "))#poner catch
        force_loc = aprox_diff(float(input("Ingrese en donde quiere localizar su fuerza puntual: ")), diff_x)
        index = get_idx(force_loc, diff_x)
        print(index)
        f_mat[1, index] += np.array([force_value]) 
    

    distr_num = int(input("¿Cuántas cargas distribuidas desea poner?: "))
    while distr_num:
        distr_num -= 1  

        #Display de funciones
        math_expr = str_to_function(str(input("Ingrese la función de su carga distribuida, con la debida notación de Python: ")))
        begin_fun = aprox_diff(float(input("Ingrese la posicion inicial de la carga distribuida: ")), diff_x)
        end_fun = aprox_diff(float(input("Ingrese la posición final de la carga distribuida: ")), diff_x)
        dom_f = aprox_diff(float(input("Ingrese el valor inicial del dominio de la función: ")), diff_x)
        dist_diff = (end_fun-begin_fun)
        row_diff = np.arange(dom_f ,dom_f + dist_diff + diff_x, diff_x)
        y_images = math_expr(row_diff)
        ini_idx = get_idx(begin_fun, diff_x)


        j=0
        for i in range(ini_idx, len(y_images)+ini_idx):
            f_mat[1, i] += y_images[j]
            j+=1

        print(y_images)
        print(f_mat)
    return f_mat

def calc_sheer_forces(mat, bar_len, diff):
    sheer_mat = gen_beam_mat(bar_len, diff)
    acum_integral = 0
    for i in range(0, len(mat[0])):
        acum_integral += mat[1, i]
        sheer_mat[1, i] = acum_integral

    sheer_mat[1] *= (-1)
    return sheer_mat

def riemann_sum(diff_x, y_images):
    integral = np.zeros(len(y_images))
    Area = 0
    for i in range(0, len(y_images)):
        Area += (diff_x)*(y_images[i])
        integral[i] = Area
    i += 1
    return integral





#force_value=float(input("Ingrese los Newtons de dicha fuerza puntual: "))

#arr = np.array([[20, 43, 32, 32, 31],[2,1,4,2,3]])
#arr[1][2] += 3
#print(get_idx(2,))



def main():
    pass



