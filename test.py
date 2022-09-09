import numpy as np
from functions import *

b = 20
dis = 0.1

def process_forces_old(bar_len, diff_x):
    f_mat = gen_beam_mat(bar_len, diff_x)

    
    punt_num = int(input("Ingrese la cantidad de cargas puntuales: "))
    i = 1
    while punt_num:
        punt_num -= 1   
        force_loc = aprox_diff(float(input(f"Ingrese en donde quiere localizar su fuerza puntual {i}: ")), diff_x)
        force_value = float(input(f"Ingrese la fuerza puntual {i} en Newtons: "))#poner catch
        index = get_idx(force_loc, diff_x) 
        f_mat[1, index] += np.array([force_value]) 
        i += 1

    distr_num = int(input("¿Cuántas cargas distribuidas desea poner?: "))
    k = 1
    while distr_num:
        distr_num -= 1  
        #Display de funciones
        math_expr = str_to_function(str(input(f"Ingrese la función la carga distribuida {k} (en N/m) , con la debida notación de Python: ")))
        begin_fun = aprox_diff(float(input(f"Ingrese la posicion inicial de la carga distribuida {k}: ")), diff_x)
        end_fun = aprox_diff(float(input(f"Ingrese la posición final de la carga distribuida {k}: ")), diff_x)
        dom_f = aprox_diff(float(input(f"Ingrese el valor inicial del dominio de la función {k}: ")), diff_x)
        dist_diff = end_fun - begin_fun
        row_diff = np.arange(dom_f ,dom_f + dist_diff + diff_x, diff_x)
        y_images = math_expr(row_diff)
        ini_idx = get_idx(begin_fun, diff_x)
        k += 1


        j=0
        for i in range(ini_idx, len(y_images)+ini_idx):
            f_mat[1, i] += y_images[j]
            j+=1

    return f_mat

def process_forces(bar_len, diff_x):
    f_mat = gen_beam_mat(bar_len, diff_x)
    
    punt_num = int(input("Ingrese la cantidad de cargas puntuales: "))
    inf_punt = np.zeros(shape=(2, punt_num), dtype=float)
    i = 1
    while punt_num:
        punt_num -= 1   
        force_loc = aprox_diff(float(input(f"Ingrese en donde quiere localizar su fuerza puntual {i}: ")), diff_x)
        force_value = float(input(f"Ingrese la fuerza puntual {i} en Newtons: "))#poner catch
        index = get_idx(force_loc, diff_x) 
        f_mat[1, index] += np.array([force_value]) 

        inf_punt[0, i-1] = force_loc
        inf_punt[1, i-1] = force_value
        i += 1
     

    distr_num = int(input("¿Cuántas cargas distribuidas desea poner?: "))
    inf_dist = [[],[],[],[]]
    k = 1
    while distr_num:
        distr_num -= 1  

        #Display de funciones
        math_expr_str = input(f"Ingrese la función la carga distribuida {k} (en N/m) , con la debida notación de Python: ")
        math_expr = str_to_function(math_expr_str)
        math_expr_wo_dis = str_to_function_wo_dis(math_expr_str)
        begin_fun = aprox_diff(float(input(f"Ingrese la posicion inicial de la carga distribuida {k}: ")), diff_x)
        end_fun = aprox_diff(float(input(f"Ingrese la posición final de la carga distribuida {k}: ")), diff_x)
        dom_f = aprox_diff(float(input(f"Ingrese el valor inicial del dominio de la función {k}: ")), diff_x)
        dist_diff = end_fun - begin_fun
        row_diff = np.arange(dom_f ,dom_f + dist_diff + diff_x, diff_x)
        y_images = math_expr(row_diff)
        ini_idx = get_idx(begin_fun, diff_x)

        inf_dist[0].append(math_expr_wo_dis)
        inf_dist[1].append(begin_fun)
        inf_dist[2].append(end_fun)
        inf_dist[3].append(dom_f)

        k += 1


        j=0
        for i in range(ini_idx, len(y_images)+ini_idx):
            f_mat[1, i] += y_images[j]
            j+=1

    return f_mat, inf_punt, inf_dist


def shift_mat(mat, s_value, diff):
    shift = get_idx(s_value, diff)
    new_ini_p = (-1) * diff * shift
    new_fnl_p = mat[-1] + new_ini_p
    shifted_mat = np.arange(new_ini_p, new_fnl_p + diff, diff)
    return shifted_mat

def shift_mat_mod(mat, s_value):
    func = lambda x: x - s_value
    shifted_mat = func(mat)
    return shifted_mat


row_diff = np.arange(0 ,b + dis, dis)

p1 = shift_mat_mod(row_diff, 4)
p2 = shift_mat_mod(row_diff, 10)
p3 = shift_mat_mod(row_diff, 0)

p = process_forces(b, dis)
print(p[0])
print(p[1])
print(p[2])

f = p[2][0][0]

print(f(5))