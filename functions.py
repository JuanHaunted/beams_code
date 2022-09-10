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
        return lambda x: 0.01 * (eval(math_expr))
    else:
        return lambda x: (x*0) + (eval(math_expr) * 0.01)

def str_to_function_wo_dis(math_expr):
    if "x" in math_expr:
        return lambda x:  (eval(math_expr))
    else:
        return lambda x: (x*0) + (eval(math_expr))
    
def get_idx(xi, diff):
    return int(round(xi/diff))

def gen_beam_mat(b_len, x_diff):
    int_len = b_len/x_diff
    f = lambda x: x*x_diff
    row_diff = np.arange(0 ,int_len + 1, 1)
    row_diff = f(row_diff)
    empty_row = np.zeros(len(row_diff))
    return np.stack([row_diff, empty_row])

def shift_mat(mat, s_value):
    func = lambda x: x - s_value
    shifted_mat = func(mat)
    return shifted_mat

def moment_sum(mat, react_indx, mat_pure_moments = np.zeros(0)): #Editar para que funcione con empotradas
    #Rf = (sumatoria(momentos)+sumatoria(momentos_puros))/dist_rf
    dist_rf = mat[0, react_indx]
    neg_mat_pure_moments = mat_pure_moments[1] * (-1)
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
        sup_i_idx = -1
        react_y = force_sum_y(beam_mat[1])
        react_ang = moment_sum(beam_mat, sup_i_idx, pure_moments)
        #pure_moments[1, -1] += react_ang  
        beam_mat[1, -1] -= react_y


        return [beam_mat, pure_moments]

    else:
        sup_i_idx = get_idx(sup_i, diff)
        sup_f_idx = get_idx(sup_f, diff)
        dist_mat = shift_mat(beam_mat[0], sup_i)
        print(len(beam_mat[0]))
        print(len(dist_mat))
        sh_beam_mat = np.stack((dist_mat, beam_mat[1]))
        print(sh_beam_mat)
        react_sf = moment_sum(sh_beam_mat, sup_f_idx, pure_moments)
        beam_mat[1, sup_f_idx] -= react_sf 
        react_si = force_sum_y(beam_mat[1]) 
        beam_mat[1, sup_i_idx] -= react_si
        print(np.array([react_si, react_sf]))

        return [beam_mat, pure_moments]

def initial_values():
    bar_len = float(input("Ingrese la longitud de la viga: "))
    diff = bar_len/(bar_len*10)
    bar_type = int(input("Ingrese el tipo de viga \n-0 => Con dos soportes \n-1 => Empotrada \nSu selección: "))
    if bar_type == 0:

        sup_1 = aprox_diff(float(input("Donde quiere localizar su soporte 1: ")), diff)
        while sup_1<0 or sup_1>bar_len:
            print("Su localización de soporte está por fuerza de la barra")
            sup_1 = aprox_diff(float(input("Donde quiere localizar su soporte 1: ")), diff)

        sup_2 = aprox_diff(float(input("Donde quiere localizar su soporte 2: ")), diff)
        while sup_2<0 or sup_2>bar_len:
            print("Su localización de soporte está por fuerza de la barra")
            sup_2 = aprox_diff(float(input("Donde quiere localizar su soporte 2: ")), diff)

        if sup_1 > sup_2:
            temp = sup_1
            sup_1 = sup_2
            sup_2 = temp

        return np.array([sup_1,sup_2,bar_len, bar_type]) 

    if bar_type == 1:
        print("El análisis se realizará automáticamente con la barra empotrada a la derecha")
        sup_1 = 0

        return np.array([sup_1, bar_len, bar_type])

def pure_moments(bar_len, diff_x):
    moment_mat = gen_beam_mat(bar_len, diff_x)
    flag = input("¿Desea ingresar momentos puros? \n-y => Si \n-n => No \nSu selección: ")
    if flag == 'n':
        return moment_mat
    n = int(input("¿Cuántos momentos puros desea ingresar?: "))

    c = 1
    while n:
        n -= 1
        loc = aprox_diff(float(input(f"Ingresa la localización del momento {c} : ")), diff_x)
        moment_idx = get_idx(loc, diff_x)
        value = float(input(f"Ingrese la magnitud del momento {c} (N.m): "))
        moment_mat[1, moment_idx] += value
        c =+ 1

    return moment_mat

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

def add_pure_moments(moment_mat, bending_mat, b_len, x_diff):
    row_diff = np.arange(0 ,b_len + x_diff, x_diff)
    bending_mat = np.stack([row_diff, bending_mat])
    sum = 0
    for i in range(len(bending_mat[0])):
        sum -= moment_mat[1, i]
        bending_mat[1, i] += sum
    return bending_mat[1]

def aditional_info():
    decision = input("Desea considerar el perfil de la viga? \n -y => Si \n -n => No \n Su selección: ")

    if decision == "n":
        my_S = 0
        my_Q = 0 
        my_I = 0
        my_t = 0 
        my_r = 0
        decision2 = 0
    elif decision == "y":
        decision2 = input("Ingrese el perfil que quiere para su viga \n -r => Rectangular \n -w => Tipo H \n -c => Circular \n Su selección:")
        if decision2 == "w":
            print("El perfil a considerar será de tipo W, el cual es una viga tipo H. \n Revise el anexo de excel llamado SpecsTipo_H. ")
            my_I = float(input("Ingrese el momento de inercia respecto a x (I) (Columna DQ): "))
            my_Q = float(input("Ingrese el primer momento respecto al eje centroidal del área de la sección transversal (Q) (Columna EI): "))
            my_t = float(input("Ingrese el espesor del alma de la viga (t) (Columna CU): "))
            my_S = float(input("Ingrese el módulo de sección (S) (Columna DS): "))
            my_r = 0
        if decision2 == "c":
            print("El perfil a considerar será de tipo circular. ")
            my_I = 1
            my_Q = 1
            my_t = 1
            my_r = float(input("Ingrese el radio del círculo: "))
            my_S = float(input("Ingrese el módulo de sección (S): "))
        elif decision2 == "r":
            decision3 = input("El perfil a considerar será de tipo rectangular \n ¿Quiere calcular usted mismo las especificaciones de la viga? \n -y => Si \n -n => No \n Su selección:")
            if decision3 == "y":
                print("Revise el anexo tipo PDF llamado especificaciones \n Siendo t el espesor de la viga y h su altura, calcule como se muestra. ")
                my_t = float(input("Ingrese el espesor de la viga (t): "))
                my_h = float(input("Ingrese la altura de la viga (h): "))
                my_I = float(input("Ingrese el momento de inercia respecto a x (I): "))
                my_Q = float(input("Ingrese el primer momento respecto al eje centroidal del área de la sección transversal (Q): "))
                my_S = float(input("Ingrese el módulo de sección de la viga (S): "))
                my_r = 0
            elif decision3 == "n":
                print("Revise el anexo tipo PDF llamado especificaciones \n Siendo t el espesor de la viga y h su altura, ingrese los datos. ")
                my_t = float(input("Ingrese el espesor de la viga (t): "))
                my_h = float(input("Ingrese la altura de la viga (h): "))
                my_I = ((my_t)*(my_h**3))/12
                my_S = my_I/(my_h/2)
                my_r = 0
                print("El primer momento respecto al eje centroidal lo calculamos nosotros esta vez ;) ")
                my_Q = (my_t*(my_h**2))/8
    return np.array([my_S, my_Q, my_I, my_t, my_r, decision, decision2]) 

def max_sigma(moments, section):
    moment = np.zeros(2)
    moment[0] = abs(max(moments))
    moment[1] = abs(min(moments))
    max_moment = max(moment)
    sigma_max = (max_moment)/(section)
    return sigma_max
    

def max_tao(sheers, Q, I, t, type, r):
    sheer = np.zeros(2)
    sheer[0] = abs(max(sheers))
    sheer[1] = abs(min(sheers))
    max_sheer = max(sheer)
    if type == "r" or "w":
        tao_max = (max_sheer*Q)/(I*t)
    elif type == "c":
        tao_max = (4*max_sheer)/(pi*(r**2))
    return tao_max
