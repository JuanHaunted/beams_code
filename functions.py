from operator import mod
from math import *
from tkinter import X
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

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
    row_diff = np.arange(0, int_len + 1, 1)
    row_diff = f(row_diff)
    empty_row = np.zeros_like(row_diff)
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
        while sup_1<0 or sup_1>bar_len+0.0001:
            print("Su localización de soporte está por fuerza de la barra")
            sup_1 = aprox_diff(float(input("Donde quiere localizar su soporte 1: ")), diff)

        sup_2 = aprox_diff(float(input("Donde quiere localizar su soporte 2: ")), diff)
        while sup_2<0 or sup_2>bar_len+0.0001:
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
    inf_dist = [[], [], [], [], []]
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
        int_len = bar_len/diff_x
        int_dom = dom_f/diff_x
        int_dist = dist_diff/diff_x
        f = lambda x: x*diff_x
        row_diff = np.arange(int_dom, int_dom + int_dist + 1, 1)
        row_diff = f(row_diff)
        #row_diff = np.arange(dom_f ,dom_f + dist_diff + diff_x, diff_x)
        y_images = math_expr(row_diff)
        ini_idx = get_idx(begin_fun, diff_x)

        inf_dist[0].append(math_expr_wo_dis)
        inf_dist[1].append(begin_fun)
        inf_dist[2].append(end_fun)
        inf_dist[3].append(dom_f)
        inf_dist[4].append(math_expr_str)
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
    int_len = b_len/x_diff
    f = lambda x: x*x_diff
    row_diff = np.arange(0, int_len + 1, 1)
    row_diff = f(row_diff)
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
            decision3 = input("El perfil a considerar será de tipo circular\n ¿Quiere calcular usted mismo las especificaciones de la viga? \n -y => Si \n -n => No \n Su selección:")
            if decision3 == "y":
                print("Revise el anexo tipo PDF llamado especificaciones \n Siendo r el radio del círculo. Ingrese:  ")
                my_Q = 1
                my_I = 1
                my_Q = 1
                my_t = 1
                my_r = float(input("Ingrese el radio del círculo: "))
                my_S = float(input("Ingrese el módulo de sección (S): "))
            elif decision3 == "n":
                print("Revise el anexo tipo PDF llamado especificaciones \n Siendo r el radio del círculo. Ingrese:  ")
                my_Q = 1
                my_t = 1
                my_r = float(input("Ingrese el radio del círculo: "))
                my_I = ((1/8)*pi)*(my_r**4)
                my_S = my_I/my_r
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
                print("Revise el anexo tipo PDF llamado especificaciones \n Siendo t el espesor de la viga y h su altura, ingrese según la tabla. ")
                my_t = float(input("Ingrese el espesor de la viga (t): "))
                my_h = float(input("Ingrese la altura de la viga (h): "))
                my_I = ((my_t)*(my_h**3))/12
                my_S = my_I/(my_h/2)
                my_r = 0
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
    tao_max = 0
    if type == "r" or "w":
        tao_max = (max_sheer*Q)/(I*t)
    elif type == "c":
        tao_max = (4*max_sheer)/(pi*(r**2))
    return tao_max

# Función para hallar los puntos inferiores de los soportes, para graficar
def find_point(z):
    bajo_der = [z[0] + 1, z[1] - 3]
    bajo_iz = [z[0] - 1, z[1] - 3]
    return bajo_der, bajo_iz


# Función para la gráfica ilustrativa de fuerzas en una viga
def gen_graph(bar_l, bar_type, x_sup, der_sup, iz_sup, x2_sup, der_sup2, iz_sup2, forc_m, distrib, moments):
    with plt.style.context('ggplot'):
        if bar_type == 0:

            # Diseño gráfica
            plt.xlim([-2, bar_l + 2])
            plt.ylim([-8, 30])
            plt.suptitle("Gráfica de cargas sobre una barra apoyada", fontsize=12, fontweight='bold',
                         color="darkslateblue")
            plt.xlabel("Distancia (m)")

            # Se dibuja la barra, modulo patches rectangle
            # dar coordenada de inicio, largo en y, longitud en x y grosor. zorder pone adelante la figura
            rect = mpatches.Rectangle((0, -3), bar_l, 3, zorder=2, color="lightslategrey")
            plt.gca().add_patch(rect)

            # Soporte1, modulo patches crea figuras de n cantidad de lados
            list_trian = [x_sup, der_sup, iz_sup]
            trian = mpatches.Polygon(list_trian, zorder=2, color="darkseagreen")
            plt.gca().add_patch(trian)

            # soporte2
            list_trian2 = [x2_sup, der_sup2, iz_sup2]
            trian2 = mpatches.Polygon(list_trian2, zorder=2, color="darkseagreen")
            plt.gca().add_patch(trian2)

            # Variables fuerzas puntuales
            loc_f = forc_m[0]
            value_f = forc_m[1]

            # Variables distribuidas

            math_exp = distrib[0]
            begin_f = distrib[1]
            end_f = distrib[2]
            domain = distrib[3]
            str_m = distrib[4]

            # Variables momento

            loc_mom = moments[0]
            value_mom = moments[1]

            # Dibujar fuerzas puntuales (flechas)
            for i in range(0, len(loc_f)):
                col = (np.random.random(), np.random.random(), np.random.random())  # genera colores random

                # dibuja una flecha con coordenadas iniciales a finales
                puntual = plt.arrow(loc_f[i], 5, 0, -3.5, zorder=2, color=col, width=0.3, label=f"{value_f[i]} N")

            # Dibujar distribuidas

            for x in range(0, len(math_exp)):
                # condicional para graficar dominios diferentes a 0
                if domain[x] != 0:
                    # se añade la traslación
                    xpts = np.linspace(begin_f[x] + domain[x], end_f[x] + domain[x], 50)
                    expr = math_exp[x]
                    evalu = expr(xpts)
                    col = (np.random.random(), np.random.random(), np.random.random())

                else:
                    xpts = np.linspace(begin_f[x], end_f[x], 50)
                    expr = math_exp[x]
                    evalu = expr(xpts)
                    col = (np.random.random(), np.random.random(), np.random.random())  # genera colores random

                plt.plot(xpts, evalu, zorder=2, color=col, label=str_m[x])
                plt.fill_between(xpts, 0, evalu, color=col, alpha=0.2)

            # Dibujar momentos

            # ellip = mpatches.FancyArrowPatch((4, 3), (4, 5), connectionstyle="arc3,rad=.5")
            # plt.gca().add_patch(ellip)

            plt.legend(fancybox=True, loc='upper right', shadow=True)
            plt.subplots_adjust(left=0.135,
                                bottom=0.125,
                                right=0.905,
                                top=0.835,
                                wspace=0.2,

                                hspace=0.205)
        plt.minorticks_on()  # se activan rayas en los ejes para dar presición a las medidas
        plt.tick_params(which='major',  # Options for both major and minor ticks
                        bottom='off')  # turn off bottom ticks

        plt.grid(color='white')

        plt.show()


def cool_graphs(mat, mat_m, dis, bar_l):
    x = mat[0]
    y = mat[1]

    # Encontrar coordenadas xy y del cortante max
    max_ref = np.max(y)
    abs_cort = np.absolute(y)
    max_cort = np.max(abs_cort)

    if max_ref < max_cort:
        max_cort *= -1

    punto_xc = np.where(y == max_cort)
    punto_x = x[punto_xc[0]][0]
    integral_list = riemann_sum(dis, y)

    if np.sum(mat_m[1]) != 0:
        integral_list = add_pure_moments(mat_m, integral_list, bar_l, dis)

    # Encontrar coordenadas x y del flexionante máximo
    max_refc = np.max(integral_list)
    abs_flex = np.absolute(integral_list)
    max_flex = np.max(abs_flex)

    if max_refc < max_flex:
        max_flex *= -1

    punto_xf = np.where(integral_list == max_flex)
    punto_xflex = x[punto_xf[0]][0]


    with plt.style.context('ggplot'):
        # se crea la ventana de subplots
        fig, axs = plt.subplots(2)
        fig.suptitle('Gráficas de análisis de viga', fontsize=12, fontweight='bold', color="mediumslateblue")

        # grafica cortantes
        func = axs[0].plot(np.concatenate([[0], x]), np.concatenate([[0], y]), color="red", label="función cortante")
        axs[0].set_xlabel('Puntos de la viga (m)')
        axs[0].set_ylabel('Cortantes (N)')
        axs[0].grid(which='major', linestyle='-', linewidth='1.2', color='white')
        axs[0].set_title("Gráfica de fuerzas cortantes", fontsize=9, fontweight='bold', color="darkslategrey")
        axs[0].minorticks_on()
        axs[0].tick_params(which='major',  # Options for both major and minor ticks
                           left='on',  # turn off left ticks
                           bottom='off')  # turn off bottom ticks
        # Dibujar punto máximo
        maxi = round(max_cort, 1)
        axs[0].annotate((punto_x, round(max_cort, 2)), xy=(punto_x, max_cort), xytext=(punto_x, round(max_cort, 2)))
        axs[0].plot(punto_x, max_cort, marker="D", color="green", label=f'Cortante máximo\n{maxi}')
        axs[0].fill_between(x, 0, y, color="red", alpha=0.2)  # colorear entre función y eje x
        axs[0].legend(func, [punto_x, max_cort])
        axs[0].legend(fancybox=True, loc='center left', shadow=True, bbox_to_anchor=((1, 0.5)))

        # grafica flexionantes
        axs[1].plot(x, integral_list, color="navy", label="función flex.")
        axs[1].set_xlabel('Puntos de la viga (m)')
        axs[1].set_ylabel('Flexionantes (N)')
        axs[1].grid(which='major', linestyle='-', linewidth='1.2', color='white')
        axs[1].set_title("Gráfica de fuerzas flexionantes", fontsize=9, fontweight='bold', color="darkslategrey")
        axs[1].minorticks_on()
        axs[1].tick_params(which='major',  # Options for both major and minor ticks
                           left='on',  # turn off left ticks
                           bottom='off')  # turn off bottom ticks

        axs[1].fill_between(x, 0, integral_list, color="blue", alpha=0.2)

        # punto max flexionante
        maxi_f = round(max_flex, 1)
        axs[1].annotate((punto_xflex, round(max_flex, 2)), xy=(punto_xflex, max_flex),
                        xytext=(punto_xflex, round(max_flex, 2)))
        axs[1].plot(punto_xflex, max_flex, marker="D", color="salmon", label=f'flex. máximo\n{maxi_f}')
        axs[1].legend(func, [punto_xflex, max_flex])
        axs[1].legend(fancybox=True, loc='center left', shadow=True, bbox_to_anchor=((1, 0.5)))

        plt.subplots_adjust(left=0.13,
                            bottom=0.15,
                            right=0.7,
                            top=0.87,
                            wspace=0.675,
                            hspace=0.695)

        plt.show(block=False)
        plt.pause(9900)