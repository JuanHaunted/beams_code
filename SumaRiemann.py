import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import sympy as smp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import cumulative_trapezoid



def f(x):
    y = (1/2)*(x**2)
    return (y)

def Riemann(f,xi,xf,n):
    x = np.linspace(xi, xf, n)             #Valores de x en los n intervalos
    A = 0                                  #Aproximación del area bajo la curva. Inicia en cero
    a = []                                 #Vector que guarda el area del rectangulo en cada subintervalo
    f1 = []                                #Valores de la funcion a la izq de cada subintervalo
    f2 = []
    xbar = []
    ybar = []

    for i in range (1,n):
        #Suma de RIemman
        f1.append(f(x[i-1]))
        f2.append(f(x[i]))
        fx = np.amin([f(x[i-1]), f(x[i])])
        deltax = x[i] - x[i-1]
        a.append(fx*deltax)
        A += fx*deltax
        print(A)
        #Datos gráficos de Barras:
        #1) Esquina inferior izquierda
        xbar.append(x[i-1])
        ybar.append(0)
        #2) Esquina superior izquierda
        xbar.append(x[i-1])
        ybar.append(fx)
        #3) Esquina superior derecha
        xbar.append(x[i])
        ybar.append(fx)
        #4) Esquina inferior derecha
        xbar.append(x[i])
        ybar.append(0)

    return(A, xbar, ybar)


xi = -20
xf = 20
n1 = 50

A1, xbar, ybar = Riemann(f, xi, xf, n1)

x = np.linspace(xi, xf, 500)

print("Grafica de la funcion")
plt.plot(x, f(x), 'k', label=("f(x"))
plt.plot(xbar,ybar, 'b:', label=("Suma de Riemann Inferior"))
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Grafica de f(x)")
plt.legend()
plt.show()
