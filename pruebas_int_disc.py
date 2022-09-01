from re import A, X
from math import sin,cos,pi
import numpy as np
import scipy as sp
import sympy as smp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import cumulative_trapezoid


#Datos daods
a = 1              #Limite inferior
b = 10              #Limite superior
n = 10           #Particiones 
dif_x = (b-a)/(n)   #Diferencial de x

print(dif_x)

#Estructura de datos
covid = np.zeros((2,n)).astype(float)
my_x=np.linspace(a, b, n)


#Relleno mis listas
for i in range (0, n):
    covid[0,i]= my_x[i]

    i+=i

for i in range (0, n):
    covid[1,i]= (-1/2)*(my_x[i])

    i+=i


print(covid[0,:])  #Es lo mismo q printear  print(my_x)
print(covid[1,:])


#Graficación e integración con trapezoides
plt.plot(covid[0,:], covid[1,:])
plt.xlabel('Variable independiente')
plt.ylabel('Variable dependiente')
plt.show()
integral = cumulative_trapezoid(covid[0,:],covid[1,:], initial=0)
plt.plot(covid[0,:], integral)
plt.xlabel('Variable independiente')
plt.ylabel('Variable dependiente acumulada')
plt.grid()
plt.show()


#Valor integral
A = 0

#for i in range(1, 100):




#Prueba real con simpy y porcentaje de error
x = smp.symbols('x', real=True)
f = (x+1)**(1/2)
smp.integrate(f, (x, a, b))

error_ = (abs((smp.integrate(f, (x, a, b))) - A)/(smp.integrate(f, (x, a, b))))*100


#Print para pruebas
print(integral)
print(error_)