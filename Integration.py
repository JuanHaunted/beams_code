import numpy as np
import scipy as sp
import sympy as smp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import cumulative_trapezoid

# https://github.com/lukepolson/youtube_channel/blob/main/Python%20Tutorial%20Series/integrals1.ipynb
# https://www.youtube.com/watch?v=2I44Y9hfQ4Q


# SOLVABLE INTEGRALS 
x = smp.symbols('x', real=True)
f = smp.sin(x)**3 * smp.exp(-5*x)
smp.integrate(f, x)


a, b = smp.symbols('a b', real=True, positive=True)
f = smp.cos(b*x)* smp.exp(-a*x)
smp.integrate(f, x).simplify()


f = (1+smp.sqrt(x))**smp.Rational(1,3) / smp.sqrt(x)
smp.integrate(f, x).simplify()


f = smp.exp(x) / smp.sqrt(smp.exp(2*x) + 9)
smp.integrate(f, (x, 0, smp.log(4)))


f = 16*smp.atan(x) / (1+x**2)
smp.integrate(f, (x, 0, smp.oo))



# UNSOLVABLE INTEGRALS
f = smp.exp(-smp.sin(x))
#smp.integrate(f, (x, 1, 2)).simplify()
#So we use scipy's quad function to integrate numerically

f = lambda x: np.exp(-np.sin(x))
f = lambda x: np.exp(-np.sin(x))
quad(f, 1, 2)

#Ver bien el video, muy teso



# NUMERICAL CASE

x, y = np.loadtxt('../data/coviddata.txt')

plt.plot(x,y)
plt.xlabel('Day')
plt.ylabel('Cases per Day')
plt.show()

integral = cumulative_trapezoid(y,x, initial=0)

plt.plot(x,integral)
plt.xlabel('Day')
plt.ylabel('Cumulative Cases')
plt.grid()
plt.show()