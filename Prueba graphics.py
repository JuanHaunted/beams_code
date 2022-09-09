import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import style
import matplotlib.markers
#import seaborn as sns

y = [-4,2,-2,-1,0,1,-3,3,4,5]
x = np.arange(0,10,1)

print(x)
print(y)

with plt.style.context('bmh'):
    plt.plot(x,y,label='Función cortante')
    plt.xlabel('Puntos de la viga')
    plt.ylabel('Cortantes')
    plt.suptitle("Gráfica de cortantes", )
    plt.minorticks_on()

    plt.grid(which='major', linestyle='-', linewidth='1.5', color='white')
    plt.grid(which='minor', linestyle=':', linewidth='1.0', color='white')


    plt.tick_params(which='both', # Options for both major and minor ticks
                    left='on', # turn off left ticks
                    bottom='off') # turn off bottom ticks


    maximo = abs(max(x))
    punto_y = y[maximo]




    plt.annotate((maximo, punto_y) ,xy =(maximo, punto_y), xytext = (maximo, punto_y))
    plt.plot(maximo,punto_y, marker="D", color="green", label='Cortante máximo' )

    #plt.legend([maximo, punto_y],label="Cortante máximo")
    plt.legend(fancybox=True, loc='upper left', shadow=True)

    plt.show()

