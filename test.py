import numpy as np
from functions import *

test = empty_row = np.zeros(101)
moments = np.array([[0, 2.7],[90, 50]])

help = add_pure_moments(moments, test, 10, 0.1)
print(help)