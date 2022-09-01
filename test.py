import numpy as np
from functions import *

test = gen_beam_mat(10, 0.1)

test[0, 4] = 5
print(test)