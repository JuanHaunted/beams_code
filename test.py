import numpy as np
from functions import *

b = 20
dis = 0.1

def shift_mat_mod(mat, s_value):
    func = lambda x: x - s_value
    shifted_mat = func(mat)
    return shifted_mat


row_diff = np.arange(0 ,b + dis, dis)
print(len(row_diff))
p1 = shift_mat_mod(row_diff, 4)
p2 = shift_mat_mod(row_diff, 10)
p3 = shift_mat_mod(row_diff, 0)

print(row_diff)
print(p1)

print(len(p1))
print(len(p2))
print(len(p3))