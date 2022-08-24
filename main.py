from operator import mod
import numpy as np

def aprox_diff(num, diff_x):
    n_diff = num//diff_x
    residual = num % diff_x
    diff_mean = diff_x/2
    if residual >= diff_mean:
        return (n_diff*diff_x) + diff_x
    else:
        return n_diff * diff_x

n = 100.2
x_diff = 0.5

print(aprox_diff(n, x_diff))