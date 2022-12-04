import numpy as np
from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
import matplotlib.pyplot as plt
import scipy.signal
import scipy.interpolate as interpolate
from bisect import bisect_left
def get_quadric(radius,c_R):
    posi_max = np.argmax(c_R)
    c_R_left = c_R[:posi_max]
    radius_left = radius[:posi_max]
    c_R_right = c_R[posi_max:]
    radius_right = radius[posi_max:]

    coef1 = np.polyfit(radius_left,c_R_left,2)
    coef2 = np.polyfit(radius_right,c_R_right,2)
    
    return coef1,coef2,posi_max
def f3(t, x):
    i = bisect_left(t, x)
    if t[i] - x > 0.5:
        i-=1
    return i
a = 0.10591694287633399 
b = [0.03, 0.04, 0.05, 0.06, 0.07 ,0.08, 0.09, 0.1 , 0.11 ,0.12 ,0.13, 0.14, 0.15]

c = f3(b,a)
# print(a.shape,r.shape)
# t,c,_ = interpolate.splrep(r,a,s=0,k=4)
print(c)