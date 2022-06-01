import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root, minimize
from mpmath import findroot
from pylab import *

c_0 = 450
k = 3.1416006868531317+0.10033638047403654j
# k = 562.5014383502627+17.965159789940266j
f = c_0*k/(2*np.pi)
Z = -10

L = 0.4
h = 0.1

def twoD_nicoud(k, Z=Z, L=L, h=h, m=1):
    k_y = np.sqrt(k**2-((m*np.pi)/L)**2)
    func = np.exp(2j*k_y*h)*(k_y-k/Z)-(k_y+k/Z)
    print(func)
    return func

def scipy_solve(k, Z=Z, L=L, h=h, m=1):
    k_y = np.sqrt(k[0]**2-((m*np.pi)/L)**2)
    func = np.exp(2j*k_y*h)*(k_y-k[0]/Z)-(k_y+k[0]/Z)
    # print(func)
    return func

def mp_solve(k, Z=Z, L=L, h=h, m=1):
    k_y = sqrt(k**2-((m*np.pi)/L)**2)
    func = exp(2j*k_y*h)*(k_y-k/Z)-(k_y+k/Z)
    # print(func)
    return func


# should_be_zero = twoD_nicoud(k)

root = minimize(scipy_solve,[0.1, 0.1], method='Nelder-Mead')

# root = findroot(mp_solve,[1+1j])


print(root)






