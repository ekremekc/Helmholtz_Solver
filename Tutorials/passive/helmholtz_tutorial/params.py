"""
A note on the nomenclature:
dim ~ dimensional quantity
ref ~ reference quantity for the non-dimensionalization
in ~ inlet, same as u ~ upstream (of the flame)
out ~ outlet, same as d ~ downstream (of the flame)
"""

from math import *
import numpy as np

import dolfin as dolf

# ------------------------------------------------------------

L_ref = 1.  # [m]

# ------------------------------------------------------------

r = 287.  # [J/kg/K]
gamma = 1.4  # [/]

p_amb = 1e5  # [Pa]
rho_amb = 1.22  # [kg/m^3]

T_amb = p_amb/(r*rho_amb)  # [K]

c_amb = sqrt(gamma*p_amb/rho_amb)  # [m/s]

# ------------------------------------------------------------

rho_in_dim = rho_amb  # [kg/m^3]
rho_out_dim = 0.85  # [kg/m^3]

# print('rho_in_dim = %f' % rho_in_dim)
# print('rho_out_dim = %f' % rho_out_dim)

T_in_dim = p_amb/(r*rho_in_dim)  # [K]
T_out_dim = p_amb/(r*rho_out_dim)  # [K]

# print('T_in_dim = %f' % T_in_dim)
# print('T_out_dim = %f' % T_out_dim)

c_in_dim = sqrt(gamma*p_amb/rho_in_dim)  # [kg/m^3]
c_out_dim = sqrt(gamma*p_amb/rho_out_dim)  # [kg/m^3]

# print('c_in_dim = %f' % c_in_dim)
# print('c_out_dim = %f' % c_out_dim)

# ------------------------------------------------------------
# Reflection coefficients

R_in = - 0.975 - 0.05j  # [/]
R_out = - 0.975 - 0.05j  # [/]

# Acoustic impedance

# Z_in = rho_amb*c_amb*(1 + R_in)/(1 - R_in)
# Z_out = rho_amb*c_amb*(1 + R_out)/(1 - R_out)

# print('Z_in =', Z_in)
# print('Z_out =', Z_out)

# Specific impedance

Z_in = (1 + R_in)/(1 - R_in)
Z_out = (1 + R_out)/(1 - R_out)

# print('Z_in =', Z_in)
# print('Z_out =', Z_out)

# Specific admittance

Y_in = 1/Z_in
Y_out = 1/Z_out

# ------------------------------------------------------------

x_f_dim = np.array([[0.25, 0., 0.]])  # [m]
a_f_dim = 0.025  # [m]

# ------------------------------------------------------------
# ------------------------------------------------------------
# Non-dimensionalization

U_ref = c_amb  # [m/s]
p_ref = p_amb  # [Pa]

# ------------------------------------------------------------

rho_in = rho_in_dim*U_ref**2/p_ref
rho_out = rho_out_dim*U_ref**2/p_ref

# print('rho_in = %f' % rho_in)
# print('rho_out = %f' % rho_out)

T_in = T_in_dim*r/U_ref**2
T_out = T_out_dim*r/U_ref**2

# print('T_in = %f' % T_in)
# print('T_out = %f' % T_out)

c_in = c_in_dim/U_ref
c_out = c_out_dim/U_ref

# print('c_in = %f' % c_in)
# print('c_out = %f' % c_out)

# ------------------------------------------------------------

x_f = x_f_dim/L_ref
a_f = a_f_dim/L_ref

# ------------------------------------------------------------

c = dolf.Expression('x[0] <= x_f ? c_in : c_out', degree=0, x_f=x_f[0][0], c_in=c_in, c_out=c_out)
