# from helmholtz_tutorial.main import *
import dolfin as dolf
from helmholtz_tutorial import params
from helmholtz_tutorial.passive_flame import *
from helmholtz_tutorial.eigensolvers import *
from helmholtz_tutorial.eigenvectors import *

import numpy as np
import matplotlib.pyplot as plt
print(params.R_out)
c_0 = 450
L = 0.4

"-----------o------- NUMERICAL RESULTS -------o-----------"

N_num = 40
Z_right_b = np.linspace(-10j, 10j, N_num)
Z_right_a = np.linspace(-10, 10, N_num)

def mshr_(nx, L):
    """
    :param nx: number of elements
    :return: mesh and boundaries
    """
    
    class Left(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], 0)

    class Right(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], L)

    left = Left()
    right = Right()

    # mesh = dolf.UnitIntervalMesh(nx)
    mesh = dolf.IntervalMesh(nx, 0, L)

    boundaries = dolf.MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)

    left.mark(boundaries, 1)
    right.mark(boundaries, 2)

    return mesh, boundaries

def onedim_eig(Z, c_0, L, index=1):

    degree = 1

    mesh, boundaries = mshr_(400, L)

    boundary_conditions = {1: {'Neumann'},
                            2: {'Robin': 1/Z}}

    A, B, C = passive_flame_mixed(mesh, boundaries, boundary_conditions,
                                  c=dolf.Constant(450.0), degree=degree)
    # A, B, C = passive_flame_mixed(mesh, boundaries, boundary_conditions,
    #                               c=dolf.Constant(c_0), degree=degree)
    
    target_value = np.pi/L*450
    
    E = pep_solver(A, B, C, target_value, 2, print_results=False)
    
    vr, vi = A.createVecs()
    k = E.getEigenpair(index, vr, vi)
    
    # FOR 1D ACOUSTIC CASE
    
    # omega = c_0*k/(L_omega)
    # f = omega/(2*np.pi)
    
    omega = k
    f = omega/(2*np.pi)
    
    return f

def get_eigvector(Z, c_0, L, index=1):

    degree = 1

    mesh, boundaries = mshr_(400, L)

    boundary_conditions = {1: {'Neumann'},
                            2: {'Robin': 1/Z}}

    A, B, C = passive_flame_mixed(mesh, boundaries, boundary_conditions,
                                  c=dolf.Constant(1.0), degree=degree)
    # A, B, C = passive_flame_mixed(mesh, boundaries, boundary_conditions,
    #                               c=dolf.Constant(c_0), degree=degree)
    
    E = pep_solver(A, B, C, np.pi, 2, print_results=False)

    vr, vi = A.createVecs()
    k = E.getEigenpair(index, vr, vi)
    
    # FOR 1D ACOUSTIC CASE
    print(k)
    omega, p = normalize_eigenvector(mesh, E, 1, degree=degree)

    p_r, p_i = p.split(True)

    dolf.plot(p_i)
    plt.xlabel('x')
    plt.ylabel('p_r')
    plt.show()

# get_eigvector(-10, c_0, L)

def get_eig_freqs(Z):
    omega_real = []
    omega_imag = []
    if isinstance(Z[0], complex)==False:
        for i in range(N_num):
            if Z[i]<0:
                omega = onedim_eig(Z[i], c_0, L, index = 0)
            else:
                omega = onedim_eig(Z[i], c_0, L)
            omega_real.append(omega.real)
            omega_imag.append(omega.imag)
    else:
        for i in range(N_num):
            omega = onedim_eig(Z[i], c_0, L)
            omega_real.append(omega.real)
            omega_imag.append(omega.imag)
        
    return omega_real, omega_imag

x_num =  np.linspace(-10, 10, N_num)
numeric_b_real, numeric_b_imag = get_eig_freqs(Z_right_b)
numeric_a_real, numeric_a_imag = get_eig_freqs(Z_right_a)

"-----------o------- ANALYTICAL RESULTS -------o-----------"

N_analytic = 400
Z_b = np.linspace(-10, 10, N_analytic)
Z_a = Z_b

def reactive(Z_b, c_0, L, m=1):
    f_m = m*(c_0/(2*L)) + (c_0/((2*L)*(np.pi)))*np.arctan(-1/Z_b)
    return f_m

def resistive(Z_a, c_0, L, m=1):
    f_m = m*(c_0/(2*L)) + (c_0/((2*L)*(np.pi)))*np.arctan(-1j/Z_a)
    return f_m

x_analytic = Z_a
analytic_b_real, analytic_b_imag = reactive(Z_b, c_0, L), np.zeros(N_analytic)
analytic_a_real, analytic_a_imag = resistive(Z_a, c_0, L).real, resistive(Z_a, c_0, L).imag

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12,10))
fig.suptitle('One-Dimensional First Eigenmodes for Z=a+bi')
ax1.plot(x_analytic, analytic_b_real)
ax1.plot(x_num, numeric_b_real,'^')
ax1.set_xlabel(r"$b$ (pure imaginary, a=0)")
ax1.set_ylabel(r"Re{$f$}")

ax2.plot(Z_b, analytic_b_imag)
ax2.plot(x_num, numeric_b_imag,'^')
ax2.set_ylim([-0.1, 0.1])
ax2.set_xlabel(r"$b$ (pure imaginary, a=0)")
ax2.set_ylabel(r"Im{$f$}")


ax3.plot(Z_a, analytic_a_real)
ax3.plot(x_num, numeric_a_real,'^')
ax3.set_xlabel(r"$a$ (pure real, b=0)")
ax3.set_ylabel(r"Re{$f$}")


ax4.plot(Z_a, analytic_a_imag)
ax4.plot(x_num, numeric_a_imag,'^')
ax4.set_xlabel(r"$a$ (pure real, b=0)")
ax4.set_ylabel(r"Im{$f$}")

# plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                    wspace=0.3, hspace=0.3)
# left  = 0.125  # the left side of the subplots of the figure
# right = 0.9    # the right side of the subplots of the figure
# bottom = 0.1   # the bottom of the subplots of the figure
# top = 0.9      # the top of the subplots of the figure
# wspace = 0.2   # the amount of width reserved for blank space between subplots
# hspace = 0.2   # the amount of height reserved for white space between subplots
plt.show()

plt.savefig("1dmode.pdf")