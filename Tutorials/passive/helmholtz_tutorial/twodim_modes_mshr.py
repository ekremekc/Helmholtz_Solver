import dolfin as dolf
from helmholtz_tutorial import params
from helmholtz_tutorial.passive_flame import *
from helmholtz_tutorial.eigensolvers import *
from helmholtz_tutorial.eigenvectors import *
import matplotlib.pyplot as plt
from mshr_pkg.mshr import MeshXDMF
from mshr_pkg.rectangle import geom_rectangle

c_0 = 450
L = 0.4
h = 0.1

"-----------o------- NUMERICAL RESULTS -------o-----------"

nx,ny,N_num = 40,10, 8
Z_right_b = np.linspace(-10j, 10j, N_num)
Z_right_a = np.linspace(-10, 10, N_num)

def mshr_():
    geom_rectangle()
    geometry = MeshXDMF('MeshDir/rectangle', write_xdmf_file=True)
    geometry()
    mesh = geometry.mesh
    boundaries = geometry.boundaries
    
    return mesh, boundaries

mesh, boundaries = mshr_()

def onedim_eig(Z, c_0, L, target, index=1):
    
    boundary_conditions = {1: {'Neumann'},
                           2: {'Neumann'},
                           3: {'Robin': 1/Z},
                           4: {'Neumann'}}
    degree = 2
    
    operators = PassiveFlame(mesh, boundaries, boundary_conditions,
                             c=dolf.Constant(1.0), degree=degree)
    operators.assemble_A()
    operators.assemble_B()
    operators.assemble_C()
    A = operators.A
    B = operators.B
    C = operators.C
    target = target
    # target = np.sqrt((np.pi/L)**2)
    E = pep_solver(A, B, C, target, 2, print_results=False)
    
    vr, vi = A.createVecs()
    k = E.getEigenpair(index, vr, vi)
    # print(k)
    # FOR 2D ACOUSTIC CASE
    
    omega = k*c_0
    f = omega/(2*np.pi)
    
    return f , k.real

def get_eig_freqs(Z):
    # k0 = 8.16
    omega_real = []
    omega_imag = []
    if isinstance(Z[0], complex)==False:
        k0 = 11
        for i in range(N_num):
            if Z[i]<0:
                omega, k_ = onedim_eig(Z[i], c_0, L, k0, index = 0)
                # k0 = k_
            else:
                omega, k_ = onedim_eig(Z[i], c_0, L, k0)
                # k0 = k_
            print(omega)
            omega_real.append(omega.real)
            omega_imag.append(omega.imag)
    else:
        k0 = 7
        for i in range(N_num):
            if Z[i]<0:
                omega, k_ = onedim_eig(Z[i], c_0, L, k0, index=1)
                k0 = k_
            else:
                k0 = 6
                omega, k_ = onedim_eig(Z[i], c_0, L, k0, index=1)
                k0 = k_
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

analytic_b_real, analytic_b_imag = [], []
analytic_a_real, analytic_a_imag = [], []

f=open('analytical.txt',"r")
lines=f.readlines()
for x in lines:
    analytic_b_real.append(float(x.split(' ')[0]))
    analytic_b_imag.append(float(x.split(' ')[1]))
    analytic_a_real.append(float(x.split(' ')[2]))
    analytic_a_imag.append(float(x.split(' ')[3]))
f.close()

x_analytic = Z_a

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12,10))
# fig.suptitle('One-Dimensional Eigenmodes')
ax1.plot(x_analytic, analytic_b_real)
ax1.plot(x_num, numeric_b_real,'^')
ax1.set_xlabel(r"$b$")
ax1.set_ylabel(r"Re{$f$}")

ax2.plot(Z_b, analytic_b_imag)
ax2.plot(x_num, numeric_b_imag,'^')
ax2.set_ylim([-0.1, 0.1])
ax2.set_xlabel(r"$b$")
ax2.set_ylabel(r"Im{$f$}")


ax3.plot(Z_a, analytic_a_real)
ax3.plot(x_num, numeric_a_real,'^')
ax3.set_xlabel(r"$a$")
ax3.set_ylabel(r"Re{$f$}")


ax4.plot(Z_a, analytic_a_imag)
ax4.plot(x_num, numeric_a_imag,'^')
ax4.set_xlabel(r"$a$")
ax4.set_ylabel(r"Im{$f$}")

# plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                    wspace=0.3, hspace=None)
# left  = 0.125  # the left side of the subplots of the figure
# right = 0.9    # the right side of the subplots of the figure
# bottom = 0.1   # the bottom of the subplots of the figure
# top = 0.9      # the top of the subplots of the figure
# wspace = 0.2   # the amount of width reserved for blank space between subplots
# hspace = 0.2   # the amount of height reserved for white space between subplots
plt.show()

plt.savefig("2dmode.pdf")
