from helmholtz_tutorial.main import *
import matplotlib.pyplot as plt

c_0 = 450
L = 0.4

"-----------o------- NUMERICAL RESULTS -------o-----------"

N_num = 10
Z_right_b = np.linspace(-10j, 10j, N_num)
Z_right_a = np.linspace(-10, 10, N_num)

def onedim_eig(Z, c_0, L):

    degree = 1

    mesh, boundaries = mshr(400)

    boundary_conditions = {1: {'Neumann'},
                            2: {'Robin': 1/Z}}

    A, B, C = passive_flame_mixed(mesh, boundaries, boundary_conditions,
                                  c=dolf.Constant(1.0), degree=degree)
    
    E = pep_solver(A, B, C, np.pi, 2, print_results=True)

    vr, vi = A.createVecs()
    k = E.getEigenpair(1, vr, vi)
    
    # FOR 1D ACOUSTIC CASE
    
    print("Z = ",Z)
    omega = c_0*k/(2*np.pi*L)
    print("Eigenvalue is",k)
    print("Eigenfrequency is ", omega)
    return omega


def get_eig_freqs(Z):
    omega_real = []
    omega_imag = []
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

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.suptitle('One-Dimensional Eigenmodes')
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

plt.tight_layout()
plt.show()