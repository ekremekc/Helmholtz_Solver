import dolfin as dolf
from math import pi, sqrt
import matplotlib.pyplot as plt

from helmholtz_pkg.passive_flame import PassiveFlame
from helmholtz_pkg.flame_transfer_function import n_tau
from helmholtz_pkg.active_flame import ActiveFlame
from helmholtz_pkg.eigensolvers import fixed_point_iteration_pep
from helmholtz_pkg.eigenvectors import normalize_eigenvector
from helmholtz_pkg.petsc4py_utils import vector_matrix_vector
import params


def mshr(el):

    mesh = dolf.UnitIntervalMesh(el)

    def l_boundary_func(x, on_boundary):
        x = x[0]
        return on_boundary and dolf.near(x, 0.)

    def r_boundary_func(x, on_boundary):
        x = x[0]
        return on_boundary and dolf.near(x, 1.)

    boundaries = dolf.MeshFunction('size_t', mesh, mesh.topology().dim() - 1)

    l_boundary = dolf.AutoSubDomain(l_boundary_func)
    r_boundary = dolf.AutoSubDomain(r_boundary_func)

    l_boundary.mark(boundaries, 1)
    r_boundary.mark(boundaries, 2)

    # ________________________________________________________________________________

    def fl_subdomain_func(x):
        x = x[0]
        x_f = params.x_f[0][0]
        a_f = params.a_f
        return x_f - a_f - dolf.DOLFIN_EPS <= x <= x_f + a_f + dolf.DOLFIN_EPS

    subdomains = dolf.MeshFunction('size_t', mesh, mesh.topology().dim())

    subdomains.set_all(1)

    fl_subdomain = dolf.AutoSubDomain(fl_subdomain_func)
    fl_subdomain.mark(subdomains, 0)

    return mesh, boundaries, subdomains


def run():

    degree = 1

    mesh, boundaries, subdomains = mshr(400)

    boundary_conditions = {1: {'Robin': params.Y_in},  # inlet
                           2: {'Robin': params.Y_out}}  # outlet

    foo = PassiveFlame(mesh, boundaries, boundary_conditions,
                       c=params.c,
                       degree=degree)
    foo.assemble_A()
    foo.assemble_B()
    foo.assemble_C()

    ftf = n_tau(params.n, params.tau)

    D = ActiveFlame(mesh, subdomains,
                    params.x_f, params.x_r, params.rho_in, 1., 1., ftf,
                    degree=degree)
                    
    D.assemble_submatrices(problem_type='direct')
    D.assemble_submatrices(problem_type='adjoint')

    E_dir = fixed_point_iteration_pep(foo, D, pi, nev=2, i=0, problem_type='direct')
    E_adj = fixed_point_iteration_pep(foo, D, pi, nev=2, i=1, problem_type='adjoint')

    omega_dir, p_dir = normalize_eigenvector(mesh, E_dir, i=0, degree=degree, which='right')
    omega_adj, p_adj = normalize_eigenvector(mesh, E_adj, i=1, degree=degree, which='left')

    # p_dir_r, p_dir_i = p_dir.split(True)
    # p_adj_r, p_adj_i = p_adj.split(True)

    # BASE STATE SENSITIVITY ------------------------------------------------------------------------------------------

    p_dir_vec = p_dir.vector().vec()
    p_adj_vec = p_adj.vector().vec()
    
    dLds = -D.get_derivative(omega_dir) + foo.B + foo.assemble_zC(2 * omega_dir)

    denominator = -vector_matrix_vector(p_adj_vec , dLds , p_dir_vec )
    #print(denominator)

    
    D.assemble_matrix(omega_dir)
    # sensitivity for n
    dNdn = D.matrix / params.n
    nom_n = -vector_matrix_vector(p_adj_vec , dNdn , p_dir_vec )
    dsdn = nom_n / denominator
    print("Sensitivity of ds/dn: ",dsdn)

    # sensitivity for tau
    D.assemble_matrix(1j*omega_dir)
    dNdtau = D.matrix 
    nom_tau = -vector_matrix_vector(p_adj_vec , dNdtau , p_dir_vec )
    dsdtau = nom_tau / denominator
    print("Sensitivity of ds/dtau: ",dsdtau)

    dtau = 0.01
    ftf = n_tau(params.n, params.tau + dtau)

    D = ActiveFlame(mesh, subdomains,
                    params.x_f, params.x_r, params.rho_in, 1., 1., ftf,
                    degree=degree)
                    
    D.assemble_submatrices(problem_type='direct')
    D.assemble_submatrices(problem_type='adjoint')

    E_dir = fixed_point_iteration_pep(foo, D, pi, nev=2, i=0, problem_type='direct')
    E_adj = fixed_point_iteration_pep(foo, D, pi, nev=2, i=1, problem_type='adjoint')

    omega_dir_1, p_dir_1 = normalize_eigenvector(mesh, E_dir, i=0, degree=degree, which='right')
    omega_adj_1, p_adj_1 = normalize_eigenvector(mesh, E_adj, i=1, degree=degree, which='left')

    print(omega_dir_1, " ?= ", omega_dir+dsdtau*dtau)
    print("n0 = ", params.tau, " n1 = ", params.tau+dtau)


    # pl, ax = plt.subplots(); fig = plt.gcf(); fig.set_size_inches(10, 3.6)
    # plt.subplot(1, 2, 1); p1 = dolf.plot(p_adj_r); plt.xlabel('x');plt.title('Real Part')
    # plt.subplot(1, 2, 2); p2 = dolf.plot(-p_adj_i); plt.xlabel('x');plt.title('Imaginary Part')
    # plt.tight_layout()
    # plt.show()

if __name__ == '__main__':

    import time

    t0 = time.time()

    run()

    print('eigenvalues computed in {:.3f} sec'.format(time.time() - t0))
