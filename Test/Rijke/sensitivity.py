import dolfin as dolf
from math import pi, sqrt
import matplotlib.pyplot as plt

from helmholtz_pkg.passive_flame import PassiveFlame
from helmholtz_pkg.flame_transfer_function import n_tau
from helmholtz_pkg.active_flame import ActiveFlame
from helmholtz_pkg.eigensolvers import fixed_point_iteration_pep
from helmholtz_pkg.fixed_point_iteration import fixed_point_iteration
from helmholtz_pkg.eigenvectors import normalize_adjoint, normalize_eigenvector

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
                       c=params.c_,
                       degree=degree)
    foo.assemble_A()
    foo.assemble_B()
    foo.assemble_B_astuple()
    foo.assemble_C()
    foo.assemble_C_astuple()

    ftf = n_tau(params.n, params.tau)

    D = ActiveFlame(mesh, subdomains,
                    params.x_f, params.x_r, params.rho_in, 1., 1., ftf,
                    degree=degree)
    D.assemble_submatrices()

    E = fixed_point_iteration(foo, D, pi**2, nev=2, i=0, two_sided=True)

    omega_dir, p_dir = normalize_eigenvector(mesh, E, i=0, degree=degree)
    omega_adj, p_adj = normalize_eigenvector(mesh, E, i=1, degree=degree)

    p_adj = normalize_adjoint(omega_dir, p_dir, p_adj, foo, D)


def ufl_a(mesh, degree):

    CG = dolf.FiniteElement('CG', mesh.ufl_cell(), degree = degree)
    W  = dolf.FunctionSpace(mesh, CG*CG)
    dofmap = W.dofmap()

    
    (u_1, u_2) = dolf.TrialFunction(W)
    (v_1, v_2) = dolf.TestFunction(W)
    dx = dolf.Measure('dx', domain = mesh)
    da_11 = -2 * params.c_ * dolf.dot(dolf.grad(v_1), dolf.grad(u_1)) * dx
    da_22 = -2 * params.c_ * dolf.dot(dolf.grad(v_2), dolf.grad(u_2)) * dx
    
    a = da_11 + da_22

    return dofmap, a

def sensitivities(mesh, degree, p_dir, p_adj):

    dofmap, a = ufl_a(mesh, degree)




if __name__ == '__main__':

    import time

    t0 = time.time()

    run()

    print('eigenvalues computed in {:.3f} sec'.format(time.time() - t0))
