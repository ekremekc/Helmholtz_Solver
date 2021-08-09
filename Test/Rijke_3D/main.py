import dolfin as dolf
from math import pi
import numpy as np
from helmholtz_solver.helmholtz_pkg.passive_flame import PassiveFlame
from helmholtz_solver.helmholtz_pkg.flame_transfer_function import n_tau
from helmholtz_solver.helmholtz_pkg.active_flame import ActiveFlame
from helmholtz_solver.helmholtz_pkg.eigensolvers import fixed_point_iteration_pep
from helmholtz_solver.helmholtz_pkg.eigenvectors import normalize_eigenvector
from helmholtz_solver.helmholtz_pkg.mshr import MeshXDMF
import params 


def mshr():

    geometry = MeshXDMF('MeshDir/rijke3d', write_xdmf_file=True)
    geometry()
    mesh = geometry.mesh
    boundaries = geometry.boundaries
    subdomains = geometry.subdomains

    return mesh, boundaries, subdomains


def run():

    degree = 1

    mesh, boundaries, subdomains = mshr()

    boundary_conditions = {1: {'Robin': params.Y_in},  # inlet
                           2: {'Robin': params.Y_out}, # outlet
                           3: {'Neumann'}}             # wall

    foo = PassiveFlame(mesh, boundaries, boundary_conditions,
                       c=params.c,
                       degree=degree)
    foo.assemble_A()
    foo.assemble_B()
    foo.assemble_C()

    ftf = n_tau(params.n, params.tau)

    D = ActiveFlame(mesh, subdomains,
                     params.x_r, params.rho_in, 1., 1., ftf,
                    degree=degree)
    D.assemble_submatrices()

    E = fixed_point_iteration_pep(foo, D, pi, nev=2, i=0)

    omega, p = normalize_eigenvector(mesh, E, i=0, degree=degree)

    p_r, p_i = p.split(True)
    
    
    p_r_file   = dolf.File("Results/p_r.pvd")
    p_i_file   = dolf.File("Results/p_i.pvd")
    p_r_file   << p_r
    p_i_file   << p_i

    p_abs = p_r.copy()
    p_abs.vector()[:] = np.abs(p_r.vector()[:]+1j*p_i.vector()[:])
    p_abs_file = dolf.File("Results/p_abs.pvd")
    p_abs_file << p_abs

if __name__ == '__main__':

    import time

    t0 = time.time()

    run()

    print('eigenvalues computed in {:.3f} sec'.format(time.time() - t0))
