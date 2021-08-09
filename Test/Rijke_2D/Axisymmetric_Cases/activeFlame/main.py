import dolfin as dolf
import matplotlib.pyplot as plt
from helmholtz_pkg.mshr import MeshXDMF
import numpy as np

from helmholtz_solver.helmholtz_pkg.passive_flame import PassiveFlame
from helmholtz_solver.helmholtz_pkg.active_flame import ActiveFlame
from helmholtz_solver.helmholtz_pkg.flame_transfer_function import n_tau
from helmholtz_solver.helmholtz_pkg.eigensolvers import  fixed_point_iteration_pep
from helmholtz_solver.helmholtz_pkg.eigenvectors import normalize_eigenvector
from helmholtz_solver.helmholtz_pkg.solver_utils import plot_2d

import params

geometry = MeshXDMF('MeshDir/rijke')
geometry()
mesh = geometry.mesh
boundaries = geometry.boundaries
subdomains = geometry.subdomains

degree = 1

boundary_conditions = {1: {'Neumann'},
                       2: {'Robin': params.Y_out},
                       3: {'Neumann'},
                       4: {'Robin': params.Y_in}}

matrices = PassiveFlame(mesh, boundaries, boundary_conditions,
                   c=params.c,
                   degree=degree)

matrices.assemble_A()
matrices.assemble_B()
matrices.assemble_C()

ftf = n_tau(params.n, params.tau)

D = ActiveFlame(mesh, subdomains, params.x_r, params.rho_in, 1, 1, ftf, degree = degree)

D.assemble_submatrices()

E = fixed_point_iteration_pep(matrices, D, target = np.pi, nev=2, i=0, print_results=False)

omega, p = normalize_eigenvector(mesh, E, i=0, degree=degree) 

p_r, p_i = p.split(True)

p1 = plot_2d(p_r, mesh, orientation = "H", save = True)
p2 = plot_2d(p_i, mesh, orientation = "H")

