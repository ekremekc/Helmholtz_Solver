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

# p1 = plot_2d(p_r, mesh, orientation = "H", save = True)
# p2 = plot_2d(p_i, mesh, orientation = "H")

import matplotlib.tri as tri
from mpl_toolkits.axes_grid1 import make_axes_locatable

def mesh2triang(mesh):
        xy = mesh.coordinates()
        return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())
fig = plt.figure(figsize=(8,4))

ax1 = fig.add_subplot(121)
C1 = p_r.compute_vertex_values(mesh)
tpc1 = ax1.tripcolor(mesh2triang(mesh), C1)
divider = make_axes_locatable(ax1)
cax = divider.new_vertical(size="5%", pad=0.7, pack_start=True)
fig.add_axes(cax)
fig.colorbar(tpc1, cax=cax, orientation="horizontal")

ax2 = fig.add_subplot(122)
C2 = p_i.compute_vertex_values(mesh)
tpc2 = ax2.tripcolor(mesh2triang(mesh), C2)
divider2 = make_axes_locatable(ax2)
cax2 = divider2.new_vertical(size="5%", pad=0.7, pack_start=True)
fig.add_axes(cax2)
fig.colorbar(tpc2, cax=cax2, orientation="horizontal")

fig.tight_layout()
fig.savefig("2D.png", dpi=300)

