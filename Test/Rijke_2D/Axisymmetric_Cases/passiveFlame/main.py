import dolfin as dolf
import matplotlib.pyplot as plt
from helmholtz_pkg.mshr import MeshXDMF, MeshXML
import numpy as np

from helmholtz_solver.helmholtz_pkg.passive_flame import PassiveFlame
from helmholtz_solver.helmholtz_pkg.eigensolvers import pep_solver
from helmholtz_solver.helmholtz_pkg.eigenvectors import normalize_eigenvector
from helmholtz_solver.helmholtz_pkg.petsc4py_utils import vector_matrix_vector
import params

geometry = MeshXDMF('MeshDir/rijke')
geometry()
mesh = geometry.mesh
boundaries = geometry.boundaries

degree = 1

boundary_conditions = {1: {'Neumann'},
                       2: {'Robin': params.Y_out},
                       3: {'Neumann'},
                       4: {'Robin': params.Y_in}}

matrices = PassiveFlame(mesh, boundaries, boundary_conditions,
                   c=params.c_,
                   degree=degree)

matrices.assemble_A()
matrices.assemble_B()
matrices.assemble_C()

E = pep_solver(matrices.A, matrices.B, matrices.C, target = np.pi, nev=2, print_results=True)

omega, p = normalize_eigenvector(mesh, E, i=0, degree=degree) 

p_r, p_i = p.split(True)

plot1 = dolf.plot(p_r) ; plt.colorbar(plot1)
plt.show()
plot2 = dolf.plot(p_i) ; plt.colorbar(plot2)
plt.show()