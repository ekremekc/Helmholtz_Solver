from dolfin import *

from helmholtz_pkg.mshr import MeshXDMF, MeshXML
import numpy as np
from rijke_geom import geom_rectangle
from petsc4py import PETSc
from helmholtz_pkg.passive_flame import PassiveFlame, passive_flame
from helmholtz_pkg.eigensolvers import fixed_point_iteration_pep, pep_solver
from helmholtz_pkg.eigenvectors import normalize_eigenvector
from helmholtz_pkg.petsc4py_utils import vector_matrix_vector
# Circle Generation by GMSH Python API

geometry = MeshXML('MeshDir/rijkeax', write_xml_file=True)
geometry()
mesh = geometry.mesh
boundaries = geometry.boundaries
subdomains = geometry.subdomains

r = Expression("x[1]", degree = 1)
V = FunctionSpace(mesh,'CG',1)
area = interpolate(Constant(1),V)

dx = Measure('dx', subdomain_data=subdomains)
ds = Measure('ds', domain = mesh, subdomain_data=boundaries)
volume = assemble(2*np.pi*area*r*dx(0))
volume_all = assemble(2*np.pi*area*r*dx(1))



surface_area = assemble( 2 * np.pi * Constant(1) * r * ds(2))