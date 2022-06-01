from mshr_pkg.mshr import MeshXDMF
from mshr_pkg.rectangle import geom_rectangle
import dolfin as dolf
geom_rectangle('meshDir/rectangle', fltk=True)

geometry = MeshXDMF('meshDir/rectangle', write_xdmf_file=True)
geometry()
mesh = geometry.mesh
boundaries = geometry.boundaries

dolf.plot(mesh)