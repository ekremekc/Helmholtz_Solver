from mshr_pkg.mshr import MeshXDMF
from mshr_pkg.rectangle import geom_rectangle
import dolfin as dolf
geom_rectangle('MeshDir/rectangle', fltk=False)

geometry = MeshXDMF('MeshDir/rectangle', write_xdmf_file=True)
geometry()
mesh = geometry.mesh
boundaries = geometry.boundaries

dolf.plot(mesh)