from mshr_pkg.mshr import MeshXDMF
from mshr_pkg.square import geom_square
from mshr_pkg.pipe import geom_pipe
from mshr_pkg.rectangle import geom_rectangle


# The first argument is the relative path of the msh file, without extension.
# fltk=True for visualization

# geom_pipe('MeshDir/pipe', fltk=True)

# geom_square('MeshDir/square', fltk=True)

geom_rectangle('MeshDir/rectangle', fltk=True)

# You can use MeshXDMF as a black box.
# The first argument is the relative path of the msh file, without extension (same as geom_pipe).
# If write_xdmf_file=False, MeshXDMF will check if 'MeshDir/pipe.xdmf' exists.

geometry = MeshXDMF('MeshDir/rectangle', write_xdmf_file=True)
geometry()
mesh = geometry.mesh
boundaries = geometry.boundaries
# subdomains = geometry.subdomains
