from math import pi, cos, sin
import gmsh
import os 
if not os.path.exists('MeshDir'):
    os.makedirs('MeshDir')
# Fix file path
os.chdir(os.path.dirname(os.path.realpath(__file__)))

def fltk_options():

    # Type of entity label (0: description,
    #                       1: elementary entity tag,
    #                       2: physical group tag)
    gmsh.option.setNumber("Geometry.LabelType", 2)

    gmsh.option.setNumber("Geometry.PointNumbers", 0)
    gmsh.option.setNumber("Geometry.LineNumbers", 0)
    gmsh.option.setNumber("Geometry.SurfaceNumbers", 2)
    gmsh.option.setNumber("Geometry.VolumeNumbers", 2)

    # Mesh coloring(0: by element type, 1: by elementary entity,
    #                                   2: by physical group,
    #                                   3: by mesh partition)
    gmsh.option.setNumber("Mesh.ColorCarousel", 2)

    gmsh.option.setNumber("Mesh.Lines", 0)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 0) # CHANGE THIS FLAG TO 0 TO SEE LABELS

    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    gmsh.option.setNumber("Mesh.VolumeFaces", 0)
    
# def merge(file="MeshDir/merge", fltk=False):

file="MeshDir/merge"
fltk=True

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add(__name__)

gmsh.open("MeshDir/upstream.msh")
gmsh.model.geo.synchronize()

gmsh.merge("MeshDir/flame.msh")
gmsh.model.addPhysicalGroup(3, [2], tag=0)
gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.geo.synchronize()

gmsh.merge("MeshDir/downstream.msh")
gmsh.model.geo.synchronize()
gmsh.model.mesh.removeDuplicateNodes()
gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2)
gmsh.write("{}.msh".format(file))

if fltk:
    fltk_options()
    gmsh.fltk.run()

gmsh.finalize()





# if __name__ == '__main__':

#     merge(fltk=True)