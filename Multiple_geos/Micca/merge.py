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
    

file="MeshDir/Micca"
fltk=True

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
# gmsh.option.setNumber("Mesh.IgnorePeriodicity", 0)
gmsh.model.add(__name__)

gmsh.merge("MeshDir/flame.msh")

flame_vol_tags=gmsh.model.getEntities(dim=3)

for i in range(0, int(len(flame_vol_tags)/2)):
    gmsh.model.addPhysicalGroup(3, [flame_vol_tags[2*i][1],flame_vol_tags[2*i+1][1]], tag=i)
 
gmsh.model.geo.synchronize()

gmsh.merge("MeshDir/plenum.msh")
gmsh.model.geo.synchronize()

gmsh.merge("MeshDir/burner.msh")
gmsh.model.geo.synchronize()

gmsh.merge("MeshDir/injection.msh")
gmsh.model.geo.synchronize()

gmsh.merge("MeshDir/chamber.msh")
gmsh.model.geo.synchronize()

gmsh.model.mesh.removeDuplicateNodes()

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