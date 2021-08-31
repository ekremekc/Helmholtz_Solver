from math import pi, cos, sin
import gmsh
import os 
if not os.path.exists('MeshDir'):
    os.makedirs('MeshDir')


def geom_pipe(file="MeshDir/downstream", fltk=False):

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(__name__)

    add_elementary_entities()
    add_physical_entities()

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.write("{}.msh".format(file))

    if fltk:
        fltk_options()
        gmsh.fltk.run()

    gmsh.finalize()


def add_elementary_entities():

    
    R = 0.047/2
    lc = 5e-3

    L_flame_end   = 0.275
    L_total = 1.
    
    
    geom = gmsh.model.geo

   

    # CIRCLE 3 ________________________________________________________________________________

    p11 = geom.addPoint(0, 0, L_flame_end, lc)

    p12 = geom.addPoint(R * cos(0), R * sin(0), L_flame_end, lc)
    p13 = geom.addPoint(R * cos(pi/2), R * sin(pi/2), L_flame_end, lc)
    p14 = geom.addPoint(R * cos(pi), R * sin(pi), L_flame_end, lc)
    p15 = geom.addPoint(R * cos(3*pi/2), R * sin(3*pi/2), L_flame_end, lc)

    l13 = geom.addCircleArc(p12, p11, p13)
    l14 = geom.addCircleArc(p13, p11, p14)
    l15 = geom.addCircleArc(p14, p11, p15)
    l16 = geom.addCircleArc(p15, p11, p12)

    ll7 = geom.addCurveLoop([l13, l14, l15, l16])

    s7 = geom.addPlaneSurface([ll7], tag=7)


    # CIRCLE 4 ________________________________________________________________________________

    p16 = geom.addPoint(0, 0, L_total, lc)

    p17 = geom.addPoint(R * cos(0), R * sin(0), L_total, lc)
    p18 = geom.addPoint(R * cos(pi/2), R * sin(pi/2), L_total, lc)
    p19 = geom.addPoint(R * cos(pi), R * sin(pi), L_total, lc)
    p20 = geom.addPoint(R * cos(3*pi/2), R * sin(3*pi/2), L_total, lc)

    l21 = geom.addCircleArc(p17, p16, p18)
    l22 = geom.addCircleArc(p18, p16, p19)
    l23 = geom.addCircleArc(p19, p16, p20)
    l24 = geom.addCircleArc(p20, p16, p17)

    ll12 = geom.addCurveLoop([l21, l22, l23, l24])

    s12 = geom.addPlaneSurface([ll12], tag=12)

    # SHELL 3 AND VOLUME 3 (DOWNSTREAM) ________________________________________________________________________________

    l25 = geom.addLine(p12, p17)
    l26 = geom.addLine(p13, p18)
    l27 = geom.addLine(p14, p19)
    l28 = geom.addLine(p15, p20)

    ll13 = geom.addCurveLoop([l13, l26, -l21, -l25])
    ll14 = geom.addCurveLoop([l14, l27, -l22, -l26])
    ll15 = geom.addCurveLoop([l15, l28, -l23, -l27])
    ll16 = geom.addCurveLoop([l16, l25, -l24, -l28])

    s13  = geom.addSurfaceFilling([ll13], tag=13)
    s14  = geom.addSurfaceFilling([ll14], tag=14)
    s15 = geom.addSurfaceFilling([ll15], tag=15)
    s16 = geom.addSurfaceFilling([ll16], tag=16)

    sl3 = geom.addSurfaceLoop([s7, s13, s14, s15, s16, s12])
    vol3 = geom.addVolume([sl3], tag=3)

def add_physical_entities():

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [12], tag=2)
    gmsh.model.addPhysicalGroup(2, [13, 14, 15, 16], tag=3)

    gmsh.model.addPhysicalGroup(3, [3], tag=999) # Upstream and Downstream



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

    gmsh.option.setNumber("Mesh.VolumeEdges", 1)
    gmsh.option.setNumber("Mesh.VolumeFaces", 1)


if __name__ == '__main__':

    geom_pipe(fltk=True)