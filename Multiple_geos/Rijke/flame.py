from math import pi, cos, sin
import gmsh
import os 
if not os.path.exists('MeshDir'):
    os.makedirs('MeshDir')


def geom_pipe(file="MeshDir/flame", fltk=False):

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
    
    L_flame_start = 0.225
    L_flame_end   = 0.275
    
    geom = gmsh.model.geo

    # CIRCLE 2 ________________________________________________________________________________

    p6 = geom.addPoint(0, 0, L_flame_start, lc)

    p7 = geom.addPoint(R * cos(0), R * sin(0), L_flame_start, lc)
    p8 = geom.addPoint(R * cos(pi/2), R * sin(pi/2), L_flame_start, lc)
    p9 = geom.addPoint(R * cos(pi), R * sin(pi), L_flame_start, lc)
    p10 = geom.addPoint(R * cos(3*pi/2), R * sin(3*pi/2), L_flame_start, lc)

    l5 = geom.addCircleArc(p7, p6, p8)
    l6 = geom.addCircleArc(p8, p6, p9)
    l7 = geom.addCircleArc(p9, p6, p10)
    l8 = geom.addCircleArc(p10, p6, p7)

    ll2 = geom.addCurveLoop([l5, l6, l7, l8])
    
    s2 = geom.addPlaneSurface([ll2], tag=2)

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

    

    # SHELL 2 AND VOLUME 2 (FLAME) ________________________________________________________________________________

    l17 = geom.addLine(p7, p12)
    l18 = geom.addLine(p8, p13)
    l19 = geom.addLine(p9, p14)
    l20 = geom.addLine(p10, p15)

    ll8 = geom.addCurveLoop([l5, l18, -l13, -l17])
    ll9 = geom.addCurveLoop([l6, l19, -l14, -l18])
    ll10 = geom.addCurveLoop([l7, l20, -l15, -l19])
    ll11 = geom.addCurveLoop([l8, l17, -l16, -l20])

    s8  = geom.addSurfaceFilling([ll8], tag=8)
    s9  = geom.addSurfaceFilling([ll9], tag=9)
    s10 = geom.addSurfaceFilling([ll10], tag=10)
    s11 = geom.addSurfaceFilling([ll11], tag=11)

    sl2 = geom.addSurfaceLoop([s2, s8, s9, s10, s11, s7])
    vol2 = geom.addVolume([sl2], tag=2)

def add_physical_entities():

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [8, 9, 10, 11], tag=3)

    gmsh.model.addPhysicalGroup(3, [2], tag=0)     # Flame region


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