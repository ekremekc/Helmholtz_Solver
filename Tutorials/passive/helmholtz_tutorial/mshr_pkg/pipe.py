from math import pi, cos, sin
import gmsh


def geom_pipe(file="MeshDir/pipe", fltk=False):

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

    L = 1.
    R = 0.025
    lc = 1e-2

    geom = gmsh.model.geo

    # ________________________________________________________________________________

    p1 = geom.addPoint(0, 0, 0, lc)

    p2 = geom.addPoint(R * cos(0), R * sin(0), 0, lc)
    p3 = geom.addPoint(R * cos(pi/2), R * sin(pi/2), 0, lc)
    p4 = geom.addPoint(R * cos(pi), R * sin(pi), 0, lc)
    p5 = geom.addPoint(R * cos(3*pi/2), R * sin(3*pi/2), 0, lc)

    l1 = geom.addCircleArc(p2, p1, p3)
    l2 = geom.addCircleArc(p3, p1, p4)
    l3 = geom.addCircleArc(p4, p1, p5)
    l4 = geom.addCircleArc(p5, p1, p2)

    ll1 = geom.addCurveLoop([l1, l2, l3, l4])

    s1 = geom.addPlaneSurface([ll1])

    # ________________________________________________________________________________

    p6 = geom.addPoint(0, 0, L, lc)

    p7 = geom.addPoint(R * cos(0), R * sin(0), L, lc)
    p8 = geom.addPoint(R * cos(pi/2), R * sin(pi/2), L, lc)
    p9 = geom.addPoint(R * cos(pi), R * sin(pi), L, lc)
    p10 = geom.addPoint(R * cos(3*pi/2), R * sin(3*pi/2), L, lc)

    l5 = geom.addCircleArc(p7, p6, p8)
    l6 = geom.addCircleArc(p8, p6, p9)
    l7 = geom.addCircleArc(p9, p6, p10)
    l8 = geom.addCircleArc(p10, p6, p7)

    ll2 = geom.addCurveLoop([l5, l6, l7, l8])

    s2 = geom.addPlaneSurface([ll2])

    # ________________________________________________________________________________

    l9 = geom.addLine(p2, p7)
    l10 = geom.addLine(p3, p8)
    l11 = geom.addLine(p4, p9)
    l12 = geom.addLine(p5, p10)

    ll3 = geom.addCurveLoop([l1, l10, -l5, -l9])
    ll4 = geom.addCurveLoop([l2, l11, -l6, -l10])
    ll5 = geom.addCurveLoop([l3, l12, -l7, -l11])
    ll6 = geom.addCurveLoop([l4, l9, -l8, -l12])

    s3 = geom.addSurfaceFilling([ll3])
    s4 = geom.addSurfaceFilling([ll4])
    s5 = geom.addSurfaceFilling([ll5])
    s6 = geom.addSurfaceFilling([ll6])

    sl1 = geom.addSurfaceLoop([s1, s2, s3, s4, s5, s6])
    vol1 = geom.addVolume([sl1])


def add_physical_entities():

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [1], tag=1)
    gmsh.model.addPhysicalGroup(2, [2], tag=2)
    gmsh.model.addPhysicalGroup(2, [3, 4, 5, 6], tag=3)

    gmsh.model.addPhysicalGroup(3, [1], tag=1)


def fltk_options():

    # Type of entity label (0: description,
    #                       1: elementary entity tag,
    #                       2: physical group tag)
    gmsh.option.setNumber("Geometry.LabelType", 2)

    gmsh.option.setNumber("Geometry.PointNumbers", 0)
    gmsh.option.setNumber("Geometry.LineNumbers", 0)
    gmsh.option.setNumber("Geometry.SurfaceNumbers", 0)
    gmsh.option.setNumber("Geometry.VolumeNumbers", 0)

    # Mesh coloring(0: by element type, 1: by elementary entity,
    #                                   2: by physical group,
    #                                   3: by mesh partition)
    gmsh.option.setNumber("Mesh.ColorCarousel", 2)

    gmsh.option.setNumber("Mesh.Lines", 0)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 1)

    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    gmsh.option.setNumber("Mesh.VolumeFaces", 0)


if __name__ == '__main__':

    geom_pipe(fltk=True)
