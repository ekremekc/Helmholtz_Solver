from math import pi, radians, cos, sin
import gmsh
import os

# Fix file path
os.chdir(os.path.dirname(os.path.realpath(__file__)))

def reflection_matrix():
    return [1,  0,  0,  0,
            0, -1,  0,  0,
            0,  0,  1,  0,
            0,  0,  0,  1]


def rotation_matrix(angle):
    c, s = cos(angle), sin(angle)
    return [c, -s,  0,  0,
            s,  c,  0,  0,
            0,  0,  1,  0,
            0,  0,  0,  1]


def plenum(file, fltk=False, **kwargs):

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(__name__)
    
    # geom = gmsh.model.geo

    add_elementary_entities(**kwargs)
    
    apply_symmetry()
    
    apply_rotation()

    add_physical_entities()
   
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    # gmsh.option.setNumber("Mesh.FirstNodeTag", 500.0)
    # gmsh.option.setNumber("Mesh.FirstElementTag", 500.0)
    # gmsh.model.geo.synchronize()
    # gmsh.model.mesh.renumberNodes()
    # gmsh.model.mesh.renumberElements()
    
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4)
    gmsh.write("{}.msh".format(file))

    if fltk:
        fltk_options()
        gmsh.fltk.run()

    gmsh.finalize()


def add_elementary_entities(**kwargs):
    """
    inner : in
    middle:
    outer : out
    p : plenum
    b : burner
    pp: injection
    f : flame
    cc: combustion chamber
    """

    params = {'R_in_p': .14,
              'R_out_p': .21,
              'l_p': .07,
              'h_b': .0165,
              'l_b': .014,
              'h_pp': .00945,
              'l_pp': .006,
              'h_f': .018,
              'l_f': .006,
              'R_in_cc': .15,
              'R_out_cc': .2,
              'l_cc': .2,
              'l_ec': .041,
              'lc_1': 5e-2,
              'lc_2': 1e-2
              }

    for key, value in kwargs.items():
        if key in params.keys():
            params[key] = value

    R_in_p = params['R_in_p']
    R_out_p = params['R_out_p']
    l_p = params['l_p']

    h_b = params['h_b']
    l_b = params['l_b']

    h_pp = params['h_pp']
    l_pp = params['l_pp']

    h_f = params['h_f']
    l_f = params['l_f']

    R_in_cc = params['R_in_cc']
    R_out_cc = params['R_out_cc']
    l_cc = params['l_cc']
    l_ec = params['l_ec']  # end correction

    R_mid = .175  # mid plane

    lc_1 = params['lc_1']
    lc_2 = params['lc_2']
    
    theta = 22.5 / 2
    theta = radians(theta)

    geom = gmsh.model.geo

    # ________________________________________________________________________________
    # PLENUM

    p1 = geom.addPoint(0., 0., - l_pp - l_b - l_p, lc_1)

    p2 = geom.addPoint(R_in_p, 0., - l_pp - l_b - l_p, lc_1)
    p3 = geom.addPoint(R_out_p, 0., - l_pp - l_b - l_p, lc_1)
    p4 = geom.addPoint(R_in_p * cos(theta), R_in_p * sin(theta), - l_pp - l_b - l_p, lc_1)
    p5 = geom.addPoint(R_out_p * cos(theta), R_out_p * sin(theta), - l_pp - l_b - l_p, lc_1)

    l1 = geom.addLine(p2, p3)
    l2 = geom.addCircleArc(p3, p1, p5)
    l3 = geom.addLine(p5, p4)
    l4 = geom.addCircleArc(p4, p1, p2)

    ll1 = geom.addCurveLoop([-l1, -l4, -l3, -l2])
    s1 = geom.addPlaneSurface([ll1])

    p6 = geom.addPoint(0., 0., - l_pp - l_b, lc_1)

    p7 = geom.addPoint(R_in_p, 0., - l_pp - l_b, lc_1)
    p8 = geom.addPoint(R_out_p, 0., - l_pp - l_b, lc_1)
    p9 = geom.addPoint(R_in_p * cos(theta), R_in_p * sin(theta), - l_pp - l_b, lc_1)
    p10 = geom.addPoint(R_out_p * cos(theta), R_out_p * sin(theta), - l_pp - l_b, lc_1)

    l5 = geom.addLine(p2, p7)
    l6 = geom.addLine(p3, p8)
    l7 = geom.addLine(p5, p10)
    l8 = geom.addLine(p4, p9)

    p11 = geom.addPoint(R_mid, 0., - l_pp - l_b, lc_2)

    p12 = geom.addPoint(R_mid + h_b, 0., - l_pp - l_b, lc_2)
    p13 = geom.addPoint(R_mid, h_b, - l_pp - l_b, lc_2)
    p14 = geom.addPoint(R_mid - h_b, 0., - l_pp - l_b, lc_2)

    l9 = geom.addLine(p7, p14)
    l10 = geom.addLine(p14, p12)
    l11 = geom.addLine(p12, p8)
    l12 = geom.addCircleArc(p8, p6, p10)
    l13 = geom.addLine(p10, p9)
    l14 = geom.addCircleArc(p9, p6, p7)

    l15 = geom.addCircleArc(p12, p11, p13)
    l16 = geom.addCircleArc(p13, p11, p14)

    ll2 = geom.addCurveLoop([-l11, -l10, -l9, -l5, l1, l6])
    s2 = geom.addPlaneSurface([ll2])

    ll3 = geom.addCurveLoop([-l6, l2, l7, -l12])
    s3 = geom.addSurfaceFilling([ll3])

    ll4 = geom.addCurveLoop([-l13, -l7, l3, l8])
    s4 = geom.addPlaneSurface([ll4])

    ll5 = geom.addCurveLoop([l5, -l14, -l8, l4])
    s5 = geom.addSurfaceFilling([ll5])

    ll6 = geom.addCurveLoop([l9, -l16, -l15, l11, l12, l13, l14])
    s6 = geom.addPlaneSurface([ll6])

    ll7 = geom.addCurveLoop([l10, l15, l16])
    s7 = geom.addPlaneSurface([ll7])

    sl1 = geom.addSurfaceLoop([s1, s2, s3, s4, s5, s6, s7])
    vol1 = geom.addVolume([sl1])


    
def apply_symmetry():

    symmetry = gmsh.model.geo.symmetrize

    gmsh.model.geo.synchronize()
    old_entities = [ent[1] for ent in gmsh.model.getEntities(dim=2)]

    symmetry(gmsh.model.geo.copy(gmsh.model.getEntities(dim=2)), 0, 1, 0, 0)
    symmetry(gmsh.model.geo.copy(gmsh.model.getEntities(dim=3)), 0, 1, 0, 0)

    gmsh.model.geo.synchronize()
    
    new_entities = [ent[1] for ent in gmsh.model.getEntities(dim=2)]
    new_entities = new_entities[len(old_entities):]
    old_entities.pop(1) # exctract intersection plane
    
    gmsh.model.mesh.setPeriodic(2, new_entities, old_entities, reflection_matrix())
    # print(old_entities, new_entities)
    gmsh.model.geo.synchronize()

def apply_rotation():

    theta = 22.5 / 2
    theta *= pi / 180

    rotate = gmsh.model.geo.rotate

    gmsh.model.geo.synchronize()
    my_surf = gmsh.model.getEntities(dim=2)
    my_surf_old = gmsh.model.getEntities(dim=2)
    
    old_entities = gmsh.model.getEntities(dim=2)  # for periodicity
    old_entities.pop(-4) #Substract intersection plane(bottom of first sample)
    old_entities = [x[1] for x in old_entities]

    my_vol = gmsh.model.getEntities(dim=3)

    for i in range(1, 16):

        angle = 2 * i * theta

        rotate(gmsh.model.geo.copy(my_surf), 0, 0, 0, 0, 0, 1, angle)

        gmsh.model.geo.synchronize()
        my_surf_new = gmsh.model.getEntities(dim=2)
        
        # Difference between two lists
        tmp = set(my_surf_old)
        new_entities = [x for x in my_surf_new if x not in tmp]
        new_entities = [x[1] for x in new_entities]

        rotate(gmsh.model.geo.copy(my_vol), 0, 0, 0, 0, 0, 1, angle)

        if i == 15:
            old_entities.pop(3) #Extract intersection plane(top of first sample)

        gmsh.model.mesh.setPeriodic(2, new_entities, old_entities, rotation_matrix(angle))
        my_surf_old = my_surf_new  
        gmsh.model.geo.synchronize()

def add_physical_entities():
    a=gmsh.model.getEntities(dim=3)
    pl_rear = []
    pl_front = []
    pl_outer = []
    pl_inner = []
    for volume in a:
        value = gmsh.model.getAdjacencies(3, volume[1])
        pl_rear.append(value[1][5])
        pl_front.append(value[1][0])
        pl_outer.append(value[1][2])
        pl_inner.append(value[1][4])
        
        
    gmsh.model.addPhysicalGroup(2, pl_rear, tag=1)
    gmsh.model.addPhysicalGroup(2, pl_front, tag=4)
    gmsh.model.addPhysicalGroup(2, pl_outer, tag=2)
    gmsh.model.addPhysicalGroup(2, pl_inner, tag=3)
    
    volume_tags = [x[1] for x in a]
    gmsh.model.addPhysicalGroup(3, volume_tags, tag=99)
    


def fltk_options():

    # Type of entity label (0: description,
    #                       1: elementary entity tag,
    #                       2: physical group tag)
    gmsh.option.setNumber("Geometry.LabelType", 2)

    gmsh.option.setNumber("Geometry.PointNumbers", 0)
    gmsh.option.setNumber("Geometry.LineNumbers", 0)
    gmsh.option.setNumber("Geometry.SurfaceNumbers", 2)
    gmsh.option.setNumber("Geometry.VolumeNumbers", 0)

    # Mesh coloring(0: by element type, 1: by elementary entity,
    #                                   2: by physical group,
    #                                   3: by mesh partition)
    gmsh.option.setNumber("Mesh.ColorCarousel", 2)

    gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 0)

    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    gmsh.option.setNumber("Mesh.VolumeFaces", 0)

    # my_dict = {'One': (255, 0, 0)}
    # for key, value in my_dict.items():
    #     gmsh.option.setColor('Mesh.{}'.format(key), *value)


if __name__ == '__main__':

    plenum('MeshDir/plenum', fltk=True, l_ec=0.041)