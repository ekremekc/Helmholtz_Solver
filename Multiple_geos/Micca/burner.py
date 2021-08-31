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


def burner(file, fltk=False, **kwargs):

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    print(gmsh.option.getNumber("Mesh.FirstElementTag"))
    gmsh.model.add(__name__)

    # geom = gmsh.model.geo

    add_elementary_entities(**kwargs)
    
    apply_symmetry()
    
    apply_rotation()

    add_physical_entities()
    
    gmsh.model.geo.synchronize()

    
    
    gmsh.model.mesh.generate(3)
    
    
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4)
    gmsh.write("{}.msh".format(file))
    print(gmsh.model.getEntities(dim=2))
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
              'lc_1': 1e-1,
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

    p11 = geom.addPoint(R_mid, 0., - l_pp - l_b, lc_2,5000)
    p12 = geom.addPoint(R_mid + h_b, 0., - l_pp - l_b, lc_2)
    p13 = geom.addPoint(R_mid, h_b, - l_pp - l_b, lc_2)
    p14 = geom.addPoint(R_mid - h_b, 0., - l_pp - l_b, lc_2)

    l10 = geom.addLine(p14, p12,5000)
    l15 = geom.addCircleArc(p12, p11, p13)
    l16 = geom.addCircleArc(p13, p11, p14)

    ll7 = geom.addCurveLoop([l10, l15, l16])
    s7 = geom.addPlaneSurface([ll7],5000)


    
    # BURNER

    p15 = geom.addPoint(R_mid, 0., - l_pp, lc_2)

    p16 = geom.addPoint(R_mid + h_b, 0., - l_pp, lc_2)
    p17 = geom.addPoint(R_mid, h_b, - l_pp, lc_2)
    p18 = geom.addPoint(R_mid - h_b, 0., - l_pp, lc_2)

    p19 = geom.addPoint(R_mid + h_pp, 0., - l_pp, lc_2)
    p20 = geom.addPoint(R_mid, h_pp, - l_pp, lc_2)
    p21 = geom.addPoint(R_mid - h_pp, 0., - l_pp, lc_2)

    l17 = geom.addLine(p14, p18)
    l18 = geom.addLine(p12, p16)
    l19 = geom.addLine(p13, p17)

    l20 = geom.addLine(p18, p21)
    l21 = geom.addLine(p21, p19)
    l22 = geom.addLine(p19, p16)

    l23 = geom.addCircleArc(p16, p15, p17)
    l24 = geom.addCircleArc(p17, p15, p18)

    l25 = geom.addCircleArc(p19, p15, p20)
    l26 = geom.addCircleArc(p20, p15, p21)

    ll8 = geom.addCurveLoop([-l22, -l21, -l20, -l17, l10, l18])
    s8 = geom.addPlaneSurface([ll8])

    ll9 = geom.addCurveLoop([-l23, -l18, l15, l19])
    s9 = geom.addSurfaceFilling([ll9])

    ll10 = geom.addCurveLoop([-l24, -l19, l16, l17])
    s10 = geom.addSurfaceFilling([ll10])

    ll11 = geom.addCurveLoop([l20, -l26, -l25, l22, l23, l24])
    s11 = geom.addPlaneSurface([ll11])

    ll12 = geom.addCurveLoop([l21, l25, l26])
    s12 = geom.addPlaneSurface([ll12])

    sl2 = geom.addSurfaceLoop([s7, s8, s9, s10, s11, s12])
    vol2 = geom.addVolume([sl2],5000)

    
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

        gmsh.model.mesh.setPeriodic(2, new_entities, old_entities, rotation_matrix(angle))
        my_surf_old = my_surf_new  
        gmsh.model.geo.synchronize()

def add_physical_entities():
    a=gmsh.model.getEntities(dim=3)
    b_lateral = []
    b_front = []

    for volume in a:
        value = gmsh.model.getAdjacencies(3, volume[1])
        b_lateral.append(value[1][2])
        b_lateral.append(value[1][3])
        b_front.append(value[1][4])
        
    gmsh.model.addPhysicalGroup(2, b_lateral, tag=5)
    gmsh.model.addPhysicalGroup(2, b_front, tag=6)

    
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
    gmsh.option.setNumber("Geometry.VolumeNumbers", 2)

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

    burner('MeshDir/burner', fltk=True, l_ec=0.041)