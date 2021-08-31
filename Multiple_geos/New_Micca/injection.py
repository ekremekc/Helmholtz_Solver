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


def injection(file, fltk=False, **kwargs):

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
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
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
    
    # BURNER

    p15 = geom.addPoint(R_mid, 0., - l_pp, lc_2,tag=10000)

    p19 = geom.addPoint(R_mid + h_pp, 0., - l_pp, lc_2)
    p20 = geom.addPoint(R_mid, h_pp, - l_pp, lc_2)
    p21 = geom.addPoint(R_mid - h_pp, 0., - l_pp, lc_2)




    l21 = geom.addLine(p21, p19,tag=10000)
    

    l25 = geom.addCircleArc(p19, p15, p20)
    l26 = geom.addCircleArc(p20, p15, p21)

    

    ll12 = geom.addCurveLoop([l21, l25, l26])
    s12 = geom.addPlaneSurface([ll12],tag=10000)

    
    
    # INJECTION

    p22 = geom.addPoint(R_mid, 0., 0., lc_2)

    p23 = geom.addPoint(R_mid + h_pp, 0., 0., lc_2)
    p24 = geom.addPoint(R_mid, h_pp, 0., lc_2)
    p25 = geom.addPoint(R_mid - h_pp, 0., 0., lc_2)

    l27 = geom.addLine(p21, p25)
    l28 = geom.addLine(p19, p23)
    l29 = geom.addLine(p20, p24)

    l30 = geom.addLine(p25, p23)
    l31 = geom.addCircleArc(p23, p22, p24)
    l32 = geom.addCircleArc(p24, p22, p25)

    ll13 = geom.addCurveLoop([-l30, -l27, l21, l28])
    s13 = geom.addPlaneSurface([ll13])

    ll14 = geom.addCurveLoop([-l31, -l28, l25, l29])
    s14 = geom.addSurfaceFilling([ll14])

    ll15 = geom.addCurveLoop([-l32, -l29, l26, l27])
    s15 = geom.addSurfaceFilling([ll15])

    ll16 = geom.addCurveLoop([l30, l31, l32])
    s16 = geom.addPlaneSurface([ll16])

    sl3 = geom.addSurfaceLoop([s12, s13, s14, s15, s16])
    vol3 = geom.addVolume([sl3],tag=10000)

    
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
    p_lateral = []
    

    for volume in a:
        value = gmsh.model.getAdjacencies(3, volume[1])
        print(value)
        p_lateral.append(value[1][2])
        p_lateral.append(value[1][3])

        
    gmsh.model.addPhysicalGroup(2, p_lateral, tag=7)


    
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

    injection('MeshDir/injection', fltk=True, l_ec=0.041)