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


def chamber(file, fltk=False, **kwargs):

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(__name__)
    
    # geom = gmsh.model.geo

    add_elementary_entities(**kwargs)
    
    apply_symmetry()
    print(gmsh.model.getEntities(dim=2))
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
    # FLAME
    p22 = geom.addPoint(R_mid, 0., 0., lc_2, tag=20000)
    p26 = geom.addPoint(R_mid + h_f, 0., 0., lc_2)
    p27 = geom.addPoint(R_mid, h_f, 0., lc_2)
    p28 = geom.addPoint(R_mid - h_f, 0., 0., lc_2)

    l35 = geom.addCircleArc(p26, p22, p27, tag=20000)
    l36 = geom.addCircleArc(p27, p22, p28)


    p29 = geom.addPoint(R_mid, 0., 0. + l_f, lc_2)

    p30 = geom.addPoint(R_mid + h_f, 0., 0. + l_f, lc_2)
    p31 = geom.addPoint(R_mid, h_f, 0. + l_f, lc_2)
    p32 = geom.addPoint(R_mid - h_f, 0., 0. + l_f, lc_2)

    l37 = geom.addLine(p28, p32)
    l38 = geom.addLine(p26, p30)
    l39 = geom.addLine(p27, p31)

    l40 = geom.addLine(p32, p30)
    l41 = geom.addCircleArc(p30, p29, p31)
    l42 = geom.addCircleArc(p31, p29, p32)


    ll19 = geom.addCurveLoop([-l41, -l38, l35, l39])
    s19 = geom.addSurfaceFilling([ll19], tag=20000)

    ll20 = geom.addCurveLoop([-l42, -l39, l36, l37])
    s20 = geom.addSurfaceFilling([ll20])

    ll21 = geom.addCurveLoop([l40, l41, l42])
    s21 = geom.addPlaneSurface([ll21])


    # ________________________________________________________________________________
    # COMBUSTION CHAMBER

    p33 = geom.addPoint(0., 0., 0., lc_1)

    p34 = geom.addPoint(R_in_cc, 0., 0., lc_1)
    p35 = geom.addPoint(R_out_cc, 0., 0., lc_1)
    p36 = geom.addPoint(R_in_cc * cos(theta), R_in_cc * sin(theta), 0., lc_1)
    p37 = geom.addPoint(R_out_cc * cos(theta), R_out_cc * sin(theta), 0., lc_1)

    l43 = geom.addLine(p34, p28)
    l44 = geom.addLine(p26, p35)

    l45 = geom.addCircleArc(p35, p33, p37)
    l46 = geom.addLine(p37, p36)
    l47 = geom.addCircleArc(p36, p33, p34)

    ll22 = geom.addCurveLoop([-l44, l35, l36, -l43, -l47, -l46, -l45])
    s22 = geom.addPlaneSurface([ll22])

    p38 = geom.addPoint(0., 0., l_cc + l_ec, lc_1)
    p39 = geom.addPoint(R_in_cc, 0., l_cc + l_ec, lc_1)
    p40 = geom.addPoint(R_out_cc, 0., l_cc + l_ec, lc_1)
    p41 = geom.addPoint(R_in_cc * cos(theta), R_in_cc * sin(theta), l_cc + l_ec, lc_1)
    p42 = geom.addPoint(R_out_cc * cos(theta), R_out_cc * sin(theta), l_cc + l_ec, lc_1)

    l48 = geom.addLine(p34, p39)
    l49 = geom.addLine(p35, p40)
    l50 = geom.addLine(p37, p42)
    l51 = geom.addLine(p36, p41)

    l52 = geom.addLine(p39, p40)
    l53 = geom.addCircleArc(p40, p38, p42)
    l54 = geom.addLine(p42, p41)
    l55 = geom.addCircleArc(p41, p38, p39)

    ll23 = geom.addCurveLoop([-l52, -l48, l43, l37, l40, -l38, l44, l49])
    s23 = geom.addPlaneSurface([ll23])

    ll24 = geom.addCurveLoop([-l49, l45, l50, -l53])
    s24 = geom.addSurfaceFilling([ll24])

    ll25 = geom.addCurveLoop([-l54, -l50, l46, l51])
    s25 = geom.addPlaneSurface([ll25])

    ll26 = geom.addCurveLoop([l48, -l55, -l51, l47])
    s26 = geom.addSurfaceFilling([ll26])

    ll27 = geom.addCurveLoop([l52, l53, l54, l55])
    s27 = geom.addPlaneSurface([ll27])

    sl5 = geom.addSurfaceLoop([s19, s20, s21, s22, s23, s24, s25, s26, s27])
    vol5 = geom.addVolume([sl5], tag=20000)


    
def apply_symmetry():

    symmetry = gmsh.model.geo.symmetrize

    gmsh.model.geo.synchronize()
    old_entities = [ent[1] for ent in gmsh.model.getEntities(dim=2)]

    symmetry(gmsh.model.geo.copy(gmsh.model.getEntities(dim=2)), 0, 1, 0, 0)
    symmetry(gmsh.model.geo.copy(gmsh.model.getEntities(dim=3)), 0, 1, 0, 0)

    gmsh.model.geo.synchronize()
    
    new_entities = [ent[1] for ent in gmsh.model.getEntities(dim=2)]
    new_entities = new_entities[len(old_entities):]
    old_entities.pop(4) # exctract intersection plane
    # print(old_entities, new_entities)
    gmsh.model.mesh.setPeriodic(2, new_entities, old_entities, reflection_matrix())
    gmsh.model.geo.synchronize()

def apply_rotation():

    theta = 22.5 / 2
    theta *= pi / 180

    rotate = gmsh.model.geo.rotate

    gmsh.model.geo.synchronize()
    my_surf = gmsh.model.getEntities(dim=2)
    my_surf_old = gmsh.model.getEntities(dim=2)
    
    old_entities = gmsh.model.getEntities(dim=2)  # for periodicity
    print(gmsh.model.getEntities(dim=2))
    old_entities.pop(-3) #Substract intersection plane(bottom of first sample)
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
            old_entities.pop(6) #Extract intersection plane(top of first sample)

        gmsh.model.mesh.setPeriodic(2, new_entities, old_entities, rotation_matrix(angle))
        my_surf_old = my_surf_new  
        gmsh.model.geo.synchronize()

def add_physical_entities():
    a=gmsh.model.getEntities(dim=3)
    cc_rear = []   
    cc_outer = []
    c_inner = []
    cc_front = []
    for volume in a:
        value = gmsh.model.getAdjacencies(3, volume[1])
        # print(value)
        cc_rear.append(value[1][3])
        cc_outer.append(value[1][5])
        c_inner.append(value[1][7])
        cc_front.append(value[1][8])
        
        
    gmsh.model.addPhysicalGroup(2, cc_rear, tag=8)
    gmsh.model.addPhysicalGroup(2, cc_outer, tag=9)
    gmsh.model.addPhysicalGroup(2, c_inner, tag=10)
    gmsh.model.addPhysicalGroup(2, cc_front, tag=11)
    
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

    chamber('MeshDir/chamber', fltk=True, l_ec=0.041)