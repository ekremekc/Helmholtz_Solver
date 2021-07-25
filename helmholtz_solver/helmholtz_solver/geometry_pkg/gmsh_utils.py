import gmsh

def create_geom(points, lcar, lines, physical_lines, file, visualization=False):
    """
    creates meshs using gmsh API
    
    Parameters
    ----------
    points : list
    whole points along the boundary of the geometry
    
    lcar : float
    characteristic length between the points on the boundary
    
    lines : list of lists
    defined lines by function f() in geometry module 
    files(e.g. module_tswc)
    
    physical_lines : list of lists
    defined physical lines by function f() in geometry module files(e.g. module_tswc)
    
    file : string
        filename of the geometry.
    visualization : boolean, optional
        Gmsh visulation tool. The default is False.

    Returns
    -------
    None.

    """
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0) # Prints log into console
    gmsh.model.add(__name__)

    geom = gmsh.model.geo
    model = gmsh.model
    # if refinement:
    #     from pygmsh_file_1 import mesh_refinement
    #     geom = mesh_refinement(geom)

    ############################################################################
    ### elementary entities
    ############################################################################

    # geom_points
    # geom_lines
    # geom_line_loops
    # geom_plane_surfaces

    # points
    geom_points = []

    for p, q in zip(points, lcar):
        geom_points.append(geom.addPoint(p[0], p[1], 0.0, q))

    # lines
    geom_lines = []

    for pts in lines:
        geom_lines.append(geom.addLine(geom_points[pts[0]], geom_points[pts[1]]))

    # line loops
    geom_line_loops = geom.addCurveLoop(geom_lines)

    # surfaces
    geom_plane_surfaces = geom.addPlaneSurface([geom_line_loops])

    ############################################################################
    ### physical groups
    ############################################################################
    gmsh.model.geo.synchronize()
    # lines (boundaries)
    i = 1
    for ll in physical_lines:
        geom.addPhysicalGroup(1, [geom_lines[l] for l in ll], i)
        i += 1

    # surfaces (subdomains)
    geom.addPhysicalGroup(2, [geom_plane_surfaces],i)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    print("Mesh is generated.")
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    
    gmsh.write("{}.msh".format(file))
    print("Mesh file {} is saved.".format(file+".msh"))
    
    gmsh.write("{}.geo_unrolled".format(file))
    print("Mesh file {} is saved.".format(file+".geo_unrolled"))

    if visualization:
        fltk_options()
        gmsh.fltk.run()

    gmsh.finalize()

    
    

def fltk_options():

    # Type of entity label (0: description,
    #                       1: elementary entity tag,
    #                       2: physical group tag)
    gmsh.option.setNumber("Geometry.LabelType", 2)

    gmsh.option.setNumber("Geometry.PointNumbers", 1)
    gmsh.option.setNumber("Geometry.LineNumbers", 2)
    gmsh.option.setNumber("Geometry.SurfaceNumbers", 0)
    gmsh.option.setNumber("Geometry.VolumeNumbers", 0)

    # Mesh coloring(0: by element type, 1: by elementary entity,
    #                                   2: by physical group,
    #                                   3: by mesh partition)
    gmsh.option.setNumber("Mesh.ColorCarousel", 2)

    gmsh.option.setNumber("Mesh.Lines", 0)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 0)

    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    gmsh.option.setNumber("Mesh.VolumeFaces", 0)


if __name__ == '__main__':
    from main_geometry import *
    import module_5
    from dolfin import *
    
    inst = Geometry(module_5.f)
    name = "mesh_dir/mesh_6"
    mesh, boundaries = inst(name)
    
    
    # pts, ll, physical_ll, ctrl_pts, lcar = module_2.f()
    points = inst.points
    lcar = inst.lcar1
    lines = inst.lines
    physical_lines = inst.physical_lines
    
    # create_geom(points, lcar, lines, physical_lines, name, visualization=True)
