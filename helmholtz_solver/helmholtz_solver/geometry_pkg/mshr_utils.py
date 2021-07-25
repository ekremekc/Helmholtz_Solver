from geometry_pkg_ekrem import mshr

def mesh_and_boundaries(filename):
    """
    Computes the all the mesh dependencies and boundaries of given geometry by using mshr package

    Parameters
    ----------
    filename : Name of the file which contains requested output formats of mesh
    
    Returns
    -------
    mesh : dolfin.cpp.mesh.Mesh
        generated mesh of the given geometry
    boundaries : dolfin.cpp.mesh.MeshFunctionSizet
        boundary data of the generated mesh

    """
    geometry = mshr.MeshXDMF(filename, write_xdmf_file=True)
    geometry()
    mesh = geometry.mesh
    boundaries = geometry.boundaries
    # subdomains = geometry.subdomains

    return mesh, boundaries

