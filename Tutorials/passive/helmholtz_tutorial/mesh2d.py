from helmholtz_tutorial.main import *
import matplotlib.pyplot as plt

c_0 = 450
L = 0.4
h = 0.1

"-----------o------- NUMERICAL RESULTS -------o-----------"

nx,ny,N_num = 40,10,40
Z_right_b = np.linspace(-10j, 10j, N_num)
Z_right_a = np.linspace(-10, 10, N_num)

def mshr_(L, h, nx, ny):
    """
    :param nx: number of elements
    :return: mesh and boundaries
    """

    class Left(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], 0)

    class Right(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], L)

    class Bottom(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], 0)

    class Top(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], h)

    left = Left()
    right = Right()
    bottom = Bottom()
    top = Top()

    mesh = dolf.RectangleMesh(dolf.Point(0.0, 0.0), dolf.Point(L, h), nx, ny)
    # mesh = dolf.IntervalMesh(nx, 0, 0.4)

    boundaries = dolf.MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)

    left.mark(boundaries, 1)
    top.mark(boundaries, 2)
    right.mark(boundaries, 3)
    bottom.mark(boundaries, 4)

    return mesh, boundaries

mesh, dum =mshr_(L, h, nx, ny) 

dolf.plot(mesh)
