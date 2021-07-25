from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

from itertools import chain

from geometry_pkg import bspline_module
from geometry_pkg import pygmsh_utils
from geometry_pkg import unstructured_mesh


class Geometry:
    def __init__(self, f):

        self.pts, self.ll, self.physical_ll, self.ctrl_pts, self.lcar = f()

        # introduced in create_bspline_instance_and_get_points method
        self.bsplines = dict()
        self.bsplines_pts = dict()

        # introduced in points_lines_and_physical_lines method
        self.points = None
        self.lines = None
        self.physical_lines = None
        self.lcar1 = None

        # introduced in create_gmsh_geom method
        self.geom = None

        # introduced in __call__ method
        self.mesh = None
        self.boundaries = None

        # introduced method_1 method
        self.bspline = None
        self.control_pts = None
        self.indices = None
        self.boundary_pts = None

    def update_control_points(self, t3, t4):

        t3[0] = 0.
        t3[-1] = 0.
        u3 = [item.tolist() for item in t3]

        for p, q in zip(self.ctrl_pts[3], u3):
            p[0] += q[0]
            p[1] += q[1]

        t4[0] = 0.
        t4[-1] = 0.
        u4 = [item.tolist() for item in t4]

        for p, q in zip(self.ctrl_pts[4], u4):
            p[0] += q[0]
            p[1] += q[1]

    def create_bspline_instance_and_get_points(self):

        for i in (3, 4):
            self.bsplines[i] = bspline_module.BSpline(self.ctrl_pts[i], lcar=self.lcar, p=3)
            self.bsplines[i].create_curve()
            self.bsplines_pts[i] = self.bsplines[i].get_pts_as_list()

    def points_lines_and_physical_lines(self):

        self.bsplines_pts[3] = self.bsplines_pts[3][1:-1]
        self.bsplines_pts[4] = self.bsplines_pts[4][1:-1]

        ###

        n_pts = len(self.pts)
        n_bspl_pts_3 = len(self.bsplines_pts[3])
        n_bspl_pts_4 = len(self.bsplines_pts[4])

        bspl_ll_3 = []
        bspl_ll_3.append([0, n_pts])
        for i in range(n_pts, n_pts + n_bspl_pts_3 - 1):
            bspl_ll_3.append([i, i + 1])
        bspl_ll_3.append([n_pts + n_bspl_pts_3 - 1, 1])

        bspl_ll_4 = []
        bspl_ll_4.append([2, n_pts + n_bspl_pts_3])
        for i in range(n_pts + n_bspl_pts_3, n_pts + n_bspl_pts_3 + n_bspl_pts_4 - 1):
            bspl_ll_4.append([i, i + 1])
        bspl_ll_4.append([n_pts + n_bspl_pts_3 + n_bspl_pts_4 - 1, 3])

        ###

        n_ll = len(self.ll)

        n_bspl_ll_3 = len(bspl_ll_3)
        bspl_physical_ll_3 = [n_ll + i for i in range(n_bspl_ll_3)]
        bspl_physical_ll_3 = [bspl_physical_ll_3]

        n_bspl_ll_4 = len(bspl_ll_4)
        bspl_physical_ll_4 = [n_ll + n_bspl_ll_3 + i for i in range(n_bspl_ll_4)]
        bspl_physical_ll_4 = [bspl_physical_ll_4]

        ###

        self.points = self.pts + self.bsplines_pts[3] + self.bsplines_pts[4]
        self.lines = self.ll + bspl_ll_3 + bspl_ll_4
        self.physical_lines = self.physical_ll + bspl_physical_ll_3 + bspl_physical_ll_4

        # lazy (temporary) solution
        self.lcar1 = [self.lcar for _ in self.points]  # _ is a placeholder

    def create_gmsh_geom(self):

        self.geom = pygmsh_utils.create_geom(self.points, self.lcar1, self.lines, self.physical_lines)
        # print('self.geom', type(self.geom))
        # <class 'pygmsh.built_in.geometry.Geometry'>

    def __call__(self, filename):

        self.create_bspline_instance_and_get_points()
        self.points_lines_and_physical_lines()
        self.create_gmsh_geom()

        self.mesh, self.boundaries = unstructured_mesh.Unstructured(self.geom, filename)

        return self.mesh, self.boundaries

    def method_1(self, i):

        self.bspline = self.bsplines[i]
        self.control_pts = self.bspline.P

        indices = []
        k = 0  # counter
        for f in facets(self.mesh):
            if self.boundaries[f] == i:
                indices.append([])
                for v in vertices(f):
                    indices[k].append(v.index())
                k += 1

        # split into separate lines for readability
        self.indices = set(chain.from_iterable(indices))
        self.indices = list(self.indices)
        coordinates = self.mesh.coordinates()
        self.boundary_pts = coordinates[self.indices]

    def get_curvature_field(self):

        V = FunctionSpace(self.mesh, "CG", 1)
        vertex_to_dof_V = vertex_to_dof_map(V)
        dofs_V = vertex_to_dof_V[self.indices]

        a = np.array([self.bspline.get_curvature_from_point(p) for p in self.boundary_pts])

        curvature = Function(V)
        curvature.vector()[dofs_V] = a

        return curvature

    def get_displacement_field(self, j):

        W = VectorFunctionSpace(self.mesh, "CG", 1)

        vertex_to_dof_W = vertex_to_dof_map(W)
        vertex_to_dof_W = vertex_to_dof_W.reshape((-1, self.mesh.geometry().dim()))

        dofs_W = vertex_to_dof_W[self.indices]
        dofs_W = dofs_W.reshape(-1)

        b = np.zeros((len(self.boundary_pts), 2))
        c = np.zeros((len(self.boundary_pts), 2))

        for (i, p) in enumerate(self.boundary_pts):
            b[i, 0] = self.bspline.get_displacement_from_point(j, p)
            c[i, 1] = self.bspline.get_displacement_from_point(j, p)

        d = b.reshape(-1)
        e = c.reshape(-1)

        V_x = Function(W)
        V_x.vector()[dofs_W] = d

        V_y = Function(W)
        V_y.vector()[dofs_W] = e

        return V_x, V_y


if __name__ == '__main__':

    import os
    # import sys
    # sys.path.insert(1, '/Users/stefanofalco/Desktop/FEniCS_dev/Rijke_opt')

    import module_5

    inst = Geometry(module_5.f)

    if not os.path.exists("my_dir"):
        os.mkdir("my_dir")

    mesh, boundaries = inst("my_dir/mesh_1")

    plot(mesh)
    plt.show()

    # from helmholtz_pkg.mshr import MeshXDMF
    #
    # geometry = MeshXDMF("my_dir/mesh_1", write_xdmf_file=True)
    # geometry()
    #
    # plot(geometry.mesh)
    # plt.show()



