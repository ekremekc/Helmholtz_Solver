from dolfin import *
import numpy as np
from ufl.operators import Dn  # directional derivative

#
#

class Class1:

    def __init__(self, geometry, boundary_conditions, p_dir, p_adj, c):
        """
        Shape derivative class for geometry which has Robin BC

        Parameters
        ----------
        geometry : Geometry
            Geometry object produced by this module
        boundary_conditions : dict
            Boundary conditions dictionary
        p_dir : Coefficient
            Direct pressure eigenfunction (not splitted)
        p_adj : Coefficient
            Adjoint pressure eigenfunction (not splitted)
        c : Coefficient
            Speed of sound, it should be interpolated on defined function space!

        Returns
        -------
        None.

        """

        self.c = c
        
        self.p_dir = p_dir
        self.p_adj = p_adj

        self.boundary_conditions = boundary_conditions

        self.geometry = geometry

        self.mesh = self.geometry.mesh
        self.boundaries = self.geometry.boundaries
        self.n = FacetNormal(self.mesh)
        self.ds = Measure("ds", domain=self.mesh, subdomain_data=self.boundaries)

    #
    #

    def integral_1(self):

        c = self.c

        (p_dir_real, p_dir_imag) = self.p_dir.split()
        (p_adj_real, p_adj_imag) = self.p_adj.split()

        real_part = - c**2 * (   Dn(p_adj_real) * Dn(p_dir_real)
                            + Dn(p_adj_imag) * Dn(p_dir_imag) )

        imag_part = - c**2 * (   Dn(p_adj_real) * Dn(p_dir_imag)
                            - Dn(p_adj_imag) * Dn(p_dir_real) )

        return real_part, imag_part

    #
    #

    def integral_2(self):

        c = self.c

        (p_dir_real, p_dir_imag) = self.p_dir.split()
        (p_adj_real, p_adj_imag) = self.p_adj.split()

        real_part = c**2 * (   p_adj_real * Dn(Dn(p_dir_real))
                          + p_adj_imag * Dn(Dn(p_dir_imag)) )

        imag_part = c**2 * (   p_adj_real * Dn(Dn(p_dir_imag))
                          - p_adj_imag * Dn(Dn(p_dir_real)) )

        return real_part, imag_part

    #
    #

    def common_term(self):

        c = self.c

        (p_dir_real, p_dir_imag) = self.p_dir.split()
        (p_adj_real, p_adj_imag) = self.p_adj.split()

        n = self.n

        real_part = c**2 * (   p_adj_real * ( grad(p_dir_real) - Dn(p_dir_real)*n )
                          + p_adj_imag * ( grad(p_dir_imag) - Dn(p_dir_imag)*n ) )

        imag_part = c**2 * (   p_adj_real * ( grad(p_dir_imag) - Dn(p_dir_imag)*n )
                          - p_adj_imag * ( grad(p_dir_real) - Dn(p_dir_real)*n ) )

        return real_part, imag_part

    #

    def integral_3(self):
        
        (real_part, imag_part) = self.common_term()

        curvature = self.curvature

        n = self.n

        real_part = curvature * dot(real_part, n) + div(real_part) - dot(dot(grad(real_part), n), n)
        imag_part = curvature * dot(imag_part, n) + div(imag_part) - dot(dot(grad(imag_part), n), n)

        return real_part, imag_part

    #
    #

    def get_Dirichlet(self):

        (r1, c1) = self.integral_1()

        Dir_real = r1
        Dir_imag = c1

        return Dir_real, Dir_imag

    #
        
    def get_Neumann(self):

        (r2, c2) = self.integral_2()
        (r3, c3) = self.integral_3()

        Neu_real = r2 + r3
        Neu_imag = c2 + c3

        return Neu_real, Neu_imag

    #

    def get_Robin(self):

        (r1, c1) = self.integral_1()
        (r2, c2) = self.integral_2()
        (r3, c3) = self.integral_3()

        Rob_real = r1 + r2 + r3
        Rob_imag = c1 + c2 + c3

        return Rob_real, Rob_imag

    #
    #

    def __call__(self, i):

        """
        i
        3 : bottom
        4 : top
        """

        c = self.c
        n = self.n
        ds = self.ds

        self.geometry.indices_of_boundary_points(i)

        ctrl_pts = self.geometry.control_pts  # ndarray

        self.curvature = self.geometry.get_curvature_field(i)

        # def foo(V, expr):  # i is a global variable
        #     mylist = []
        #     for V_ in V:
        #         for expr_ in expr:
        #             mylist.append(assemble(dot(V_, n) * expr_ * ds(i)))
            # return mylist

        shape_derivatives = np.zeros((len(ctrl_pts), 2), dtype=complex)

        # for j in range(len(ctrl_pts)):

        #     V = self.geometry.get_displacement_field(i, j)

        #     if "Dirichlet" in self.boundary_conditions[i]:
        #         expr = self.get_Dirichlet()
        #     elif "Neumann" in self.boundary_conditions[i]:
        #         expr = self.get_Neumann()
        #     elif "Robin" in self.boundary_conditions[i]:
        #         expr = self.get_Robin()

        #     mylist = foo(V, expr)

        #     shape_derivatives[(j, 0)] = np.complex(mylist[0], mylist[1])
        #     shape_derivatives[(j, 1)] = np.complex(mylist[2], mylist[3])

        for j in range(len(ctrl_pts)):

            V_x, V_y = self.geometry.get_displacement_field(i, j)

            if "Dirichlet" in self.boundary_conditions[i]:

                (Dir_real, Dir_imag) = self.get_Dirichlet()

                real_x = assemble( dot(V_x, n) * Dir_real * ds(i) )
                imag_x = assemble( dot(V_x, n) * Dir_imag * ds(i) )

                real_y = assemble( dot(V_y, n) * Dir_real * ds(i) )
                imag_y = assemble( dot(V_y, n) * Dir_imag * ds(i) )

            elif "Neumann" in self.boundary_conditions[i]:

                (Neu_real, Neu_imag) = self.get_Neumann()

                real_x = assemble( dot(V_x, n) * Neu_real * ds(i) )
                imag_x = assemble( dot(V_x, n) * Neu_imag * ds(i) )

                real_y = assemble( dot(V_y, n) * Neu_real * ds(i) )
                imag_y = assemble( dot(V_y, n) * Neu_imag * ds(i) )

            elif "Robin" in self.boundary_conditions[i]:

                (Rob_real, Rob_imag) = self.get_Robin()

                real_x = assemble( dot(V_x, n) * Rob_real * ds(i) )
                imag_x = assemble( dot(V_x, n) * Rob_imag * ds(i) )

                real_y = assemble( dot(V_y, n) * Rob_real * ds(i) )
                imag_y = assemble( dot(V_y, n) * Rob_imag * ds(i) )

            shape_derivatives[j, 0] = np.complex(real_x, imag_x)
            shape_derivatives[j, 1] = np.complex(real_y, imag_y)

        return shape_derivatives

#
#

# example
# instance_1 = shape_derivatives_complex.Class1(geometry, boundary_conditions, p_dir, p_adj, v=v)
# shape_der_3 = instance_1(3)
# shape_der_4 = instance_1(4)
