from dolfin import *
import numpy as np
from ufl.operators import Dn # directional derivative

class Class1:

    def __init__(self, geometry, boundary_conditions, p_dir, p_adj, v=1.0):

        self.v = v

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

        p_dir = self.p_dir
        p_adj = self.p_adj

        return - Dn(p_adj) * Dn(p_dir)

    #
    #

    def integral_2(self):

        p_dir = self.p_dir
        p_adj = self.p_adj

        return p_adj * Dn(Dn(p_dir))

    #
    #

    def common_term(self):

        v     = self.v
        p_dir = self.p_dir
        p_adj = self.p_adj
        n     = self.n

        return ( p_adj * ( grad(p_dir) - Dn(p_dir)*n ) )

    #
    #

    def integral_3(self):

        common_term = self.common_term()
        curvature   = self.curvature       
        n           = self.n

        return - curvature * dot(common_term, n) + div(common_term) - dot(dot(grad(common_term), n), n)

    #
    #

    def get_Dirichlet(self):

        r1 = self.integral_1()

        return r1

    #
        
    def get_Neumann(self):

        r2 = self.integral_2()
        r3 = self.integral_3()

        return r2 + r3

    #

    def get_Robin(self):

        r1 = self.integral_1()
        r2 = self.integral_2()
        r3 = self.integral_3()

        return r1 + r2 + r3

    #
    #

    def __call__(self, i):

        v = self.v
        n = self.n
        ds = self.ds

        self.geometry.indices_of_boundary_points(i)

        control_pts = self.geometry.control_pts  # ndarray

        self.curvature = self.geometry.get_curvature_field(i)

        shape_derivatives = np.zeros((len(control_pts),2))

        for j in range(len(control_pts)):

            Vx, Vy = self.geometry.get_displacement_field(i,j) 
            # V = interpolate(V, self.V)       

            if "Dirichlet" in self.boundary_conditions[i]:
                Dir_real = self.get_Dirichlet()
                shape_derivatives[j][0] = assemble( dot(Vx, n) * v *  Dir_real * ds(i) )
                shape_derivatives[j][1] = assemble( dot(Vy, n) * v *  Dir_real * ds(i) )

            elif "Neumann" in self.boundary_conditions[i]:
                Neu_real = self.get_Neumann()
                shape_derivatives[j][0]  = assemble( dot(Vx, n) * v * Neu_real * ds(i) )
                shape_derivatives[j][1]  = assemble( dot(Vy, n) * v * Neu_real * ds(i) )
            elif "Robin" in self.boundary_conditions[i]:
                Rob_real = self.get_Robin()
                shape_derivatives[j][0] = assemble( dot(Vx, n) * v *  Rob_real * ds(i) )
                shape_derivatives[j][1] = assemble( dot(Vy, n) * v *  Rob_real * ds(i) )

        return shape_derivatives
