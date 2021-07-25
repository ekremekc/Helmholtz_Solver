from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

from itertools import chain

from geometry_pkg_ekrem import gmsh_utils
from geometry_pkg_ekrem import mshr_utils
from geometry_pkg_ekrem import bsplines

class Geometry:

    def __init__(self, f, visualize = False):
        """
        Parameters
        ----------
        f : function
            the geometry function which can be callable from the geometry module 
            e.g. module_tswc.f()
        visualize : boolean, optional
            Visualizes the generated mesh by using gmsh python API. The default is True.

        """
        

        self.pts, self.ll, self.physical_ll, self.ctrl_pts, self.lcar = f()

        # introduced in bsplines module
        self.bspline = None
        self.curve = None
        # Bspline- point list
        self.bspline_loops = None
        self.bspl_pts = None
        # introduced in points_lines_and_physical_lines method
        self.points = None
        self.lines = None
        self.physical_lines = None
        self.lcar1 = None

        # introduced in __call__ method
        self.mesh = None
        self.boundaries = None

        # introduced indices_of_boundary_points method
        self.control_pts = None
        self.indices = None
        self.boundary_pts = None
        
        self.visualize = visualize
        
    def create_bspline_instance_and_get_points(self):
        if type(self.ctrl_pts)==list:
            
            self.bspline = bsplines.BSpline(self.ctrl_pts)
            self.curve = self.bspline.curve
            self.bspl_pts = self.bspline.pts
            
        elif type(self.ctrl_pts) == dict:
            self.bspl_pts = dict()
            self.bspline = dict()
            for key in self.ctrl_pts:
                self.bspline[key] = bsplines.BSpline(self.ctrl_pts[key])
                self.bspl_pts[key] = self.bspline[key].pts
        else:
            raise ValueError("Check the type of control points!")
    
    def update_control_points(self, t4):

        t4[0] = 0.
        t4[-1] = 0.
        u4 = [item.tolist() for item in t4]

        for p, q in zip(self.ctrl_pts, u4):
            p[0] += q[0]
            p[1] += q[1]
    def move_edge(self, increment):
        """
        Moves top edge for taylor test

        Returns
        -------
        None.

        """
        # elementary points
        y=.0235
        
        y_top = y+increment
        
        p0 = [0., - y]
        p1 = [1., - y]
        p2 = [1., y_top]
        p3 = [0., y_top]
    
        # p0 = [0., -0.125]
        # p1 = [1, -0.125]
        # p2 = [1, .125]
        # p3 = [0., .125]
    
        self.pts  = [p0, p1, p2, p3]
    
    
        # elementary lines
        
        l0 = [3, 0]  # inlet 
        l1 = [1, 2]  # outlet
    
        """ Order of the l0, l1, l2.. is not important
        """
        self.ll = [l0, l1]
        
        # physical_lines
    
        self.physical_ll = [[0], [1]]
    
        
        def g(p0, p1, n):
            # n is the number of intervals
            pts = []
            for i in range(n + 1):
                x = p0[0] + i*(p1[0] - p0[0])/n
                y = p0[1] + i*(p1[1] - p0[1])/n
                pts.append([x, y])
            return pts
        # control points
        pts3 = g(p0, p1, 10)
        pts4 = g(p2, p3, 10)
    
        self.ctrl_pts = {3: pts3, 4: pts4}
        # ctrl_pts = {3: pts4, 4: pts3} # Also works
    
        self.lcar = 9.4e-3
    
        
        
            
    def update_control_point(self, i, j, component, increment):
        """
        Updates control point's element

        Parameters
        ----------
        i : integer
            Physical Tag of the edge
        j : integer
            jth control point along the edge i
        component : int
            0 for x-component, 1 for y-component
        increment : float
            change in component value

        Returns
        -------
        None.

        """
        self.ctrl_pts[i][j][component] += increment
        
    def points_lines_and_physical_lines(self):
        """
        This function determines lines, pyhsical lines and points which are 
        going to be used to build mesh for given geometry in function f()
        
        For if the Bspline is defined as consistent from one end to another end,
        the list condition is used.
        
        If the BSplines defined as dictionary which has the keys for the BSpline edge,
        the dict condition is used. It is important to note that this condition is quite complicated.
        But works for L-shaped domain and Rijke tube. But the same result with the module_2 geometry can be
        achieved by module_1_b

        """
        
        if type(self.bspl_pts)==list:
            self.bspl_pts = self.bspl_pts[1:-1]

            n_pts = len(self.pts)
            n_bspl_pts = len(self.bspl_pts)
    
            bspl_ll = []
            bspl_ll.append([n_pts - 1, n_pts])
            for i in range(n_pts, n_pts + n_bspl_pts - 1):
                bspl_ll.append([i, i + 1])
            bspl_ll.append([n_pts + n_bspl_pts - 1, 0])
    
            n_ll = len(self.ll)
            n_bspl_ll = len(bspl_ll)
    
            bspl_physical_ll = [n_ll + i for i in range(n_bspl_ll)]
            # split into two lines for readability
            bspl_physical_ll = [bspl_physical_ll]
    
            self.points = self.pts + self.bspl_pts
            self.lines = self.ll + bspl_ll
            self.physical_lines = self.physical_ll + bspl_physical_ll
    
            # lazy (temporary) solution
            self.lcar1 = [self.lcar for _ in self.points]  # _ is a placeholder
            
        elif type(self.bspl_pts) == dict:
            
            for key in self.bspl_pts:
                self.bspl_pts[key] = self.bspl_pts[key][1:-1]

            ###
    
            n_pts = len(self.pts)
            n_bspl_pts = dict()
            for key in self.bspl_pts:
                n_bspl_pts[key] = len(self.bspl_pts[key])
            
            # Required for assigning lines for each bspline
            self.bspline_loops = dict()
            for key in self.ctrl_pts:
                    self.bspline_loops[key] = [self.pts.index(self.ctrl_pts[key][0]),
                                               self.pts.index(self.ctrl_pts[key][-1])]
            
            # Assigning lines     
            bspl_ll = dict()
            for key in self.bspl_pts:
                bspl_ll[key] = []
                bspl_ll[key].append([self.bspline_loops[key][0], n_pts])
                for i in range(n_pts, n_pts + n_bspl_pts[key] - 1):
                    bspl_ll[key].append([i, i + 1])
                bspl_ll[key].append([n_pts + n_bspl_pts[key] - 1, self.bspline_loops[key][1]])
                n_pts+=n_bspl_pts[key]
    
            # Assigning physical lines
            n_ll = len(self.ll)
            bspl_physical_ll = dict()
            for key in self.bspl_pts:
                bspl_physical_ll[key] = [n_ll + i for i in range(len(bspl_ll[key]))]
                bspl_physical_ll[key] = [bspl_physical_ll[key]]
                n_ll += len(bspl_ll[key])
    
            ###
    
            self.points = self.pts 
            self.lines = self.ll 
            self.physical_lines = self.physical_ll 
            for key in self.bspl_pts:
                self.points += self.bspl_pts[key]
                self.lines  += bspl_ll[key]
                self.physical_lines += bspl_physical_ll[key]
            # lazy (temporary) solution
            self.lcar1 = [self.lcar for _ in self.points]  # _ is a placeholder
        
        else:
            raise ValueError("Check the physical lines!")

    def create_mesh_geom(self, name):
        """Calls GMSH API to generate geometry
        """

        gmsh_utils.create_geom(self.points, self.lcar1, self.lines, self.physical_lines, file=name, visualization=self.visualize)

    def __call__(self, name):
        self.create_bspline_instance_and_get_points()
        self.points_lines_and_physical_lines()
        self.create_mesh_geom(name)
        self.mesh, self.boundaries = mshr_utils.mesh_and_boundaries(name)

        return self.mesh, self.boundaries

    def indices_of_boundary_points(self, i):
        """
        Parameters
        ----------
        i : Tag of the physical curve (edge in this case)
        Returns
        -------
        indices of the mesh points and their coordinates along the 
        corresponding physical edge i
        """
        if type(self.bspline)== dict:
            self.control_pts = self.bspline[i].P
        else:
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

    def get_curvature_field(self, i):
        """
        Parameters
        ----------
        i : Tag of the physical curve (edge in this case)
        Returns
        -------
        curvature : curvature field of the edge(s) 
                    based on the corresponding boundary 
                    i's points    
        """
        self.indices_of_boundary_points(i)
        V = FunctionSpace(self.mesh, "CG", 1)
        vertex_to_dof_V = vertex_to_dof_map(V)
        dofs_V = vertex_to_dof_V[self.indices]

        a = np.array([self.bspline[i].get_curvature_from_point(p) for p in self.boundary_pts])

        curvature = Function(V)
        curvature.vector()[dofs_V] = a

        return curvature

    def get_displacement_field(self, i, j):
        """
        computes displacement field's components on x and y direction'        
        
        Parameters
        ----------
        i : Tag of the physical curve (edge in this case)
        j : jth control point 
        Returns
        -------
        V_x : Displacement field along x-direction
        V_y : Displacement field along y-direction
        
        """
        self.indices_of_boundary_points(i)
        W = VectorFunctionSpace(self.mesh, "CG", 1)

        vertex_to_dof_W = vertex_to_dof_map(W)
        vertex_to_dof_W = vertex_to_dof_W.reshape((-1, self.mesh.geometry().dim()))

        dofs_W = vertex_to_dof_W[self.indices]
        dofs_W = dofs_W.reshape(-1)

        b = np.zeros((len(self.boundary_pts), 2))
        c = np.zeros((len(self.boundary_pts), 2))

        for (k, p) in enumerate(self.boundary_pts):
            b[k, 0] = self.bspline[i].get_displacement_from_point(j, p)
            c[k, 1] = self.bspline[i].get_displacement_from_point(j, p)

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

    # You need to make sure that you start from the end of the control point
    # which is end of the bspline, then finish the line loop line by line until 
    # reaching the start of the bspline that you want to obtain
    
    """
     3__________ 2
     |          |   if you want to obtain bspline from 0 to 1,
     |          |   start looping from point 1 like;
     |__________|   pts  = [p1, p2, p3, p0]
     0           1  
    """
    
    
    import module_1
    
    inst = Geometry(module_1.f)

    if not os.path.exists("mesh_dir"):
        os.mkdir("mesh_dir")

    mesh, boundaries = inst("mesh_dir/rijke")

    # plot(mesh)
    # plt.show()

    # from helmholtz_pkg.mshr import MeshXDMF
    #
    # geometry = MeshXDMF("my_dir/mesh_1", write_xdmf_file=True)
    # geometry()
    #
    # plot(geometry.mesh)
    # plt.show()



