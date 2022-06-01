import dolfin as dolf
import numpy as np
import matplotlib.pyplot as plt

import params
from helmholtz_tutorial.passive_flame import *
from helmholtz_tutorial.flame_transfer_function import *
from helmholtz_tutorial.active_flame import *
from helmholtz_tutorial.eigensolvers import *
from helmholtz_tutorial.eigenvectors import *


def mshr(nx):
    """
    :param nx: number of elements
    :return: mesh, boundaries and subdomains
    """

    class Left(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], 0)

    class Right(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], 1)
        
    
    left = Left()
    right = Right()

    mesh = dolf.UnitIntervalMesh(nx)

    boundaries = dolf.MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)

    left.mark(boundaries, 1)
    right.mark(boundaries, 2)
    
    

    # TODO: Add subdomains such that x_f - a_f <= x <= x_f + a_f is marked with 0 and everywhere else is marked with 1.
    #       Ref.: Chapter 4 of the FEniCS Tutorial book.
    # class 
    
    class FlameRegion(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return x[0] >= params.x_f[0][0]-params.a_f+dolf.DOLFIN_EPS and x[0] <= params.x_f[0][0]+params.a_f-dolf.DOLFIN_EPS
        
    subdomains = dolf.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    flame = FlameRegion()
    subdomains.set_all(1)
    flame.mark(subdomains, 0)
    
    
    return mesh, boundaries, subdomains

def run():
    """
    Helmholtz equation with unsteady heat release
    The eigenvalue problem is nonlinear
    :return:
    """

    degree = 1

    mesh, boundaries, subdomains = mshr(400)

    boundary_conditions = {1: {'Robin': params.Y_in},
                           2: {'Robin': params.Y_out}}

    operators = PassiveFlame(mesh, boundaries, boundary_conditions,
                             c=params.c, degree=degree)
    operators.assemble_A()
    operators.assemble_B()
    operators.assemble_C()

    ftf = n_tau(params.n, params.tau)

    D = ActiveFlame(mesh, subdomains, params.x_f, params.x_r, params.rho_in, 1., 1., ftf, degree=degree)
    
    D.assemble_submatrices()

    # TODO: In eigensolvers, write fixed_point_iteration_pep, using pep_solver to solve the quadratic eigenvalue problem
    #       at each iteration. Use 'D.assemble_matrix(omega[k])' to assemble the matrix for the eigenvalue at the k-th
    #       iteration and 'D_Mat = D.matrix' to get the petsc4py.PETSc.Mat matrix.

    E = fixed_point_iteration_pep(operators, D, np.pi, nev=2, i=0)

    omega, p = normalize_eigenvector(mesh, E, 0, degree=degree)
    
    # print(p.vector()[:])
    
    p_r, p_i = p.split(True)
    
    vector = p_r.compute_vertex_values()+1j*p_i.compute_vertex_values()
    new_vector = np.abs(vector)*np.exp(1j*np.angle(vector))
    
    # plt.plot(np.abs(vector))
    # plt.plot(np.angle(vector))
    

    # new_vector *= np.exp(-1j*np.angle(vector)[0]) # Equate the phase to 0 at x = 0

    # plt.plot(np.linspace(0,1,401),new_vector.real, label="real")
    # plt.plot(np.linspace(0,1,401),new_vector.imag, label="imag")
    # plt.legend()
    # plt.show()
    
    #   Calculation of U
    
    # grad_p_real = np.gradient(vector.real)
    # grad_p_imag = np.gradient(vector.imag)
    # grad_p = grad_p_real + 1j*grad_p_imag
    
    # u = -grad_p/(params.gamma*omega)
    # plt.plot(u.real)
    
    pl, ax = plt.subplots(); fig = plt.gcf(); fig.set_size_inches(10, 3.6)
    plt.subplot(1, 2, 1); p1 = dolf.plot(p_r); plt.xlabel('x');plt.title('Real Part')
    plt.subplot(1, 2, 2); p2 = dolf.plot(-p_i); plt.xlabel('x');plt.title('Imaginary Part')
    plt.tight_layout()
    plt.show()
    
    # dolf.plot(p_r)
    # plt.xlabel('x')
    # plt.ylabel('p_r')
    # plt.show()


run()