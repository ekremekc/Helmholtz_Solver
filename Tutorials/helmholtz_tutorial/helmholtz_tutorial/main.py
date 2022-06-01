# import dolfin as dolf
import numpy as np
import matplotlib.pyplot as plt

from helmholtz_tutorial import params
from helmholtz_tutorial.passive_flame import *
from helmholtz_tutorial.eigensolvers import *
from helmholtz_tutorial.eigenvectors import *

def print_matrix(M):
    print('Size of matrix is: ',(M.getSize()))
    for i in range(M.getSize()[0]):
        for j in range(M.getSize()[1]):
            print(M.getValue(i,j), '\t\t' ,end = '')
        print('\n')
    print('---------------o------------o---------------')

def mshr(nx):
    """
    :param nx: number of elements
    :return: mesh and boundaries
    """

    class Left(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], 0)

    class Right(dolf.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and dolf.near(x[0], 1)

    left = Left()
    right = Right()

    mesh = dolf.IntervalMesh(nx, 0, 1)

    boundaries = dolf.MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)

    left.mark(boundaries, 1)
    right.mark(boundaries, 2)

    return mesh, boundaries


def run():
    """

    :param A: matrix~petsc4py.PETSc.Mat
    :param B: matrix~petsc4py.PETSc.Mat
    :param C: matrix~petsc4py.PETSc.Mat
    :param vr: vector~petsc4py.PETSc.Vec
    :param vi: vector~petsc4py.PETSc.Vec
    :param E: solver object~slepc4py.SLEPc.EPS

    Helmholtz equation with Dirichlet and/or Neumann boundary conditions
    (The eigensolver uses a shift-and-invert spectral transformation)
    
    :param mesh: mesh object~dolfin.cpp.generation.UnitIntervalMesh
    :param boundaries : boundary object~dolfin.cpp.mesh.MeshFunctionSizet
    
    :return:
    """

    degree = 2

    mesh, boundaries = mshr(100)

    boundary_conditions = {1: 'Dirichlet', 2: 'Dirichlet'}

    # A, B, C = passive_flame(mesh, boundaries, boundary_conditions,
    #                         c=dolf.Constant(1.0), degree=degree) #c=dolf.Constant(1.0)
    A, B, C = passive_flame(mesh, boundaries, boundary_conditions,
                            c=params.c, degree=degree) 
    
    # print_matrix(A)
    # print_matrix(C)
    
    E = eps_solver(A, C, np.pi**2, 2, print_results=True)

    vr, vi = A.createVecs()
    k = E.getEigenpair(0, vr, vi)
    
    omega = np.sqrt(k)

    p = normalize_1(mesh, vr, degree=degree)

    dolf.plot(-p)
    plt.xlabel('x')
    plt.ylabel('p')
    plt.show()

def run_mixed():
    """
    Helmholtz equation with Dirichlet and/or Neumann boundary conditions
    on a mixed function space (real representation of complex numbers)
    (The eigensolver uses a shift-and-invert spectral transformation)
    :return:
    """

    degree = 1

    mesh, boundaries = mshr(400)

    boundary_conditions = {1: 'Dirichlet', 2: 'Dirichlet'}

    A, B, C = passive_flame_mixed(mesh, boundaries, boundary_conditions,
                                  c=params.c, degree=degree)
    
    # print_matrix(A)
    # print_matrix(C)
    
    E = eps_solver(A, C, np.pi**2, 2, print_results=True)

    vr, vi = A.createVecs()
    k = E.getEigenpair(1, vr, vi)

    omega = np.sqrt(k)

    p_r, p_i = normalize_2(mesh, vr, vi, degree=degree)

    # dolf.plot(p_r)
    # plt.xlabel('x')
    # plt.ylabel('p_r')
    # plt.show()
    
    pl, ax = plt.subplots(); fig = plt.gcf(); fig.set_size_inches(10, 3.6)
    plt.subplot(1, 2, 1); p1 = dolf.plot(p_r); plt.xlabel('x');plt.title('Real Part')
    plt.subplot(1, 2, 2); p2 = dolf.plot(p_i); plt.xlabel('x');plt.title('Imaginary Part')
    plt.tight_layout()
    plt.show()


def run_mixed_Robin():
    """
    Helmholtz equation with Robin boundary conditions
    on a mixed function space (real representation of complex numbers)
    The eigenvalue problem is quadratic
    (The eigensolver uses a shift-and-invert spectral transformation)
    :return:
    """

    degree = 1

    mesh, boundaries = mshr(400)

    # boundary_conditions = {1: {'Robin': params.Y_in},
    #                         2: {'Robin': params.Y_out}}
    
    boundary_conditions = {1: {'Neumann'},
                            2: {'Robin': params.Y_out}}

    A, B, C = passive_flame_mixed(mesh, boundaries, boundary_conditions,
                                  c=dolf.Constant(1.0), degree=degree)
    
    # A, B, C = passive_flame_mixed(mesh, boundaries, boundary_conditions,
    #                               c=params.c, degree=degree)
    
    
    E = pep_solver(A, B, C, np.pi, 2, print_results=True)

    vr, vi = A.createVecs()
    k = E.getEigenpair(1, vr, vi)

       
    # p_r, p_i = normalize_3(mesh, vr, vi, degree=degree)
    
    # pl, ax = plt.subplots(); fig = plt.gcf(); fig.set_size_inches(10, 3.6)
    # plt.subplot(1, 2, 1); p1 = dolf.plot(p_r); plt.xlabel('x');plt.title('Real Part')
    # plt.subplot(1, 2, 2); p2 = dolf.plot(p_i); plt.xlabel('x');plt.title('Imaginary Part')
    # plt.tight_layout()
    # plt.show()


if __name__ == '__main__':

    # run()
    # run_mixed()
    run_mixed_Robin()
 