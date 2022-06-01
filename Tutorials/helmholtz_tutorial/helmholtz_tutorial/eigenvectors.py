import dolfin as dolf
from numpy.linalg import norm
import numpy as np
# from petsc4py import PETSc


def normalize_1(mesh, vr, degree=1):
    """
    real (, with real eigenvalues and eigenvectors)
    :param mesh:
    :param vr: <class 'petsc4py.PETSc.Vec'>
    :param degree:
    :return: <class 'dolfin.function.function.Function'>
    """
    
    # TODO: normalize p such that the L2-norm is 1
    # vr = vr/norm(vr, ord=2)
    print(norm(vr, ord=2))
    # vr = vr/norm(vr, np.inf)
    
    # TODO: normalize p such that the L^infty-Norm is 1
    # vr = vr/norm(vr, np.inf)
    
    V = dolf.FunctionSpace(mesh, 'CG', degree)
    p = dolf.Function(V)
    
    # vr = 
    # vr.getArray() converts a petsc4py.PETSc.Vec into a numpy.ndarray
    # p.vector() converts a dolfin.function.function.Function into a dolfin.cpp.la.PETScVector
    # set (local) values of a dolfin Function
    p.vector().set_local(vr.getArray())
    # p.vector().apply('insert')

    # TODO: normalize p such that the L^infty-Norm is 1
    # TODO: normalize p such that the L2-norm is 1
    
    dx = dolf.Measure('dx', domain=mesh)
    p = p/np.sqrt(dolf.assemble(p*p*dx))    
    print(np.sqrt(dolf.assemble(p*p*dx)))    
    # p = abs(p)/np.sqrt(dolf.assemble(p*p*dx))

    return p


def normalize_2(mesh, vr, vi, degree=1):
    """
    complex, (but) with real eigenvalues and eigenvectors
    :param mesh:
    :param vr: <class 'petsc4py.PETSc.Vec'>
    :param vi: <class 'petsc4py.PETSc.Vec'>
    :param degree:
    :return: <class 'dolfin.function.function.Function'>
    """

    vr_1 = vr.getArray()[0::2]
    vi_1 = vi.getArray()[0::2]
    # print(vr.getArray())
    # print(vi.getArray())
    # shallow copy of vr
    x = vr.copy() 
    x[0::2] = vr_1
    x[1::2] = vi_1

    # TODO: normalize p such that the L2-norm is 1
    # x = x/norm(x, ord=2)
    # x = x/norm(x, np.inf)
    
    # TODO: normalize p such that the L^infty-Norm is 1
    # x = x/norm(x, np.inf)
    
    CG = dolf.FiniteElement('CG', mesh.ufl_cell(), degree)
    W = dolf.FunctionSpace(mesh, dolf.MixedElement([CG, CG]))
    p = dolf.Function(W)
    p.vector().set_local(x.getArray())
    # p.vector().apply('insert')
    (p_r, p_i) = p.split(True)
    
    dx = dolf.Measure('dx', domain=mesh)
    p_r = p_r/np.sqrt(dolf.assemble((p_r*p_r+p_i*p_i)*dx))
    p_i = p_i/np.sqrt(dolf.assemble((p_r*p_r+p_i*p_i)*dx))
    
    return p_r, p_i

# TODO: write normalize_3

def normalize_3(mesh, vr, vi, degree=1):
    """
    complex, with complex eigenvalues and eigenvectors
    """
    v1_r = vr.getArray()[0::2]
    v2_r = vr.getArray()[1::2]
    v1_i = vi.getArray()[0::2]
    v2_i = vi.getArray()[1::2]
    
    real_part = v1_r-v2_i
    
    # # TODO: normalize p such that the L2-norm is 1
    # real_part = real_part/norm(real_part, ord=2)
    # real_part = real_part/norm(real_part, np.inf)
    # # TODO: normalize p such that the L^infty-Norm is 1
    # # imag_part = imag_part/norm(imag_part, np.inf)
    
    imag_part = v1_i+v2_r
    
    # # TODO: normalize p such that the L2-norm is 1
    # imag_part = imag_part/norm(imag_part, ord=2)
    # imag_part = imag_part/norm(imag_part, np.inf)
    # # TODO: normalize p such that the L^infty-Norm is 1
    # # imag_part = imag_part/norm(imag_part, np.inf)
    
    # shallow copy of vr
    x = vr.copy() 
    x[0::2] = real_part
    x[1::2] = imag_part
    
    # # TODO: normalize p such that the L2-norm is 1
    # x = x/norm(x, ord=2)
    # x = x/norm(x, np.inf)
    
    # TODO: normalize p such that the L^infty-Norm is 1
    # x = x/norm(x, np.inf)
    
    CG = dolf.FiniteElement('CG', mesh.ufl_cell(), degree)
    W = dolf.FunctionSpace(mesh, dolf.MixedElement([CG, CG]))
    p = dolf.Function(W)
    p.vector().set_local(x.getArray())
    # p.vector().apply('insert')
    (p_r, p_i) = p.split(True)
    
    dx = dolf.Measure('dx', domain=mesh)
    print(dolf.assemble((p_r*p_r+p_i*p_i)*dx))
    p_r = abs(p_r)/np.sqrt(dolf.assemble((p_r*p_r+p_i*p_i)*dx))
    p_i = p_i/np.sqrt(dolf.assemble((p_r*p_r+p_i*p_i)*dx))
    
    return p_r, p_i

