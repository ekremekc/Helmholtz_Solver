import dolfin as dolf
# from petsc4py import PETSc
from slepc4py import SLEPc

import numpy as np


def normalize_1(mesh, vr, degree=1, constrained_domain=None):
    """
    real (, with real eigenvalues and eigenvectors)
    normalize p such that the L2-Norm is 1
    :param mesh:
    :param vr: <class 'petsc4py.PETSc.Vec'>
    :param degree:
    :param constrained_domain: periodic boundary conditions
    :return: <class 'dolfin.function.function.Function'>
    """

    V = dolf.FunctionSpace(mesh, 'CG', degree, constrained_domain=constrained_domain)

    # print(vr.getArray())

    # index, value = abs(vr).max()
    # # print(index, value)
    #
    # vr /= value
    # # print(vr.getArray())

    p = dolf.Function(V)
    p.vector().set_local(vr.getArray())
    p.vector().apply('insert')

    dx = dolf.Measure('dx')
    meas = dolf.assemble(p * p * dx)
    meas = np.sqrt(meas)
    # print(meas)
    vr /= meas
    p.vector().set_local(vr.getArray())
    p.vector().apply('insert')
    # meas = dolf.assemble(p * p * dx)
    # meas = np.sqrt(meas)
    # print(meas)

    return p


def normalize_2(mesh, vr, vi, degree=1, constrained_domain=None):
    """
    complex, with complex eigenvalues and eigenvectors
    normalize p such that the L2-norm is 1
    :param mesh:
    :param vr: <class 'petsc4py.PETSc.Vec'>
    :param vi: <class 'petsc4py.PETSc.Vec'>
    :param degree:
    :param constrained_domain: periodic boundary conditions
    :return: <class 'dolfin.function.function.Function'>
    """

    # local_y = np.max(np.abs(y))
    # global_y = MPI.COMM_WORLD.allreduce(sendobj=local_y, op=MPI.MAX)
    # y /= global_y

    V = dolf.FunctionSpace(mesh, "CG", degree, constrained_domain=constrained_domain)
    CG = dolf.FiniteElement("CG", mesh.ufl_cell(), degree)
    W = dolf.FunctionSpace(mesh, dolf.MixedElement([CG, CG]), constrained_domain=constrained_domain)

    vr_1 = vr.getArray()[0::2]
    vr_2 = vr.getArray()[1::2]
    vi_1 = vi.getArray()[0::2]
    vi_2 = vi.getArray()[1::2]

    x = vr_1 - vi_2 + 1j * (vr_2 + vi_1)

    x_r = x.real
    x_i = x.imag

    p_r = dolf.Function(V)
    p_i = dolf.Function(V)

    p_r.vector().set_local(x_r)
    p_r.vector().apply('insert')
    p_i.vector().set_local(x_i)
    p_i.vector().apply('insert')

    dx = dolf.Measure('dx')
    meas = dolf.assemble((p_r * p_r + p_i * p_i) * dx)
    meas = np.sqrt(meas)
    # print(meas)

    x /= meas

    x_r = x.real
    x_i = x.imag

    # p_r.vector().set_local(x_r)
    # p_r.vector().apply('insert')
    # p_i.vector().set_local(x_i)
    # p_i.vector().apply('insert')
    #
    # meas = dolf.assemble((p_r * p_r + p_i * p_i) * dx)
    # meas = np.sqrt(meas)
    # # print(meas)

    # shallow copy of vr
    x = vr.copy()

    istart, iend = x.getOwnershipRange()

    x[istart:iend:2] = x_r
    x[istart+1:iend+1:2] = x_i

    p = dolf.Function(W)
    p.vector().set_local(x.getArray())
    p.vector().apply('insert')

    return p


def normalize_eigenvector(mesh, obj, i, degree=1):

    A = obj.getOperators()[0]
    vr, vi = A.createVecs()

    if isinstance(obj, SLEPc.EPS):
        eig = obj.getEigenvalue(i)
        omega = np.sqrt(eig)
        obj.getEigenvector(i, vr, vi)

    elif isinstance(obj, SLEPc.PEP):
        eig = obj.getEigenpair(i, vr, vi)
        omega = eig

    if np.isclose(omega.imag, 0):
        omega = omega.real
        p = normalize_1(mesh, vr, degree)
    else:
        p = normalize_2(mesh, vr, vi, degree)

    return omega, p
