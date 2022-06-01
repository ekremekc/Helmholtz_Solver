import dolfin as dolf
# from petsc4py import PETSc


def normalize_1(mesh, vr, degree=1):
    """
    real (, with real eigenvalues and eigenvectors)
    :param mesh:
    :param vr: <class 'petsc4py.PETSc.Vec'>
    :param degree:
    :return: <class 'dolfin.function.function.Function'>
    """

    V = dolf.FunctionSpace(mesh, 'CG', degree)
    p = dolf.Function(V)

    # vr.getArray() converts a petsc4py.PETSc.Vec into a numpy.ndarray
    # p.vector() converts a dolfin.function.function.Function into a dolfin.cpp.la.PETScVector
    # set (local) values of a dolfin Function
    p.vector().set_local(vr.getArray())
    # p.vector().apply('insert')

    # TODO: normalize p such that the L^infty-Norm is 1
    # TODO: normalize p such that the L2-norm is 1

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

    # shallow copy of vr
    x = vr.copy() 
    x[0::2] = vr_1
    x[1::2] = vi_1

    CG = dolf.FiniteElement('CG', mesh.ufl_cell(), degree)
    W = dolf.FunctionSpace(mesh, dolf.MixedElement([CG, CG]))
    p = dolf.Function(W)
    p.vector().set_local(x.getArray())
    # p.vector().apply('insert')
    (p_r, p_i) = p.split(True)
    # print(x)

    # TODO: normalize p such that the L^infty-Norm is 1
    # TODO: normalize p such that the L2-norm is 1

    return p_r, p_i
    # return p

# TODO: write normalize_3

# def normalize_3(mesh, vr, vi, degree=1):
#     """
#     complex, with complex eigenvalues and eigenvectors
#     """
#
#     return p_r, p_i

