import dolfin as dolf
from petsc4py import PETSc

dot = dolf.dot
grad = dolf.grad


def passive_flame(mesh, boundaries, boundary_conditions, c, degree=1):

    V = dolf.FunctionSpace(mesh, 'CG', degree)

    dx = dolf.Measure('dx', domain=mesh)

    bcs = []

    for i in boundary_conditions:
        if 'Dirichlet' in boundary_conditions[i]:
            bc = dolf.DirichletBC(V, 0.0, boundaries, i)
            bcs.append(bc)

    u = dolf.TrialFunction(V)
    v = dolf.TestFunction(V)

    a_ = - c ** 2 * dot(grad(u), grad(v)) * dx
    c_ = u * v * dx

    dummy = v * dx

    A, b = dolf.assemble_system(a_, dummy, bcs)
    # convert dolfin.cpp.la.Matrix into petsc4py.PETSc.Mat
    A = dolf.as_backend_type(A).mat()
    

    C, b = dolf.assemble_system(c_, dummy, bcs)
    [bc.zero(C) for bc in bcs]
    C = dolf.as_backend_type(C).mat()

    N = V.dim()  # global size
    istart, iend = V.dofmap().ownership_range()
    n = iend - istart  # local size

    # empty matrix
    B = PETSc.Mat().create()
    B.setSizes([(n, N), (n, N)])
    B.setFromOptions()
    B.setUp()
    B.assemble()

    return A, B, C


def passive_flame_mixed(mesh, boundaries, boundary_conditions, c, degree=1):

    CG = dolf.FiniteElement('CG', mesh.ufl_cell(), degree)
    W = dolf.FunctionSpace(mesh, dolf.MixedElement([CG, CG]))

    dx = dolf.Measure('dx', domain=mesh)
    ds = dolf.Measure('ds', domain=mesh, subdomain_data=boundaries)

    bcs = []

    for i in boundary_conditions:
        if 'Dirichlet' in boundary_conditions[i]:
            for j in range(2):
                bc = dolf.DirichletBC(W.sub(j), 0.0, boundaries, i)
                bcs.append(bc)

    (u_1, u_2) = dolf.TrialFunction(W)
    (v_1, v_2) = dolf.TestFunction(W)

    a_11 = - c ** 2 * dot(grad(v_1), grad(u_1)) * dx
    a_22 = - c ** 2 * dot(grad(v_2), grad(u_2)) * dx
    a_ = a_11 + a_22

    c_11 = v_1 * u_1 * dx
    c_22 = v_2 * u_2 * dx
    c_ = c_11 + c_22

    dummy = (v_1 + v_2) * dx

    A, b = dolf.assemble_system(a_, dummy, bcs)
    A = dolf.as_backend_type(A).mat()

    C, b = dolf.assemble_system(c_, dummy, bcs)
    [bc.zero(C) for bc in bcs]
    C = dolf.as_backend_type(C).mat()

    integrals_R = []

    for i in boundary_conditions:
        if 'Robin' in boundary_conditions[i]:

            Y = boundary_conditions[i]['Robin']
            print(Y)
            # TODO: write ufl forms (equivalent of lines 68-70 or 72-74)
            Y_real, Y_imag = Y.real, Y.imag

            b_11 = -Y_imag*c*v_1*u_1*ds(i)
            b_12 = -Y_real*c*v_1*u_2*ds(i)
                
            b_21 =  Y_real*c*v_2*u_1*ds(i)
            b_22 = -Y_imag*c*v_2*u_2*ds(i)

            b_1 = b_11 + b_12
            b_2 = b_21 + b_22
            
            b_ = b_1 + b_2
      
            integrals_R.append(b_)
            
            print('line worked..')

    if integrals_R:

        b_ = sum(integrals_R)
        B = dolf.assemble(b_)
        B = dolf.as_backend_type(B).mat()

    else:

        N = W.dim()  # global size
        istart, iend = W.dofmap().ownership_range()
        n = iend - istart  # local size

        # empty matrix
        B = PETSc.Mat().create()
        B.setSizes([(n, N), (n, N)])
        B.setFromOptions()
        B.setUp()
        B.assemble()

    return A, B, C


# if __name__ == '__main__':
