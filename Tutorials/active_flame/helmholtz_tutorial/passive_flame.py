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
            Y_r, Y_i = Y.real, Y.imag

            b_11 = - Y_i * c * v_1 * u_1 * ds(i)
            b_12 = - Y_r * c * v_1 * u_2 * ds(i)
            b_21 = + Y_r * c * v_2 * u_1 * ds(i)
            b_22 = - Y_i * c * v_2 * u_2 * ds(i)

            b_ = b_11 + b_12 + b_21 + b_22

            integrals_R.append(b_)

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


class PassiveFlame:

    def __init__(self, mesh, boundaries, boundary_conditions,
                 c, degree=1, constrained_domain=None):

        self.mesh = mesh
        # self.boundaries = boundaries
        self.boundary_conditions = boundary_conditions
        self.c = c
        # self.degree = degree

        self.dx = dolf.Measure('dx', domain=mesh)
        self.ds = dolf.Measure('ds', domain=mesh, subdomain_data=boundaries)

        CG = dolf.FiniteElement('CG', mesh.ufl_cell(), degree)
        W = dolf.FunctionSpace(mesh, dolf.MixedElement([CG, CG]), constrained_domain=constrained_domain)
        self.function_space = W  #

        self.u = dolf.TrialFunction(W)
        self.v = dolf.TestFunction(W)

        self.bcs = []

        for i in boundary_conditions:
            if 'Dirichlet' in boundary_conditions[i]:
                for j in range(2):
                    bc = dolf.DirichletBC(W.sub(j), 0.0, boundaries, i)
                    self.bcs.append(bc)

        self._A = None
        self._B = None
        self._C = None

    @property
    def A(self):
        return self._A

    @property
    def B(self):
        return self._B

    @property
    def C(self):
        return self._C

    def assemble_A(self):

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        a_11 = - self.c ** 2 * dot(grad(v_1), grad(u_1)) * self.dx
        a_22 = - self.c ** 2 * dot(grad(v_2), grad(u_2)) * self.dx
        a_ = a_11 + a_22

        dummy = (v_1 + v_2) * self.dx

        A, b = dolf.assemble_system(a_, dummy, self.bcs)
        A = dolf.as_backend_type(A).mat()
        self._A = A

    def assemble_B(self):

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        integrals_R = []

        for i in self.boundary_conditions:
            if 'Robin' in self.boundary_conditions[i]:

                Y = self.boundary_conditions[i]['Robin']
                Y_r, Y_i = Y.real, Y.imag

                b_11 = - Y_i * self.c * v_1 * u_1 * self.ds(i)
                b_12 = - Y_r * self.c * v_1 * u_2 * self.ds(i)
                b_21 = + Y_r * self.c * v_2 * u_1 * self.ds(i)
                b_22 = - Y_i * self.c * v_2 * u_2 * self.ds(i)

                b_ = b_11 + b_12 + b_21 + b_22

                integrals_R.append(b_)

        if integrals_R:

            b_ = sum(integrals_R)
            B = dolf.assemble(b_)
            B = dolf.as_backend_type(B).mat()

        else:

            N = self.function_space.dim()  # global size
            istart, iend = self.function_space.dofmap().ownership_range()
            n = iend - istart  # local size

            # empty matrix
            B = PETSc.Mat().create()
            B.setSizes([(n, N), (n, N)])
            B.setFromOptions()
            B.setUp()
            B.assemble()

        self._B = B

    def assemble_C(self):

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        c_11 = v_1 * u_1 * self.dx
        c_22 = v_2 * u_2 * self.dx
        c_   = c_11 + c_22

        dummy = (v_1 + v_2) * self.dx

        C, b = dolf.assemble_system(c_, dummy, self.bcs)
        [bc.zero(C) for bc in self.bcs]
        C = dolf.as_backend_type(C).mat()
        self._C = C


# if __name__ == '__main__':
