from dolfin import *
import numpy as np
from petsc4py import PETSc


def helper_1(v, d, dof_coordinates):
    # Convert 'dolfin.cpp.la.Vector' into 'numpy.ndarray'
    v = v.get_local()
    a = []
    for i, item in enumerate(dof_coordinates):
        if d(item) <= 1. + DOLFIN_EPS:
            a.append([i, v[i]])
    a = np.array(a)
    return a


def helper_2(m, a, b):

    # m, len(a) and len(b) are not the same

    nnz = len(a)*len(b)

    row = np.zeros(nnz)
    col = np.zeros(nnz)
    val = np.zeros(nnz)

    for i, c in enumerate(a):
        for j, d in enumerate(b):
            col[i * len(b) + j] = d[0]
            row[i * len(b) + j] = c[0]
            val[i * len(b) + j] = c[1] * d[1]

    row = row.astype(dtype='int32')
    col = col.astype(dtype='int32')

    indptr = np.bincount(row, minlength=m)
    indptr = np.insert(indptr, 0, 0).cumsum()
    indptr = indptr.astype(dtype='int32')

    A = PETSc.Mat().createAIJ(size=(m, m), csr=(indptr, col, val))

    return A


def mult_complex_scalar_real_matrix(z, A_11, A_12, A_21, A_22):
    """
    multiply a real matrix by a complex scalar
    """
    a = z.real
    b = z.imag

    B = a*A_11 - b*A_12 + b*A_21 + a*A_22

    return B


class ConstTimeDelay(object):
    def __init__(self, mesh, degree=1, gamma=1.4, n=None, tau=None):

        self.gamma = gamma
        self.n     = n
        self.tau   = tau

        CG = FiniteElement("CG", mesh.ufl_cell(), degree)
        W = FunctionSpace(mesh, CG * CG)
        self.r = Expression('x[1]', degree=1)
        self.dof_coordinates = W.tabulate_dof_coordinates()

        (self.p_1, self.p_2) = TrialFunction(W)
        (self.q_1, self.q_2) = TestFunction(W)

        self.a        = None
        self.b        = None
        self.D_jk_dir = None
        self.D_jk_adj = None
        self.D_jk     = None  # VarTimeDelay
        self.D_dir    = None
        self.D_adj    = None
        self.D        = None  # VarTimeDelay

    def assemble_left_vector(self, f, d):

        a_1 = self.q_1 * f * self.r * dx
        a_2 = self.q_2 * f * self.r * dx

        a_1 = assemble(a_1)
        a_2 = assemble(a_2)

        a_1 = helper_1(a_1, d, self.dof_coordinates)
        a_2 = helper_1(a_2, d, self.dof_coordinates)

        self.a = (a_1, a_2)

    def assemble_right_vector(self, g, d):

        b_1 = g * self.p_1.dx(0) * self.r * dx
        b_2 = g * self.p_2.dx(0) * self.r * dx

        b_1 = 2 * pi * assemble(b_1)
        b_2 = 2 * pi * assemble(b_2)

        b_1 = helper_1(b_1, d, self.dof_coordinates)
        b_2 = helper_1(b_2, d, self.dof_coordinates)

        self.b = (b_1, b_2)

    def assemble_submatrices(self, input_string='direct'):

        (a_1, a_2) = self.a
        (b_1, b_2) = self.b

        m = len(self.dof_coordinates)

        if input_string == 'direct':

            D_11 = helper_2(m, a_1, b_1)
            D_12 = helper_2(m, a_1, b_2)
            D_21 = helper_2(m, a_2, b_1)
            D_22 = helper_2(m, a_2, b_2)

            self.D_jk_dir = (D_11, D_12, D_21, D_22)

        elif input_string == 'adjoint':

            D_11 = helper_2(m, b_1, a_1)
            D_12 = helper_2(m, b_1, a_2)
            D_21 = helper_2(m, b_2, a_1)
            D_22 = helper_2(m, b_2, a_2)
    
            self.D_jk_adj = (D_11, D_12, D_21, D_22)

    def assemble_matrix(self, omega, input_string='direct'):

        omega_r = omega.real
        omega_i = omega.imag

        if input_string == 'direct':

            (D_11, D_12, D_21, D_22) = self.D_jk_dir

            r =  exp(- omega_i*self.tau) * cos(omega_r*self.tau)
            c =  exp(- omega_i*self.tau) * sin(omega_r*self.tau)
            z = r + 1j*c
            self.D_dir = mult_complex_scalar_real_matrix(z, D_11, D_12, D_21, D_22)
            self.D_dir = (self.gamma - 1)/self.gamma * self.n * self.D_dir

        elif input_string == 'adjoint':

            (D_11, D_12, D_21, D_22) = self.D_jk_adj

            r =   exp(  omega_i*self.tau) * cos(omega_r*self.tau)
            c = - exp(  omega_i*self.tau) * sin(omega_r*self.tau)
            z = r + 1j*c
            self.D_adj = mult_complex_scalar_real_matrix(z, D_11, D_12, D_21, D_22)
            self.D_adj = (self.gamma - 1)/self.gamma * self.n * self.D_adj

    def get_left_vector(self):
        return self.a

    def get_right_vector(self):
        return self.b

    def get_submatrices(self, input_string='direct'):
        if input_string == 'direct':
            return self.D_jk_dir
        elif input_string == 'adjoint':
            return self.D_jk_adj

    def get_matrix(self, input_string='direct'):
        if input_string == 'direct':
            return self.D_dir
        elif input_string == 'adjoint':
            return self.D_adj

    def tmp_name(self, omega):

        (D_11, D_12, D_21, D_22) = self.D_jk_dir

        omega_r = omega.real
        omega_i = omega.imag

        r = - self.tau * exp(- omega_i*self.tau) * sin(omega_r*self.tau)
        c =   self.tau * exp(- omega_i*self.tau) * cos(omega_r*self.tau)
        z = r + 1j*c
        dD_domega = mult_complex_scalar_real_matrix(z, D_11, D_12, D_21, D_22)
        dD_domega = (self.gamma - 1)/self.gamma * self.n * dD_domega

        return dD_domega

    def assemble_left_vector_local(self, f, cell):

        a_1 = self.q_1 * f * dx
        a_2 = self.q_2 * f * dx
        a   = a_1 + a_2

        a = assemble_local(a, cell)

        return a

    def assemble_right_vector_local(self, g, cell):

        b_1 = g * self.p_1.dx(0) * dx
        b_2 = g * self.p_2.dx(0) * dx
        b   = b_1 + b_2
        
        b = assemble_local(b, cell)

        return b


class VarTimeDelay(ConstTimeDelay):

    def assemble_left_vector(self, f, d):

        f_r, f_i = f  # check f is a tuple (or a list)

        a_1r = self.q_1 * f_r * dx
        a_1i = self.q_1 * f_i * dx
        a_2r = self.q_2 * f_r * dx
        a_2i = self.q_2 * f_i * dx

        a_1r = assemble(a_1r)
        a_1i = assemble(a_1i)
        a_2r = assemble(a_2r)
        a_2i = assemble(a_2i)

        a_1r = helper_1(a_1r, d, self.dof_coordinates) 
        a_1i = helper_1(a_1i, d, self.dof_coordinates)
        a_2r = helper_1(a_2r, d, self.dof_coordinates)
        a_2i = helper_1(a_2i, d, self.dof_coordinates)

        self.a = (a_1r, a_1i, a_2r, a_2i)

    def assemble_submatrices(self, input_string='direct'):

        (a_1r, a_1i, a_2r, a_2i) = self.a
        (b_1, b_2) = self.b

        m = len(self.dof_coordinates)

        if input_string == 'direct':

            D_11 = helper_2(m, a_1r, b_1)
            D_12 = helper_2(m, a_1i, b_2)
            D_21 = helper_2(m, a_2i, b_1)
            D_22 = helper_2(m, a_2r, b_2)

            self.D_jk = (D_11, D_12, D_21, D_22)

        elif input_string == 'adjoint':

            D_11 = helper_2(m, b_1, a_1r)
            D_12 = helper_2(m, b_1, a_2i)
            D_21 = helper_2(m, b_2, a_1i)
            D_22 = helper_2(m, b_2, a_2r)
    
            self.D_jk = (D_11, D_12, D_21, D_22)

    def assemble_matrix(self):

        (D_11, D_12, D_21, D_22) = self.D_jk

        self.D = D_11 - D_12 + D_21 + D_22
        self.D = (self.gamma - 1)/self.gamma * self.n * self.D

    def get_submatrices(self):
        return self.D_jk

    def get_matrix(self):
        return self.D
