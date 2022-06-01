import dolfin as dolf
import numpy as np
from petsc4py import PETSc
import params
from helmholtz_pkg.petsc4py_utils import mult_complex_scalar_real_matrix

class derivativeA():
    def __init__(self, passiveflame):

        self.mesh = passiveflame.mesh
        self.degree = passiveflame.degree
        self.function_space = passiveflame.function_space
        self.u = passiveflame.u
        self.v = passiveflame.v

        self._A_c = None
        self._A_v = None
    
    @property
    def c_form(self):
        return self._A_c
    @property
    def v_form(self):
        return self._A_v
    
    def derivative_v(self):

        dx = dolf.Measure('dx', domain = self.mesh)
        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        da_11 = - params.gamma * dolf.dot(dolf.grad(v_1), dolf.grad(u_1)) * dx
        da_22 = - params.gamma * dolf.dot(dolf.grad(v_2), dolf.grad(u_2)) * dx
        
        self._A_v = da_11 + da_22

    def derivative_c(self):

        dx = dolf.Measure('dx', domain = self.mesh)
        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        da_11 = -2 * params.c_perturbed * dolf.dot(dolf.grad(v_1), dolf.grad(u_1)) * dx
        da_22 = -2 * params.c_perturbed * dolf.dot(dolf.grad(v_2), dolf.grad(u_2)) * dx
        
        self._A_c = da_11 + da_22


class derivativeD():
    def __init__(self, AF, submatrices):

        self.comm = AF.comm
        self.mesh = AF.mesh
        self.subdomains = AF.subdomains
        self.x_f = AF.x_f
        self.x_r = AF.x_r
        self.rho_u = AF.rho_u
        self.Q = AF.Q 
        self.U = AF.U
        self.FTF = AF.FTF
        self.degree = AF.degree
        self.function_space = AF.function_space
        # self.coeff = (AF.gamma - 1) / AF.rho_u * AF.Q / AF.U
        
        self.D_kj = submatrices
        self._D_v = None
        self._D_c = None

    @property
    def v(self):
        return self._D_v

    @property
    def c(self):
        return self._D_c
    
    # def derivative_v(self, omega):

    #     (D_11, D_12, D_21, D_22) = self.D_kj
    #     z = self.FTF(omega)
    #     self._D_v = mult_complex_scalar_real_matrix(z, D_11, D_12, D_21, D_22)
    #     coeff = self.Q/self.U / (params.p_amb/params.p_amb) * (params.gamma - 1)
    #     self._D_v = coeff * self._D_v

    def derivative_c(self, omega):

        (D_11, D_12, D_21, D_22) = self.D_kj
        z = self.FTF(omega)
        self._D_c = mult_complex_scalar_real_matrix(z, D_11, D_12, D_21, D_22)
        coeff = self.Q/self.U *(2*params.c_in)/(params.p_amb/params.p_amb) * (params.gamma - 1)/params.gamma
        self._D_c = coeff * self._D_c

    # def _localassemble_left(self):

    #     (v_1, v_2) = dolf.TestFunction(self.function_space)
    #     dx = dolf.Measure('dx', domain=self.mesh)

    #     V = dolf.FunctionSpace(self.mesh, 'CG', 1)
    #     const = dolf.interpolate(dolf.Constant(1), V)
    #     V_fl = dolf.assemble(const * dx)

    #     a_1 =  dolf.local_assemble(v_1 / V_fl * dx)
    #     a_2 =  dolf.local_assemble(v_2 / V_fl * dx)

    #     a = a_1 + a_2

    #     return a

