import dolfin as dolf
import numpy as np
from petsc4py import PETSc

from helmholtz_pkg.petsc4py_utils import mult_complex_scalar_real_matrix


class ActiveFlame:

    gamma = 1.4

    def __init__(self, mesh, subdomains, x_r, rho_u, Q, U, FTF, degree=1, comm=None,
                 constrained_domain=None):
        """
        The Active Flame class which consist required matrices for solution
        of nonlinear eigenvalue problem.
        
        Example: D = ActiveFlame(mesh, subdomains,
                    x_f, x_r, rho_in, 1., 1., ftf,
                    degree=1)
        In order to obtain matrices, D.assemble_submatrices() should be called.

        Parameters
        ----------
        mesh : dolfin.cpp.mesh.Mesh
            Mesh of the entire domain.
        subdomains : dolfin.cpp.mesh.MeshFunctionSizet
            Relevant subdomain for flame region.
        x_f : float
            Not used in this class?
        x_r : np.array
            flame location vector
        rho_u : float
            flame upstream density
        Q : float
            Heat release, default is one because it has calculated in parameter n
        U : float
            Bulk velocity, default is one because it has calculated in parameter n
        FTF : function
            Flame transfer function (n*e^(iw\tau)).
        degree : int, optional
            Degree of basis functions. The default is 1.
        comm : TYPE, optional
            DESCRIPTION. The default is None.
        constrained_domain : TYPE, optional
            DESCRIPTION. The default is None.

        """

        self.comm = comm

        self.mesh = mesh
        self.subdomains = subdomains
        self.x_r = x_r
        self.rho_u = rho_u
        self.Q = Q
        self.U = U
        self.FTF = FTF
        self.degree = degree
        self.r = dolf.Expression("x[1]", degree = degree)
        self.coeff = (self.gamma - 1) / rho_u * Q / U

        # __________________________________________________

        self._a = {}
        self._b = {}
        self._D_kj = None
        self._D_kj_adj = None
        self._D = None
        self._D_adj = None

        # __________________________________________________

        CG = dolf.FiniteElement('CG', mesh.ufl_cell(), degree)
        W = dolf.FunctionSpace(mesh, dolf.MixedElement([CG, CG]), constrained_domain=constrained_domain)

        self.function_space = W

        for fl, x in enumerate(self.x_r):
            self._a[str(fl)] = self._assemble_left_vector(fl)
            self._b[str(fl)] = self._assemble_right_vector(x)

    @property
    def submatrices(self):
        return self._D_kj

    @property
    def matrix(self):
        return self._D
    @property
    def a(self):
        return self._a
    @property
    def b(self):
        return self._b        
    @property
    def adjoint_submatrices(self):
        return self._D_kj_adj

    @property
    def adjoint_matrix(self):
        return self._D_adj

    @staticmethod
    def _helper_func(v, dofmap):
        """
        Returns the nonzero indices and values of matrix a
        in order to perform efficient matrix multiplication

        Parameters
        ----------
        v : <class 'tuple'>
            includes mixed elements of a_1 and a_2
        dofmap : dolfin.cpp.fem.DofMap
            mapping that maps the dofs

        Returns
        -------
        v : np.array
            2D array which includes index and dofs of nonzero cells

        """
        indices = np.flatnonzero(v) # Returns the indices which are nonzero
        values = v[indices]
        my_list = []
        for index, value in zip(indices, values):
            my_list.append([dofmap.dofs()[index], value])
        v = np.array(my_list)
        return v

    def _assemble_left_vector(self, fl):
        """
        Assembles \int v(x) \phi_k dV

        Parameters
        ----------
        fl : int
            flame tag

        Returns
        -------
        v : <class 'tuple'>
            includes assembled elements of a_1 and a_2

        """

        (v_1, v_2) = dolf.TestFunction(self.function_space)
        dV = dolf.Measure('dx', subdomain_data=self.subdomains)

        V = dolf.FunctionSpace(self.mesh, 'CG', 1)
        area = dolf.interpolate(dolf.Constant(1), V)
        V_fl =   dolf.assemble(2*np.pi*area * self.r * dV(fl))
        # V_fl =   dolf.assemble(area * self.r * dV(fl))
        print(fl,"Flame volume: ",V_fl,
              "Theroetical volume: ", 0.05*np.pi*(0.047/2)**2)
        a_1 =  dolf.assemble(v_1 / V_fl * self.r * dV(fl))
        a_2 =  dolf.assemble(v_2 / V_fl * self.r * dV(fl))
        
        
        dofmap = self.function_space.dofmap()

        a_1 = self._helper_func(a_1, dofmap)
        a_2 = self._helper_func(a_2, dofmap)
        
        a = (a_1, a_2)

        return a

    def _assemble_right_vector(self, x):
        """
        Calculates degree of freedoms and indices of 
        right vector of
        
        \nabla(\phi_j(x_r)) . n
        
        which includes gradient value of test fuunction at
        reference point x_r
        
        Parameters
        ----------
        x : np.array
            flame location vector

        Returns
        -------
        np.array
            Array of degree of freedoms and indices as vector b.

        """

        v = np.array([[0, 0, 1]])  # row
        dimension = self.mesh.geometric_dimension()
        if dimension == 1:
            v = np.array([[1]])
        elif dimension == 2:
            v = np.array([[1, 0]])
        # else:  # elif dimension == 3:
        #     pass
        v = v.transpose()  # column

        b = [np.array([]), np.array([])]

        cell_index = self.mesh.bounding_box_tree().compute_first_entity_collision(dolf.Point(*x))
        # print("cell index is: ",cell_index, "max cell number: ",self.mesh.num_entities(self.mesh.topology().dim()))
        if cell_index <= self.mesh.num_entities(self.mesh.topology().dim()):
            # print("Cell is in Mesh")
            cell = dolf.Cell(self.mesh, cell_index)

            b = []

            for j in range(2):

                dofmap = self.function_space.sub(j).dofmap()
                cell_dofs = dofmap.cell_dofs(cell_index)

                element = self.function_space.sub(j).element()
                d_dx = element.evaluate_basis_derivatives_all(1, x, cell.get_vertex_coordinates(), cell.orientation())

                d_dx = d_dx.reshape((len(cell_dofs), -1))
                # print(d_dx)
                d_dv = np.dot(d_dx, v)
                # print(d_dv)
                d_dv = d_dv[:, 0]
                # print(d_dv)
                my_list = []

                for i, dof in enumerate(cell_dofs):
                    my_list.append([dofmap.tabulate_local_to_global_dofs()[dof], d_dv[i]])

                my_vec = np.array(my_list)
                #print(my_vec)
                b.append(my_vec)
            #
        else:
            
            raise ValueError("The cell for the subdomain is not in the mesh, please check the mesh and subdomain.")
            
        return (*b, )

    @staticmethod
    def _csr_matrix(a, b):

        # len(a) and len(b) are not the same

        nnz = len(a) * len(b)
        

        row = np.zeros(nnz)
        col = np.zeros(nnz)
        val = np.zeros(nnz)

        for i, c in enumerate(a):
            for j, d in enumerate(b):
                row[i * len(b) + j] = c[0]
                col[i * len(b) + j] = d[0]
                val[i * len(b) + j] = c[1] * d[1]

        row = row.astype(dtype='int32')
        col = col.astype(dtype='int32')

        return row, col, val

    def assemble_submatrices(self, problem_type='direct'):
        """
        This function handles efficient cross product of the 
        vectors a and b calculated above and generates highly sparse 
        matrix D_kj which represents active flame matrix without FTF and
        other constant multiplications.

        Parameters
        ----------
        problem_type : str, optional
            Specified problem type. The default is 'direct'.
            Matrix can be obtained by selecting problem type, other
            option is adjoint.
        
        """

        num_fl = len(self.x_r)  # number of flames
        global_size = self.function_space.dim()
        local_size = len(self.function_space.dofmap().dofs())

        D_kj = dict()

        for k in range(2):
            for j in range(2):

                row = dict()
                col = dict()
                val = dict()

                for fl in range(num_fl):

                    u = None
                    v = None

                    if problem_type == 'direct':
                        u = self._a[str(fl)][k]  # column vector
                        v = self._b[str(fl)][j]  # row vector

                    elif problem_type == 'adjoint':
                        u = self._b[str(fl)][k]
                        v = self._a[str(fl)][j]

                    row[str(fl)], col[str(fl)], val[str(fl)] = self._csr_matrix(u, v)

                row = np.concatenate([row[str(fl)] for fl in range(num_fl)])
                col = np.concatenate([col[str(fl)] for fl in range(num_fl)])
                val = np.concatenate([val[str(fl)] for fl in range(num_fl)])

                i = np.argsort(row)

                row = row[i]
                col = col[i]
                val = val[i]

                indptr = np.bincount(row, minlength=local_size)
                indptr = np.insert(indptr, 0, 0).cumsum()
                indptr = indptr.astype(dtype='int32')

                mat = PETSc.Mat().create(comm=self.comm)
                mat.setSizes([(local_size, global_size), (local_size, global_size)])
                mat.setType('aij')
                mat.setUp()

                mat.setValuesCSR(indptr, col, val)
                mat.assemblyBegin()
                mat.assemblyEnd()

                D_kj['{0}{1}'.format(k + 1, j + 1)] = mat

        if problem_type == 'direct':
            self._D_kj = (D_kj['11'], D_kj['12'], D_kj['21'], D_kj['22'])
        elif problem_type == 'adjoint':
            self._D_kj_adj = (D_kj['11'], D_kj['12'], D_kj['21'], D_kj['22'])

    def assemble_matrix(self, omega, problem_type='direct'):
        """
        This function handles the multiplication of the obtained matrix D
        with Flame Transfer Function and other constants.
        The operation is
        
        D = (gamma - 1) / rho_u * Q / U * D_kj       
        
        At the end the matrix D is petsc4py.PETSc.Mat
        Parameters
        ----------
        omega : complex
            eigenvalue that found by solver
        problem_type : str, optional
            Specified problem type. The default is 'direct'.
            Matrix can be obtained by selecting problem type, other
            option is adjoint.

        Returns
        -------
        petsc4py.PETSc.Mat

        """

        if problem_type == 'direct':

            (D_11, D_12, D_21, D_22) = self._D_kj

            z = self.FTF(omega)
            self._D = mult_complex_scalar_real_matrix(z, D_11, D_12, D_21, D_22)
            self._D = self.coeff * self._D

        elif problem_type == 'adjoint':

            (D_11, D_12, D_21, D_22) = self._D_kj_adj

            z = np.conj(self.FTF(np.conj(omega)))
            self._D_adj = mult_complex_scalar_real_matrix(z, D_11, D_12, D_21, D_22)
            self._D_adj = self.coeff * self._D_adj

    def get_derivative(self, omega):

        """ Derivative of the unsteady heat release (flame) operator D
        wrt the eigenvalue (complex angular frequency) omega."""

        (D_11, D_12, D_21, D_22) = self._D_kj
        z = self.FTF(omega, k=1)
        dD_domega = mult_complex_scalar_real_matrix(z, D_11, D_12, D_21, D_22)
        dD_domega = self.coeff * dD_domega

        return dD_domega


# if __name__ == '__main__':
