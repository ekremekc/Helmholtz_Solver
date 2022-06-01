from slepc4py import SLEPc
import numpy as np

def results(E):

    print()
    print("******************************")
    print("*** SLEPc Solution Results ***")
    print("******************************")
    print()

    its = E.getIterationNumber()
    print("Number of iterations of the method: %d" % its)

    eps_type = E.getType()
    print("Solution method: %s" % eps_type)

    nev, ncv, mpd = E.getDimensions()
    print("Number of requested eigenvalues: %d" % nev)

    tol, maxit = E.getTolerances()
    print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))

    nconv = E.getConverged()
    print("Number of converged eigenpairs %d" % nconv)

    A = E.getOperators()[0]
    vr, vi = A.createVecs()

    if nconv > 0:
        print()
    for i in range(nconv):
        k = E.getEigenpair(i, vr, vi)
        print("%15f, %15f" % (k.real, k.imag))
    print()


def eps_solver(A, C, target, nev, print_results=False):

    E = SLEPc.EPS().create()

    C = - C
    E.setOperators(A, C)

    # Generalized non-Hermitian
    # E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)

    # spectral transformation
    st = E.getST()
    st.setType('sinvert')

    E.setTarget(target)
    E.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)  # TARGET_REAL or TARGET_IMAGINARY

    E.setDimensions(nev, SLEPc.DECIDE)
    E.setTolerances(1e-15)
    E.setFromOptions()

    E.solve()

    if print_results:
        results(E)

    return E


def pep_solver(A, B, C, target, nev, print_results=False):

    Q = SLEPc.PEP().create()

    operators = [A, B, C]
    Q.setOperators(operators)

    # Q.setProblemType(SLEPc.PEP.ProblemType.GENERAL)

    # spectral transformation
    st = Q.getST()
    st.setType('sinvert')

    Q.setTarget(target)
    Q.setWhichEigenpairs(SLEPc.PEP.Which.TARGET_MAGNITUDE)  # TARGET_REAL or TARGET_IMAGINARY

    Q.setDimensions(nev, SLEPc.DECIDE)
    Q.setTolerances(1e-15)
    Q.setFromOptions()

    Q.solve()

    if print_results:
        results(Q)

    return Q

def fixed_point_iteration_pep(operators, D, target, nev=2, i=0):
    
    # Starting parameters
    omega_0 = 0
    eta = 1
    k = 1
    max_iter = 50
    
    while eta>10e-12 and k<=max_iter:
        
        # Create Matrix D 
        D.assemble_matrix(omega_0)
        D_mat = D.matrix
        
        
        # Create Matrix (A-D) in order to use pep solver
        nlinA = operators.A-D_mat
        
        # Passing Matrices to the pep solver
        E = pep_solver(nlinA, operators.B, operators.C, target, nev, print_results=False)
        
        # Dum = E.getOperators()[0]
        vr, vi = nlinA.createVecs()

        eig = E.getEigenpair(i, vr, vi)
        omega = eig
        # Extract Eigenvalue
        # omega, c, rx, cx = E.get_eigenpair(i)
        
        # Evaluate the error 
        # eta = abs(omega-omega_0)
        eta = abs((abs(omega)-abs(omega_0)))/abs(omega)
        
        # Update eigenvalue
        omega_0 = omega
        print("Omega: ", omega,"\n",
              "Error: ", eta)
        
        # Update iteration number
        k+=1
    print("Final omega:", omega)
    
    return E 

def fixed_point_iteration_pep1(operators, D, target, nev=2, i=0,
                              tol=1e-8, maxiter=50,
                              print_results=False,
                              problem_type='direct'):

    A = operators.A
    C = operators.C
    B = operators.B
    if problem_type == 'adjoint':
        B = operators.B_adj

    omega = np.zeros(maxiter, dtype=complex)
    f = np.zeros(maxiter, dtype=complex)
    alpha = np.zeros(maxiter, dtype=complex)

    E = pep_solver(A, B, C, target, nev, print_results=print_results)
    vr, vi = A.getVecs()
    eig = E.getEigenpair(i, vr, vi)

    omega[0] = eig
    alpha[0] = .5

    domega = 2 * tol
    k = - 1

    # formatting
    s = "{:.0e}".format(tol)
    s = int(s[-2:])
    s = "{{:+.{}f}}".format(s)

    while abs(domega) > tol:

        k += 1

        D.assemble_matrix(omega[k], problem_type)
        D_Mat = D.matrix
        if problem_type == 'adjoint':
            D_Mat = D.adjoint_matrix

        nlinA = A - D_Mat

        E = pep_solver(nlinA, B, C, target, nev, print_results=print_results)
        eig = E.getEigenpair(i, vr, vi)

        f[k] = eig

        if k != 0:
            alpha[k] = 1 / (1 - ((f[k] - f[k-1]) / (omega[k] - omega[k-1])))

        omega[k+1] = alpha[k] * f[k] + (1 - alpha[k]) * omega[k]

        domega = omega[k+1] - omega[k]

        print('iter = {:2d},  omega = {}  {}j,  |domega| = {:.2e}'.format(
            k + 1, s.format(omega[k + 1].real), s.format(omega[k + 1].imag), abs(domega)
        ))

    return E
# if __name__ == '__main__':
