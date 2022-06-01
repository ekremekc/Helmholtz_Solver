from slepc4py import SLEPc


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


# if __name__ == '__main__':
