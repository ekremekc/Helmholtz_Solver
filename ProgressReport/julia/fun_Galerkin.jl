include("./fun_normalize.jl")
include("./fun_nondim.jl")

using LinearAlgebra

#using .parameters

using .nondim_param:param # importing non-dimensionless variables

include("./structures.jl")
using .structures  # Importing emodes and ds structures




function fun_Galerkin(param = param, s0 = 1*im*pi, NG = 40)

    #  Galerkin model, solved as a nonlinear eigenvalue problem,
    # assuming p=0 at both ends
    #
    ##
    function fun_L(s,D2,I,t,F)
        L = s^2*I + exp(-s*t) * F - D2
        return L
    end
    ##
    function fun_dLds(s,I,t,F)
        dLds = 2*s*I - t * exp(-s*t) * F
        return dLds
    end
    ##
    function fun_dLdF(s,t)
        dLdF = exp(-s*t)
        return dLdF
    end
    ##
    function fun_dLdt(s,t,F)
        dLdt = -s*F*exp(-s*t)
        return dLdt
    end
    ##
    # INPUTS
    # param.x_m     measurement position (nondimensional)
    # param.x_f     heat releae position (nondimensional)
    # param.n       heat release index (nondimensional)
    # param.tau 	heat release time delay (nondimensional)
    # param.gam 	ratio of specific heats
    # param.an      parameter governing the width of heat and measurement zones
    #
    # scheme.s0     starting value of s
    # scheme.NG     number of Galerkin modes
    #
    # OUTPUTS:
    # emode.s       eigenvalue
    # emode.x       abscissa for eigenfunction
    # emode.P       pressure eigenfunction
    # emode.U       velocity eigenfunction
    #
    # ds.x_m        d(s)/d(param.x_m)
    # ds.x_f        d(s)/d(param.x_f)
    # ds.n       	d(s)/d(param.n)
    # ds.tau 		d(s)/d(param.tau)
    # ds.L          d(s)/d(L)

    # Unpack the number of Galerkin modes
    N = NG
    # Create a col vector pi*(1:N)
    kpi = pi*(1:N)
    # Create the second order differentiation matrix
    D2 = -Diagonal(kpi.^2)
    # Create the identity matrix
    Id = I(N)



    ## Calculate the internal parameters
    # Create the heat release envelope v(x)
    v = (2*kpi) .* sin.(kpi*param.x_f)
    # Create the measurement envelope w(x)
    w = cos.(kpi*param.x_m)
    # Wrap v, w, gamma, param.n into a matrix, F
    F = (v*w') * (param.gam-1)/param.gam * param.n
    # Extract tau
    t = param.tau

    ## Iterate to find s
    # Set tolerance, initial s, and dummy dels
    tol = 1e-8
    s = s0
    dels = 2*tol
    while abs(dels) > tol
        # Evaluate L
        L = fun_L(s,D2,Id,t,F)
        # Evaluate dL/ds
        dLds = fun_dLds(s,Id,t,F)
        # evaluate new s with Jacobi's formula
        dels = - 1/tr(L\dLds)
        # Update s
        s = s + dels
    end
    # Find the corresponding left and right eigenvectors for U
    L = fun_L(s,D2,Id,t,F)
    if N == 1
        U_dir = 1
        U_adj = 1
    else
        U_dir = nullspace(L)
        U_adj = nullspace(L')
    end
    # Evaluate P (the Galerkin coefficients of the pressure)
    P_dir = -param.gam*s*(U_dir./kpi)
    # Display the dimensional eigenvalue (frequency, then growth rate)
    println("fun_Galerkin: s = ", s, " (nondim)")

    ## Return zeros if null does not work
    if isempty(U_dir)
        U_dir = zeros(NG,1)
        U_adj = U_dir
    end

    ## Evaluate the eigenfunction on a user-defined grid
    # Set the number of elements on which to evaluate the eigenfunction
    Nx = 1000;
    # Set the abscissa
    x = range(0,1,length = Nx+1)
    # Create the matrices to evaluate u(x) from U_j and p(x) from P_j
    MU = cos.(x*kpi')
    MP = sin.(x*kpi')
    # Evaluate U(x)
    U = MU*U_dir;
    # Evaluate P(x)
    P = MP*P_dir
    P = vec(P)
    # Generate mass matrix on this abscissa
    dx = x[2]-x[1]
    M = dx*I(Nx+1)
    M[1,1] = M[1,1]/2
    M[Nx+1,Nx+1] = M[Nx+1,Nx+1]/2



    # Normalize P
    P,f = fun_normalize(P,M)
    # # Normalize U by the same factor
    U = U/f;

    emode = emodes(s,x,P,U)

    ## Wrap output into structure
    # Wrap eigenvalue
    emode.s = s;
    # Wrap eigenfunction
    emode.x = x; emode.P = P; emode.U = U;

    ## Return if sensitivities are not required
    #if nargout == 1; return; end

    ## Calculate the gradients of the internal parameters w.r.t. param.*
    # dv/d(param.x_f)
    dvdx_f = (2*kpi) .* kpi .*  cos.(kpi*param.x_f);
    # dw/d(param.x_m)
    dwdx_m =            kpi .* -sin.(kpi*param.x_m);
    # dF/d(param.n)
    dFdn   = F / param.n;
    # dF/d(param.x_f)
    dFdx_f = (dvdx_f*w') * (param.gam-1)/param.gam * param.n;
    # dF/d(param.x_m)
    dFdx_m = (v*dwdx_m') * (param.gam-1)/param.gam * param.n;

    ## Calculate the gradients of L w.r.t. to the internal parameters
    # Calculate dL/ds
    dLds = fun_dLds(s,I,t,F);
    # Calcualte dL/dF
    dLdF = fun_dLdF(s,t);
    # Calculate dL/dt
    dLdt = fun_dLdt(s,t,F);

    ## Calculate the gradients of L w.r.t. param.*
    dLdn   = dLdF * dFdn;
    dLdx_m = dLdF * dFdx_m;
    dLdx_f = dLdF * dFdx_f;

    ## Calculate the normalizing inner product
    nip = (U_adj' * dLds * U_dir);

    ds = ds_struct(1,1,1,1,1)
    ## Calculate the gradients of s w.r.t. param.*
    ds.n   = - (U_adj' * dLdn   * U_dir) / nip;
    ds.tau = - (U_adj' * dLdt   * U_dir) / nip;
    ds.x_m = - (U_adj' * dLdx_m * U_dir) / nip;
    ds.x_f = - (U_adj' * dLdx_f * U_dir) / nip;

    ## Calculate the sensitivity of s to a generic change to the operator, L
    ds.L  =  - (U_adj' * U_dir) / nip;

    return emode, ds
end
