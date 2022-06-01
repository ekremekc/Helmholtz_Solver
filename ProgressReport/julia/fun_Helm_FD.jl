#include("./parameters.jl")
include("./fun_normalize.jl")
include("./fun_nondim.jl")
include("./fun_clencurt.jl")
include("./fun_cheb.jl")

using LinearAlgebra

#using .parameters

using .nondim_param:param # importing non-dimensionless variables

include("./structures.jl")
using .structures:emodes,ds_struct  # Importing emodes and ds structures



## MAIN FUNCTION

function fun_Helm_FD(param = param, s0 = 1*im*pi, NP = 40)

    # Helmholtz Finite Difference model, solved as a nonlinear eigenvalue problem
    # assuming p=0 at both ends
    #
    ##
    function fun_L(D2,I,F,t,s)
        L  = s^2*I + exp(-s*t) * F .- D2
        return L
    end
    ##
    function  fun_dLds(I,F,t,s)
        dLds = 2*s*I - t * exp(-s*t) * F
        return dLds
    end
    ##
    function fun_dLdF(t,s)
        dLdF = exp(-s*t)
        return dLdF
    end
    ##
    function fun_dLdt(F,t,s)
        dLdt = - s * exp(-s*t) * F
        return dLdt
    end
    # INPUTS
    # param.x_m     measurement position (nondimensional)
    # param.x_f     heat releae position (nondimensional)
    # param.n       heat release index (nondimensional)
    # param.tau 	heat release time delay (nondimensional)
    # param.gam 	ratio of specific heats
    # param.an      parameter governing the width of heat and measurement zones
    #
    # scheme.s0     starting value of s
    # scheme.NP     number of Chebyshev points
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

    # Generate the grid and differentiation and mass matrices
    # Unpack the number of Chebyshev elements (NP-1)
    N = NP-1
    # Create collocation points from x = -1 to +1 and weighting vector, m
    x_base,m_base = fun_clencurt(N)
    # Create 1st order Chebyshev differentiation matrix from -1 to 1
    D_base,~ = fun_cheb(N)
    # Convert to points from x = 0 to +1
    x = (x_base.+1.0)/2
    x = vec(x)
    m = m_base/2
    D = D_base.*2

    # Create the mass matrix
    M = Diagonal(vec(m));

    # Create the 2nd order Chebyshev differentiation matrix
    D2 = D*D;

    # Create the identity matrix
    Id = I(N+1);

    ## Calculate the internal parameters
    # Generate the heat release envelope, v(x), which integrates to 1
    v = exp.(-(x.-param.x_f).^2 ./param.an^2)/sqrt(pi)/param.an
    # Generate the measurement envelope, w(x), which integrates to 1
    w = exp.(-(x.-param.x_m).^2 ./param.an^2)/sqrt(pi)/param.an
    # Wrap v, w, gamma, param.n, M, and D into a matrix, F
    F = (v*w') * (param.gam-1)/param.gam * param.n * M * D

    # Extract tau
    t = param.tau;

    ## Iterate to find s
    # Set tolerance, initial s, and dummy dels
    tol = 1e-8; s = s0; dels = 2*tol;
    while abs(dels) > tol
        # Evaluate L and apply Dirichlet boundary conditions
        L = fun_L(D2,Id,F,t,s)
        L = L[2:N,2:N]
        # Evaluate dL/ds and apply Dirichlet boundary conditions
        dLds = fun_dLds(Id,F,t,s)
        dLds = dLds[2:N,2:N]
        # evaluate new s with Jacobi's formula
        dels = - 1/tr(L\dLds)
        # Update s
        s = s + dels
    end
    # Find the corresponding left and right eigenvectors for P
    L = fun_L(D2,Id,F,t,s);
    L = L[2:N,2:N];
    P_dir = nullspace(L);
    P_dir = [0;P_dir;0];
    P_adj = nullspace(L');
    P_adj = [0;P_adj;0];
    # Display the dimensionless eigenvalue
    println("fun_Helm_FD:  s = ", s," (nondim)")

    ## Evaluate the eigenfunction on the FD grid
    # Normalize P
    P_dir = vec(P_dir)

    P,f = fun_normalize(P_dir,M);
    # Create U from P

    U = -D*P/(param.gam*s);

    ## Wrap output into structure
    emode = emodes(s,x,P,U)
    # Wrap eigenvalue
    emode.s = s;
    # Wrap eigenfunction
    emode.x = x; emode.P = P; emode.U = U;

    ## Return if sensitivities are not required
    #if nargout == 1
    #     return
    #end

    ## Calculate the gradients of the internal parameters w.r.t. param.*
    # dv/d(param.x_f)
    dvdx_f = (2*(x.-param.x_f)/param.an^2).*v;
    # dw/d(param.x_m)
    dwdx_m = (2*(x.-param.x_m)/param.an^2).*w;
    # dF/d(param.n)
    dFdn   = F / param.n;
    # dF/d(param.x_f)
    # convert 1x1 array to scalar by TYPING FIRST in order to return 1st element of 1x1
    dFdx_f = (dvdx_f*w') * (param.gam-1)/param.gam * param.n * M * D;
    # dF/d(param.x_m)
    dFdx_m = (v*dwdx_m') * (param.gam-1)/param.gam * param.n * M * D;

    ## Calculate the gradients of L w.r.t. to the internal parameters
    # Calculate dL/ds
    dLds = fun_dLds(Id,F,t,s);
    # Calcualte dL/dF
    dLdF = fun_dLdF(t,s);
    # Calculate dL/dt
    dLdt = fun_dLdt(F,t,s);

    ## Calculate the gradients of L w.r.t. param.*
    dLdn   = dLdF * dFdn;
    dLdx_m = dLdF * dFdx_m;

    dLdx_f = dLdF * dFdx_f;

    ## Calculate the normalizing inner product
    nip = (P_adj' * dLds * P_dir);

    ds = ds_struct(1,1,1,1,1)
    ## Calculate the gradients of s w.r.t. param.*
    ds.n   = - (P_adj' * dLdn   * P_dir) / nip;
    ds.tau = - (P_adj' * dLdt   * P_dir) / nip;
    ds.x_m = - (P_adj' * dLdx_m * P_dir) / nip;

    ds.x_f = - (P_adj' * dLdx_f * P_dir) / nip;

    ## Calculate the sensitivity of s to a generic change to the operator, L
    ds.L  =  - (P_adj' * P_dir) / nip;

    return emode,ds
end

#deger1,deger2 = fun_Helm_FD()
