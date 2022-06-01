include("./fun_normalize.jl")
include("./fun_nondim.jl")

using LinearAlgebra
using DSP
#using .parameters

using .nondim_param:param # importing non-dimensionless variables

include("./structures.jl")
using .structures:emodes,ds_struct  # Importing emodes and ds structures


## MAIN FUNCTION
function fun_Helm_FE(param = param, s0 = 1*im*pi, NE = 40)

    # Helmholtz Finite Element model, solved as a nonlinear eigenvalue problem
    # assuming p=0 at both ends
    #
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
    # emode.x1      abscissa for P1 functions
    # emode.x0      abscissa for P0 functions
    # emode.P       direct pressure eigenfunction
    # emode.U       direct velocity eigenfunction
    #
    # ds.x_m        d(s)/d(param.x_m)
    # ds.x_f        d(s)/d(param.x_f)
    # ds.n       	d(s)/d(param.n)
    # ds.tau 		d(s)/d(param.tau)
    # ds.L          d(s)/d(L)

    ## SUBFUNCTIONS
    function fun_L(N,D,M00,M11,t,F,s)

        A = exp(-s*t) * M11 * F + D'*M00*D; C = M11;
        # Apply Homogenous Dirichlet boundary conditions
        Abc = fun_Abc(N,A);
        Cbc = fun_Cbc(N,C);
        # Generate L
        L = s^2 * Cbc + Abc;
        return L
    end

    function fun_dLds(N,M11,t,F,s)
        dAds = -t * exp(-s*t) * M11 * F; C = M11;
        # Apply Homogenous Dirichlet boundary conditions
        dAdsbc = fun_Cbc(N,dAds);
        Cbc = fun_Cbc(N,C);
        # Generate dLds
        dLds = 2*s*Cbc + dAdsbc;
        return dLds
    end

    function fun_dLdF(M11,t,s)
        dLdF = M11*exp(-s*t)
        return dLdF
    end

    function fun_dLdt(M11,t,F,s)
        dLdt = - s*M11*F*exp(-s*t)
        return  dLdt
    end


    # Write zeros to the border of A and 1 in the top-left and bott-right elems
    function fun_Abc(N,A)
        Abc = A;
        Abc[N+1,:] = zeros(1,N+1);
        Abc[:,N+1] = zeros(N+1,1);
        Abc[N+1,N+1] = 1;
        Abc[1,:] = zeros(1,N+1);
        Abc[:,1] = zeros(N+1,1);
        Abc[1,1] = 1;
        return Abc
    end

    # Write zeros to the border of C
    function fun_Cbc(N,C)
        Cbc = C;
        Cbc[N+1,:] = zeros(1,N+1);
        Cbc[:,N+1] = zeros(N+1,1);
        Cbc[N+1,N+1] = 0;
        Cbc[1,:] = zeros(1,N+1);
        Cbc[:,1] = zeros(N+1,1);
        Cbc[1,1] = 0;
        return Cbc
    end

    ## MAIN FUNCTION

    ## Generate the grid and differentiation and mass matrices
    # Unpack the number of elements (NE)
    N = NE
    # Generate the positions of the collocation points for P1 functions
    x1 = range(+1,-1,length =N+1)
    dx = 2/N
    # Generate the positions of the collocation points for P0 functions
    ## 2D CONVOLUTION ALGORITHM
    x0 =  zeros(N)
    for i=1:1:length(x1)
        if i ==length(x1)
            break
        end
        x0[i] = 0.5*(x1[i]+x1[i+1])
    end

    #x0 = conv(x1,[0.5;0.5]) #,"valid"
    # Generate the mass matrix for two P1 functions
    M11 = 4*I(N+1) + diagm(1 => vec(ones(1,N))) + diagm(-1 => vec(ones(1,N)))
    M11[1,1] = 2
    M11[N+1,N+1] = 2
    M11 = M11*(dx/6)
    # Generate the mass matrix for two P0 functions
    M00 = I(N)*dx
    # Generate the mass matrix for a P0 function * a P1 function
    # M01 = eye(N+1) + diag(ones(1,N),+1); M01 = M01(1:N,:); M01 = M01*(dx/2);
    # Generate the mass matrix for a P1 function * a P0 function
    # M10 = M01';
    # Generate the difference matrix for a P1 function (makes a P0 function)
    D = I(N+1) - diagm(+1 => vec(ones(1,N)))
    D = D[1:N,:]
    D = D/dx
    # Generate the mean matrix for a P1 function (makes a P0 function)
    # G01 = eye(N+1) + diag(ones(1,N),+1); G01 = G01(1:N,:)/2;
    # Convert to points from x = 0 to +1
    x1 = (x1.+1.0)/2
    x0 = (x0.+1.0)/2
    M11 = M11/2
    M00 = M00/2
    D = D*2

    # M01 = M01/2; M10 = M10/2;

    ## Calculate the internal parameters
    # Generate the heat release envelope, v(x), which integrates to 1
    v = exp.(-(x1.-param.x_f).^2 ./param.an^2)/sqrt(pi)/param.an;
    # Generate the measurement envelope, w(x), which integrates to 1
    w = exp.(-(x0.-param.x_m).^2 ./param.an^2)/sqrt(pi)/param.an;
    # Wrap v, w, gamma, param.n into a matrix, F
    F = (v*w') * (param.gam-1)/param.gam * param.n * M00 * D;
    # Extract tau
    t = param.tau;

    ## Iterate to find s
    # Set tolerance, initial s, and dummy dels
    tol = 1e-8; s = s0; dels = 2*tol;
    while abs(dels) > tol
        # Evaluate L
        L = fun_L(N,D,M00,M11,t,F,s);
        # Evaluate dLds
        dLds = fun_dLds(N,M11,t,F,s);
        # evaluate new s with Jacobi's formula
        dels = - 1/tr(L\dLds);
        # Update s
        s = s + dels;
    end
    # Find the corresponding left and right eigenvectors for P
    L = fun_L(N,D,M00,M11,t,F,s)
    P_dir = nullspace(L)
    P_adj = nullspace(L')
    # Display the dimensionless eigenvalue
    println("fun_Helm_FE:  s = ",s," (nondim)")

    ## Evaluate the eigenfunction on the FD grid
    # Normalize P
    P_dir = vec(P_dir)
    P,temp_norm = fun_normalize(P_dir,M11);
    # Create U (P0 function) from P (P1 function)
    U = -D*P/(param.gam*s);
    #println(size(U),size(x0))
    ## Wrap output into structure
    emode = emodes(s,(x0,x1),P,U)
    # Wrap eigenvalue
    #emode.s = s;
    # Wrap eigenfunction
    #emode.x1 = x1; emode.x0 = x0; emode.P = P; emode.U = U;

    ## Return if sensitivities are not required
    #if nargout == 1; return; end

    ## Calculate the gradients of the internal parameters w.r.t. param.*
    # dv/d(param.x_f)
    dvdx_f = (2*(x1.-param.x_f)/param.an^2).*v;
    # dw/d(param.x_m)
    dwdx_m = (2*(x0.-param.x_m)/param.an^2).*w;
    # dF/d(param.n)
    dFdn   = F / param.n;
    # dF/d(param.x_f)
    dFdx_f = (dvdx_f*w') * (param.gam-1)/param.gam * param.n * M00 * D;
    # dF/d(param.x_m)
    dFdx_m = (v*dwdx_m') * (param.gam-1)/param.gam * param.n * M00 * D;

    ## Calculate the gradients of L w.r.t. to the internal parameters
    # Calculate dL/ds
    dLds = fun_dLds(N,M11,t,F,s);
    # Calculate dL/dF
    dLdF = fun_dLdF(M11,t,s);
    # Calculate dL/dt
    dLdt = fun_dLdt(M11,t,F,s);

    ## Calculate the gradients of L w.r.t. param.*
    # Work out dLd(param.n)
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

    ## Calculate the sensitivity to a generic change to the operator, L
    ds.L  =  - (P_adj' * P_dir) / nip;

    # # FEEDBACK SENSITIVITIES - DO THIS IF I HAVE TIME, BUT BASE ON MY NEW PAPER, NOT ON THE ARFM PAPER
    # # ## Calculate the sensitivity to changes in the v, tau, and w fields
    # # # These are defined such that conj(ds.[-]) is the influence on s of
    # # # increasing v(x), w(x), and t(x) at position x.
    # # #
    # # # define q =    (gam-1)/gam*n*exp(-s*t)
    # # # then   A = A0 - M11*q*(w.'*M00*D)
    # # # and:
    # # # dq/dn    =    (gam-1)/gam*exp(-s*t)
    # # # dq/dt    = -s*(gam-1)/gam*n*exp(-s*t)
    # # # and:
    # # # ds/dn = ds/dq * dq/dn
    # # # ds/dt = ds/dq * dq/dt
    # # dsdq = - (M11'*P_adj) * conj(w.'*M00*D*P_dir/nip);
    # # ds.nvec = dsdq .* conj((param.gam-1)/param.gam*exp(-s*t));
    # # ds.tvec = dsdq .* conj(-s*(param.gam-1)/param.gam*n.*exp(-s*t));
    # # ds.wvec = - (((param.gam-1)/param.gam* n*exp(-s*t))'*M11'*P_adj) * conj(M00*D*P_dir/nip);
    # #
    # # # NEW METHOD:
    # # # dL/dv(x) = - M11*(dF/dv)*exp(-s*t)*M00*D;
    # # # dL/dw(x) = - M11*(dF/dw)*exp(-s*t)*M00*D;
    # # # dL/dt(x) = - M11*F*(-s)*exp(-s*t)*M00*D;
    # #
    # #
    # # ## Calculate the feedback sensitivities
    # # # These are defined such that conj(ds.[-]) is the influence on s of
    # # # feedback from the [pressure/velocity] into the [mass/momentum/energy]
    # # # equation
    # # ds.mp = - conj(-s)*M11'*P_adj .* conj(P_dir/nip);             # from pressure into mass equation
    # # ds.mu = - M10'*P_adj .* conj(D*P_dir/nip);                    # from velocity into mass equation
    # # ds.fp = + M00'*D*P_adj .* conj(G01*P_dir/nip);                # from pressure into momentum equation
    # # ds.fu = + M00'*D*P_adj .* conj(D*P_dir/(-s*nip));             # from velocity into momentum equation
    # # ds.qp = - conj(-s*(gam-1)/gam)*M11'*P_adj .* conj(P_dir/nip); # from pressure into energy equation
    # # ds.qu = - conj((gam-1)/gam)*M10'*P_adj .* conj(D*P_dir/nip);  # from velocity into energy equation
    # #
    # # ## Calculate the receptivities of the energy and momentum equations
    # # ds.q = - conj(-s*(gam-1)/gam)*M11'*P_adj * (1/nip);
    # # ds.f = + M00'*D*P_adj * (1/nip);
    # #
    # # ## Wrap up M00 and M11, which are required to plot sensitivities per unit distance, rather than per gridpoint value
    # # ds.M00 = M00; ds.M11 = M11;
        return emode,ds
end
