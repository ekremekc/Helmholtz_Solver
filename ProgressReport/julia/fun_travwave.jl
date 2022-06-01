#include("./parameters.jl")
include("./fun_normalize.jl")
include("./fun_nondim.jl")

using LinearAlgebra

#using .parameters

using .nondim_param:param # importing non-dimensionless variables

include("./structures.jl")
using .structures:emodes, ds_struct  # Importing emodes and ds structures


## MAIN FUNCTION

function fun_travwave(param = param, scheme = 1*im*pi)

    ## SUBFUNCTIONS

    ## Evaluate L
    function fun_L(s,h,tu,td,tmf,t,R)

        L    = [ 1+R*exp(-s*tu)                                                      -1-R*exp(-s*td)   ;
             1-R*exp(-s*tu)+h*exp(-s*t)*(exp(-s*tmf)-R*exp(+s*tmf)*exp(-s*tu))   (1-R*exp(-s*td)) ]
        return L
    end

    ## Evaluate dL/ds
    function fun_dLds(s,h,tu,td,tmf,t,R)

        dLds = [-R*tu*exp(-s*tu)  R*td*exp(-s*td);
             R*tu*exp(-s*tu)-h*(exp(-s*tmf)-R*exp(+s*tmf)*exp(-s*tu))*t*exp(-s*t)+h*(-tmf*exp(-s*tmf)-R*tmf*exp(+s*tmf)*exp(-s*tu) +R*tu*exp(+s*tmf)*exp(-s*tu))*exp(-s*t) R*td*exp(-s*td)]
        return dLds
    end

    ## Evaluate dL/d(t)

    function fun_dLdt(s,h,tu,tmf,t,R)

        dLdt = [ 0   0 ;
                 -s*h*(exp(-s*tmf)-R*exp(+s*tmf)*exp(-s*tu))*exp(-s*t)   0]
        return dLdt
    end

    ## Evaluate dL/d(h)
    function fun_dLdh(s,h,tu,td,tmf,t,R)
        dLdh = [0  0;
                (exp(-s*tmf)-R*exp(+s*tmf)*exp(-s*tu))*exp(-s*t)   0]
        return dLdh
    end

    ## Evaluate dL/d(tu)
    function fun_dLdtu(s,h,tu,td,tmf,t,R)

        dLdtu    = [ -R*s*exp(-s*tu)                                               0;
                     +R*s*exp(-s*tu) + h*(+R*s*exp(+s*tmf)*exp(-s*tu))*exp(-s*t)   0]
        return dLdtu
    end

    ## Evaluate dL/d(td)
    function fun_dLdtd(s,h,tu,td,tmf,t,R)

        dLdtd    = [ 0   R*s*exp(-s*td) ;
                     0   R*s*exp(-s*td)]
        return dLdtd
    end

    ## Evaluate dL/d(tmf)
    function fun_dLdtmf(s,h,tu,td,tmf,t,R)

        dLdtmf    = [ 0   0 ;
                      h*(-s*exp(-s*tmf)-R*s*exp(+s*tmf)*exp(-s*tu))*exp(-s*t)   0]
        return dLdtmf
    end

    # Start of MAIN FUNCTION

    # Calculate the internal parameters
    # h = (gamma-1)/gamma * n
    h = (param.gam-1)/param.gam*param.n
    # time delay
    t = param.tau
    # tu, time for wave to travel flame -> upstream boundary -> flame
    tu = 2*(param.x_f)
    # td, time for wave to travel flame -> downstream boundary -> flame
    td = 2*(1-param.x_f)
    # tmf, time for wave to travel from flame -> measurement point
    tmf = (param.x_f - param.x_m)
    # R, reflection coefficient, assumed to be the same at both ends
    R = param.R

    if tmf < 0
        println("For the network model, x_m must be less than x_f")
        return
    end

    tol = 1e-8
    s = scheme
    dels = 2*tol

    while abs(dels)>tol
        # Evaluate L
        L = fun_L(s,h,tu,td,tmf,t,R)
        # Evaluate dL/ds
        dLds = fun_dLds(s,h,tu,td,tmf,t,R)
        # evaluate new s with Jacobi's formula
        dels = - 1/tr(L\dLds)
        # Update s
        s = s + dels
    end

    # Evaluate the right eigenvector, q_dir, s.t. L*q_dir = [0;0]
    q_dir = [+1+R*exp(-s*td) ; +1+R*exp(-s*tu)]
    # Evaluate the left eigenvector, q_adj, s.t. q_adj'*L = [0,0]
    q_adj = conj([+1-R*exp(-s*td) ; +1+R*exp(-s*td)])
    # Display the dimensionless eigenvalue
    println("fun_travwave: s = ", s," (nondim)")

    ## Evaluate the eigenfunction on a user-defined grid
    # Set the number of elements on which to evaluate the eigenfunction
    Nx = 1000;
    # Set the abscissa
    x = range(0,1,length = Nx+1)
    # Calculate the pressure and velocity eigenfunctions upstream of x_f
    Pu = (1+R*exp(-s*td))*(R*exp(-s*tu)*exp.(-s*(x.-param.x_f)) + exp.(+s*(x.-param.x_f)))
    Uu = (1+R*exp(-s*td))*(R*exp(-s*tu)*exp.(-s*(x.-param.x_f)) - exp.(+s*(x.-param.x_f)))/(param.gam)
    #Calculate the pressure and velocity eigenfunctions downstream of x_f
    Pd = (1+R*exp(-s*tu))*( exp.(-s*(x.-param.x_f)) + R*exp(-s*td)*exp.(+s*(x.-param.x_f)))
    Ud = (1+R*exp(-s*tu))*( exp.(-s*(x.-param.x_f)) - R*exp(-s*td)*exp.(+s*(x.-param.x_f)))/(param.gam)
    # Stitch together the pressure and velocity eigenfunctions
    P = Pu.*(x.<param.x_f) + Pd.*(x.>=param.x_f)

    U = Uu.*(x.<param.x_f) + Ud.*(x.>=param.x_f)
    #Generate the mass matrix on this abscissa
    dx = x[2]-x[1]
    M = dx*I(Nx+1)
    M[1,1] = M[1,1]/2
    M[Nx+1,Nx+1] = M[Nx+1,Nx+1]/2
    # Normalize P
    P,f = fun_normalize(P,M)
    # Normalize U by the same factor
    U = U/f

    emode = emodes(s,x,P,U)

    ## Return if sensitivities are not required
    ## if nargout == 1; return; end

    # Calculate the gradients of L w.r.t. to the internal parameters
    # Calculate dL/ds
    dLds   = fun_dLds(s,h,tu,td,tmf,t,R)
    # Calculate dL/dt
    dLdt   = fun_dLdt(s,h,tu,tmf,t,R)
    # Calcualte dL/dh
    dLdh   = fun_dLdh(s,h,tu,td,tmf,t,R)
    # Calculate dL/dtu
    dLdtu  = fun_dLdtu(s,h,tu,td,tmf,t,R)
    # Calculate dL/dtd
    dLdtd  = fun_dLdtd(s,h,tu,td,tmf,t,R)
    # Calculate dL/dtmf
    dLdtmf = fun_dLdtmf(s,h,tu,td,tmf,t,R)

    # Calculate the gradients of the internal parameters w.r.t. param.*
    # d(tu)/d(param.x_f)
    dtudx_f =  2
    # d(td)/d(param.x_f)
    dtddx_f = -2
    # d(tmf)/d(param.x_f)
    dtmfdx_f = +1
    # d(tmf)/d(param.x_m)
    dtmfdx_m = -1
    # d(h)/d(param.n)
    dhdn = (param.gam-1)/param.gam

    ## Calculate the gradients of L w.r.t. param.*
    dLdn   = dhdn * dLdh
    dLdx_m = dtmfdx_m * dLdtmf
    dLdx_f = dtudx_f * dLdtu + dtddx_f * dLdtd + dtmfdx_f * dLdtmf

    ## Calculate the normalizing inner product
    nip = (q_adj' * dLds * q_dir)

    ds = ds_struct(1,1,1,1,1)
    # Calculate the gradients of s w.r.t. param.*
    ds.n   = - (q_adj' * dLdn   * q_dir) / nip;
    ds.tau = - (q_adj' * dLdt   * q_dir) / nip;
    ds.x_m = - (q_adj' * dLdx_m * q_dir) / nip;
    ds.x_f = - (q_adj' * dLdx_f * q_dir) / nip;

    # Calculate the sensitivity to a generic change to the operator, L
    ds.L  =  - (q_adj' * q_dir) / nip;

    return emode, ds
end
