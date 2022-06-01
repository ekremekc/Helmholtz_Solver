module nondim_param

include("./parameters.jl")
using .parameters

pp = param_dim{Float64}()


mutable struct ref
    l_dim::Float64
    p_dim::Float64
    u_dim::Float64
end

mutable struct param_struct
    x_m::Float64
    x_f::Float64
    n::Float64
    tau::Float64
    gam::Float64
    R::Float64
    a::Float64
    an::Float64
end

temp_ref = ref(1,1,1)
param_temp = param_struct(1,1,1,1,1,1,1,1)

function fun_nondim(pp::param_dim)
    # initializing with default variables


    ## Calculate reference scales
    temp_ref.l_dim = pp.X
    temp_ref.p_dim = pp.pbar
    temp_ref.u_dim = sqrt(pp.gam*pp.pbar / pp.rhobar)


    # initializing with default variables

    ## Calculate nondimensional parameters
    param_temp.x_m = pp.x_m / temp_ref.l_dim
    param_temp.x_f = pp.x_f / temp_ref.l_dim
    param_temp.n   = pp.n   / temp_ref.p_dim
    param_temp.tau = pp.tau / temp_ref.l_dim * temp_ref.u_dim
    param_temp.gam = pp.gam
    param_temp.R   = pp.R
    param_temp.a   = pp.a / temp_ref.l_dim
    param_temp.an  = pp.an / temp_ref.l_dim

    return temp_ref, param_temp
end

temp, param = fun_nondim(pp)

export temp, param

end
