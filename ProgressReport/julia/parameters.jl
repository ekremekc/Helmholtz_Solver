

module parameters

    using Parameters

    @with_kw struct param_dim{Float64}

        rhobar::Float64 = 1.2 # kg/m3
        pbar::Float64   = 1e5 # Pascals
        gam::Float64    = 1.4 # non-dim
        X::Float64      = 1 # metres
        x_m::Float64    = 0.20 # metres
        x_f::Float64    = 0.25 # metres
        n::Float64      = 100000 # [h.rel per unit vol / velocity]
        tau::Float64    = 0.001 # seconds
        R::Float64      = -1.0  # reflection coefficient at both ends
        a::Float64      = 0.05 # radius of tube, metres
        an::Float64     = 0.01 # width of heat release and measurement envelopes
    end

    export param_dim

end

# pp = param_dim{Float64}()
# other = param_dim{Float64}(rhobar = 2.0)
# rho = pp.rhobar
# println(rho)
# println(other.rhobar)
