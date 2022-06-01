module structures

mutable struct emodes
    s
    x
    P
    U
end

mutable struct ds_struct
    n
    tau
    x_m
    x_f
    L
end

mutable struct scheme
    s0
    Np

end


export emodes, ds_struct

end
