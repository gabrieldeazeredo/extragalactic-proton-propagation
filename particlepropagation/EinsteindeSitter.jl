include("../PhysicalConsts.jl")

module EinsteindeSitter
"""
    Useful cosmology function, using EinsteindeSitterflat Universe
    model. 

    Matter density parameter = 1
    Dark energy density parameter = 0
"""
using .PhysicalConsts

export t_0, z, t, dt_dz, H

const t_0 = t(0)

"""
    Convert time to redshift.

    parameters:
        t (yr) time
"""
function z(t)
    return (t_0 / t)^(2/3) - 1
end

"""
    Convert redshift to time.

    return:
        t (yr) time
"""
function t(z)
       
    return (2 / 3 * 1 / H_0) * (1 + z)^(-3/2)
end

"""
    Compute dt/dz factor.

    parameters:
        z redshift

    return:
        dt/dz (yr)
"""
function dt_dz(z)

    return -1 / H_0 * (1 + z)^(-5/2)
end

function H(z)
    return H_0 * (1 + z)^(3/2)
end

end # EinsteindeSitter