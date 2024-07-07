include("../energyloss/energy_loss.jl")
include("EinsteindeSitter.jl")
include("runge_kutta_4.jl")

import .EinsteindeSitter

"""
    Equation defined from energy loss rate

        -1/E * /dE/dt = β(E,z) dt/dz

    used to mapping energy measured in earth, in 
    source.

    parameters:
        E (eV) proton energy
        z current redshift
"""
function f(z, E)

    return -E * beta(E, z) * EinsteindeSitter.dt_dz(z)
end

"""
    Solve equation of energy loss by Runge Kutta method
    until reach Earth (z = 0).

    parameters:
        E_src (eV) proton energy at source
        z_src redshift of source
"""
function proton_propa(E_src, z_src)

    return runge_kutta_4(f, z_src, E_src, 0, 250)
end

"""
    Solve equation of energy loss by Runge Kutta method
    until reach source's redshift.

    parameters:
        E (eV) proton energy at Earth
        z_src redshift of source
"""
function inv_proton_propa(E, z_src)

    return runge_kutta_4(f, 0, E, z_src, 50)
end

"""
    Compute E_g needed energy at source for a proton reach
    the Earth with energy E.

    parameters:
        E (eV) energy of proton at Earth
        z redshift of source
"""
function E_g(E, z_src)

    return inv_proton_propa(E, z_src)[2][end]
end
