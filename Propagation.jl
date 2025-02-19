"""
Contains main propagation functions.
"""
module Propagation

include("EnergyLoss.jl")
using .EnergyLoss

include("Cosmology.jl")
using .Cosmology

include("NumericalMethods.jl")
using .NumericalMethods

using DelimitedFiles

export dE_dz, evolve, propagate, propagate_powerlaw

"""
    dE_dz(z, E)

Compute evolution rate ``\\frac{dE}{dz}`` [eV] of proton with energy `E` [eV] and redshift `z`.
"""
function dE_dz(z, E)
    
    return -E * dE_Edt(E, z) * dt_dz(z)
end

"""
    evolve(E_src, z_src)

Return arrays with energy x redshift evolution until a proton with energy `E_src` [eV] and
emitted by a source on redshift `z_src` reach Earth.
"""
function evolve(E_src, z_src)

    return runge_kutta4(dE_dz, z_src, E_src, 0, 250)
end

"""
    propagate(E_src, z_src)

Return final energy when a proton with energy `E` [eV] at source on redshift `z` reach Earth (z = 0).
"""
function propagate(E_src, z_src)
    
    return evolve(E_src, z_src)[2][end]
end

"""
    propagate_powerlaw(E_min, E_max, gamma, N, z_src)

Propagates N samples of a power law E^{-gamma} spectrum from source at z_src redshift. Return specturm at source
and at Earth.
"""
function propagate_powerlaw(E_min, E_max, gamma, N, z_src)

    E_earth = zeros(N, 1)
    E_source = zeros(N, 1)

    for i in 1:N 

        E_source[i] = get_powerlaw_sample(E_min, E_max, gamma)
        E_earth[i] = propagate(E_source[i], z_src)
    end

    return E_source, E_earth 
end

"""
    propagate_from_file(fname, field, z_src)

Propagates a spectrum from a file with name `fname`. The parameter `field` is the number of collumn of spectrum
and `z_src` is source's redshift.
"""
function propagate_from_file(fname, nfield, z_src)
   
   
    data = readdlm("$fname")[nfield]
    E_earth = []
    E_source = []

    for E in data
        push!(E_earth, propagate(E, z_src))
        push!(E_source, E)
    end

    return E_source, E_earth
end

end # Propagation