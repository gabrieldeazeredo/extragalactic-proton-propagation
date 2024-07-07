include("../PhysicalConsts.jl")
using .PhysicalConsts

"""
    Compute 1/E * dE/dt normalized energy loss rate by pion production without
    redshift dependency(z = 0).
    
    parameters: 
        E (eV) energy of proton.

    return:
        1/E * dE/dt (yr)^-1 energy loss rate.
"""
function beta_pion(E)

    E20 = E * 10^(-20)

    return c * 3.2408E-23 / 3.17098E-8 * exp(-4 / E20) / 13.6 * (1 + 4 / E20 + 0.5 * (4 / E20)^2)
end

"""
    Shortcut function to compute unnormalized energy loss rate by pion production.

    parameters:
        E (eV) energy of proton.

    return:
        b = dE/dt (eV/yr)
"""
function b0_pion(E)
    
    return E * beta_pion(E)
end

"""
    Shortcut function to compute db/dE derivative with finite differences method.

    parameters:
        E (eV) energy of proton.
    
    return:
        db/dE (yr)^-1
"""
function db0_pion_dE(E)

    h = 1e-6 * E

    return (b0_pion(E + h) - b0_pion(E - h)) / (2 * h)
end