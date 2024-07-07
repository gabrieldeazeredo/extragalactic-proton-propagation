include("pair_production.jl")
include("pion_production.jl")
include("adiabatic_loss.jl")

"""
    Compute 1/E * dE/dt for all process, taking account redshift 
    dependecy of cosmic microwave background photon field.

    parameters:
        E (eV) energy of proton.
        z current redshift.

    return:
        1/E * dE/dt (yr)^-1

"""
function beta(E, z)
    
    return (1 + z)^3 * beta_pair(E * (1 + z)) + 
           (1 + z)^3 * beta_pion(E * (1 + z)) + 
           (1 + z)^(3/2) * beta_rsh()
end

"""
    Shortcut function to compute unnormalized energy loss rate by pair and
    pion production, taking account redshift dependency.

    parameters:
        E (eV) energy of proton.
        z current redshift.

    return:
        b = dE/dt (eV/yr)
"""
function b(E, z)

    return (1 + z)^2 * b0(E * (1 + z))
end

"""
    Shortcut function to compute unnormalized energy loss rate by pair and
    pion production.

    parameters:
        E (eV) energy of proton.

    return:
        b = dE/dt (eV/yr)
"""
function b0(E)

    return b0_pair(E) + b0_pion(E)
end

"""
    Shortcut function to compute db/dE derivative with finite differences method.

    parameters:
        E (eV) energy of proton.
    
    return:
        db/dE (yr)^-1
"""
function db0_dE(E)

    h = 1e-6 * E

    return (b0(E + h) - b0(E - h)) / (2 * h)
end