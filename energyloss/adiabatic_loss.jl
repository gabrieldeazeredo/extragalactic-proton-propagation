include("../PhysicalConsts.jl")
using .PhysicalConsts

"""
    Compute 1/E * dE/dt for energy loss by an expanding Universe.

    return:
        1/E * dE/dt (yr)^-1
"""
function beta_rsh()

    return H_0
end

"""
    Shortcut function to computes unnormalized energy loss rate by pair production.

    parameters:
        E (eV) energy of proton.

    return:
        b = dE/dt (eV/yr)
"""
function b0_rsh(E)
    
    return E * beta_rsh()
end

"""
    Shortcut function to compute db/dE derivative. Analytic db/dE = H(z = 0).
    
    return:
        db/dE (yr)^-1
"""
function db0_rsh_dE()
    
    return beta_rsh()
end