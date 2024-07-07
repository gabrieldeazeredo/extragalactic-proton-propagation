include("../PhysicalConsts.jl")
using .PhysicalConsts, QuadGK

"""
    Compute 1/E * dE/dt normalized energy loss rate by electron pair production without
    redshift dependency (z = 0).
    
    parameters: 
        E (eV) energy of proton.

    return:
        1/E * dE/dt (yr)^-1 energy loss rate.
"""
function beta_pair(E)

    if E > 1E23
        # despise from energies > 1e23.
        return 0
    end

    #   Pair-production energy rate constant. [eV/yr]
    pair_const = (alpha * r_0^2 * (joule_to_ev(m_e * c^2) * k_B * T)^2 * c)/
    (pi^2 * h_bar^3 * c^3) / 3.17098e-8 
    
    lorentz_factor =  E / joule_to_ev(m_p * c^2)
    nu = joule_to_ev(m_e * c^2) / (2 * lorentz_factor * k_B * T)

    return pair_const * f_nu(nu) * 1 / E
end

"""
    Compute auxiliar function defined by dimensionless variable (nu). 
"""
function f_nu(nu)

    return nu ^ 2 * quadgk(xi -> phi(xi) / (exp(xi * nu) - 1), 2, Inf)[1]
end

"""
    Shortcut function to computes unnormalized energy loss rate by pair production.

    parameters:
        E (eV) energy of proton.

    return:
        b = dE/dt (eV/yr)
"""
function b0_pair(E)

    return E * beta_pair(E)
end

"""
    Shortcut function to compute db/dE derivative with finite differences method.

    parameters:
        E (eV) energy of proton.
    
    return:
        db/dE (yr)^-1
"""
function db0_pair_dE(E)

    h = 1e-6 * E
    return (b0_pair(E + h) - b0_pair(E - h)) / (2 * h)
end

"""
    Compute auxiliar function phi(xi).
"""
function phi(xi)

    if xi < 25
        return pi/12 * (xi - 2)^4 / (1 + 0.8048*(xi - 2) + 0.1459*(xi - 2)^2 + 1.137E-3*(xi - 2)^3 - 3.879E-6*(xi - 2)^4)
    else
        return xi*(-86.7 + 50.96*log(xi) - 14.45*log(xi)^2 + 8/3*log(xi)^3)/(1 - (2.910*xi^(-1) + 78.35*xi^(-2) + 1837*xi^(-3)))
    end
end