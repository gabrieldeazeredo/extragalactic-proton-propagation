"""
Contains computations of energy losses rates that ultra-high energy protons
are subjected by cosmic radiation backgrond and Universe expansion.
"""
module EnergyLoss

include("PhysicalConsts.jl")
using .PhysicalConsts

using QuadGK

export dE_Edt, ee_dE_Edt, pi_dE_Edt, f_nu, phi

"""
    dE_Edt(E,z)

Compute total ``\\frac{-1}{E}\\frac{dE}{dt}`` [yr^-1] for a proton with energy `E` [eV] and redshift `z`.

Default processes taking into account:

- photopion production with CMB,
- photoapir production with CMB,
- adiabatic energy loss (flat Universe of Einstein de-Sitter).
"""
function dE_Edt(E, z)
    
    return (1 + z)^3 * ee_dE_Edt(E * (1 + z)) + 
           (1 + z)^3 * pi_dE_Edt(E * (1 + z)) + 
           ad_dE_Edt(z)
end

"""
    ee_dE_Edt(E)

Compute electron pair producttion loss rate ``\\frac{-1}{E}\\frac{dE}{dt}`` [yr^-1] for a proton with energy `E` [eV].

[1] BLUMENTHAL, G. R. Energy loss of high-energy cosmic rays in pair-producing collisions with ambient photons. 
Phys. Rev. D, American Physical Society, v. 1, p. 1596–1602, Mar 1970.

[2] Chodorowski, M. J.; Zdziarski, A. A.; Sikora, M. Reaction Rate and Energy-Loss Rate for Photopair Production 
by Relativistic Nuclei. apj, v. 400, p. 181, nov. 1992.
"""
function ee_dE_Edt(E)

    if E > 1E20
        # despise the loss to energies > 1e20.
        return 0
    end

    # multiplicative constant [eV/yr]
    pair_const = (alpha * r_0^2 * (joule_to_ev(m_e * c^2) * k_B * T)^2 * c)/(pi^2 * h_bar^3 * c^3) / 3.17098e-8 
    
    # proton lorentz factor
    gamma =  E / joule_to_ev(m_p * c^2)

    # dimensionless constant
    nu = joule_to_ev(m_e * c^2) / (2 * gamma * k_B * T)

    return pair_const *  1/E * f_nu(nu)
end

"""
    f_nu(nu)

Compute auxiliary function to pair production loss rate. The dimensionless variable `nu` is definined in [1].

[1] BLUMENTHAL, G. R. Energy loss of high-energy cosmic rays in pair-producing collisions with ambient photons. 
Phys. Rev. D, American Physical Society, v. 1, p. 1596–1602, Mar 1970.
"""
function f_nu(nu)

    return nu ^ 2 * quadgk(xi -> phi(xi) / (exp(xi * nu) - 1), 2, Inf)[1]
end

"""
    phi(xi)

Compute auxiliary function to pair production loss rate defined in [1].

[1] BLUMENTHAL, G. R. Energy loss of high-energy cosmic rays in pair-producing collisions with ambient photons. 
Phys. Rev. D, American Physical Society, v. 1, p. 1596–1602, Mar 1970.
"""
function phi(xi)

    if xi < 25
        return pi/12 * (xi - 2)^4 / (1 + 0.8048*(xi - 2) + 0.1459*(xi - 2)^2 + 1.137E-3*(xi - 2)^3 - 3.879E-6*(xi - 2)^4)
    else
        return xi*(-86.7 + 50.96*log(xi) - 14.45*log(xi)^2 + 8/3*log(xi)^3)/(1 - (2.910*xi^(-1) + 78.35*xi^(-2) + 1837*xi^(-3)))
    end
end

"""
    ad_dE_Edt(z)

Compute adiabatic energy loss rate \$\\frac{-1}{E}\\frac{dE}{dt}\$ [yr^-1] for a particle on redshift `z`.
"""
function ad_dE_Edt(z)

    return H_0 * sqrt(omega_m * (1 + z)^3 + omega_l)
end

"""
    pi_dE_Edt(E)

Compute pion production energy loss rate aproximation ``\\frac{11}{E}\\frac{dE}{dt}`` [yr^-1] for a energy `E` [eV] 
proton and redshift `z`.
    
[1] TAYLOR, A. et al. Ultra high energy cosmic ray, neutrino, and photon propagation and the multi-messenger approach. 
In: AIP Conference Proceedings. AIP, 2009. p.94–114.
"""
function pi_dE_Edt(E)

    A = 4e20

    return 2.25e-8 * exp(-A/E) * (1 + A/E + 0.5 * (A/E)^2)
end

end # EnergyLoss