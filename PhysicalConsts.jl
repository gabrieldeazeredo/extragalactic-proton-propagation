"""
    Useful physical constants.
"""
module PhysicalConsts

export joule_to_ev, c, h_bar, m_e, m_p, k_B, alpha, e, T, H_0, r_0, omega_l, omega_m

const c = 2.99792458e8         # (m/s)
const h_bar = 6.6e-16          # (eV.s)
const m_e = 9.1093837015e-31   # (Kg)
const m_p = 1.67262192369e-27  # (Kg)
const k_B = 8.6e-5             # (eV/K)
const alpha = 0.0072973525693  # dimensionless
const e = 1.602176634e-19      # (C)
const T = 2.73                 # (K)
const H_0 = 75 * 1.022e-12     # (1/yr) 
const r_0 = 2.82e-15           # (m)
const omega_m = 0.3147         # dimensionless
const omega_l = 0.6853         # dimensionless

"""
    function to convert between energy units. Joule -> eV
"""
function joule_to_ev(joule)

    return joule * 6.242e18
end


end # PhysicalConsts