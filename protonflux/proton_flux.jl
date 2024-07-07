include("../energy_loss/energy_loss.jl")
include("../particle_propagation/numerical_propagation.jl")

using SciPy

function mod_factor(E, z_g)

    gamma = 1.1

    return 1 / (sqrt(1 + z_g) - 1)^2 * 1 / (Lambda(E, z_g)^(gamma + 1)) * dEg_dE(z_g, E)
end

function dEg_dE(z_g, E)

    int_dz = SciPy.integrate.quad(integrand_dz, 0, z_g, args = (E))[1]

    return (1 + z_g) * exp(1 / H_0 * int_dz)
end

function integrand_dz(z, E)

    return (1 + z)^(1/2) * db0_dE(E_g(E, z) * (1 + z))
end

function Lambda(E, z_g)

    return E_g(E, z_g) / E
end