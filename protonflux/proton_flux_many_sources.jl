include("proton_flux.jl")

function mod_factor_many_sources(E, z_max, z_min)

    println("E = ", E)

    int_dz = SciPy.integrate.quad(integrand_many_sources, z_min, z_max, args = (E))[1]

    return int_dz
end

function integrand_many_sources(z, E)
    
    gamma = 1.1

    return 1 / (1 + z)^(5/2) * 1 / (Lambda(E, z)^(gamma + 1)) * dEg_dE(z, E)
end