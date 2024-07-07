include("../energy_loss/energy_loss.jl")
include("EinsteindeSitter.jl")
include("runge_kutta_4.jl")

import .EinsteindeSitter

function f(z, E)

    return -E * beta(E, z) * EinsteindeSitter.dt_dz(z)
end

function proton_propa(E_src, z_src)

    return runge_kutta_4(f, z_src, E_src, 0, 250)
end

function inv_proton_propa(E, z_src)

    return runge_kutta_4(f, 0, E, z_src, 50)
end

function E_g(E, z_src)

    return inv_proton_propa(E, z_src)[2][end]
end
