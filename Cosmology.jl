"""
Cosmological relations of Einstein-de Sitter's flat Universe. 
"""
module Cosmology

include("PhysicalConsts.jl")
using .PhysicalConsts

export dt_dz

"""
    dt_dz(z)

Computes dt_dz in redshift `z`.
"""
function dt_dz(z)

    return -1 / H_0 * (1 + z)^(-5/2)
end

end # Cosmology