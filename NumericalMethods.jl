module NumericalMethods

using Random

export runge_kutta4, get_powerlaw_sample

"""
    runge(kutta, f, E_min, y0, xf, h)

4th order Runge-kutta method.

# Arguments
- `f`: function to solve
- `E_min`: x bondary value
- `y0`: y bondary value
- `xf`: x end value
- `h`: step
"""
function runge_kutta4(f, E_min, y0, xf, h)

    x = collect(range(E_min, xf, h))
    y = Any[y0]

    h = x[2] - x[1]

    for x_k in x

        y_k = y[end]

        k_1 = f(x_k, y_k)
        k_2 = f(x_k + h / 2, y_k + h * k_1 / 2) 
        k_3 = f(x_k + h / 2, y_k + h * k_2 / 2)
        k_4 = f(x_k + h, y_k + h * k_3)

        y_k = y_k + h / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4) 
        x_k = x_k + h

        append!(y, y_k)
    end

    return x, y[1:end-1]
end

"""
    get_powerlaw_sample(E_min, E_max, gamma)
"""
function get_powerlaw_sample(E_min, E_max, gamma)
    
    if gamma == -1
        part1 = log(E_max)
        part2 = log(E_min)

        return exp((part1 - part2) * rand() + part2)
    else 
        part1 = E_max^(gamma + 1)
        part2 = E_min^(gamma + 1)

        return ((part1 - part2) * rand() + part2)^(1 / (gamma + 1))
    end 
end

end # NumericalMethods