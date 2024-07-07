"""
    Runge-kutta 4th method.

    parameters:
        f function to solve
        (x0, y0) boundary values
        xf final x-axis value
        h step
"""
function runge_kutta_4(f, x0, y0, xf, h)

    x = collect(range(x0, xf, h))
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