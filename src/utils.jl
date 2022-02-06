function continuous!(f, r, diam)
    dr = r[2] - r[1]
    m = 0

    for (i, j) in enumerate(r)
        if j < diam
            m = i
        end
    end

    fpd = -274.0 * f[m + 1] + 600.0 * f[m + 2] - 600.0 * f[m + 3]
    fpd += 400.0 * f[m + 4] - 150.0 * f[m + 5] + 24.0 * f[m + 6]
    fpd /= 120.0 * dr

    fppd = 225.0 * f[m + 1] - 770.0 * f[m + 2] + 1070.0 * f[m + 3]
    fppd += -780.0 * f[m + 4] + 305.0 * f[m + 5] - 50.0 * f[m + 6]
    fppd /= 60.0 * dr^2

    fd = fpd * (diam - r[m + 1]) + f[m + 1]

    fpi = 274.0 * f[m] - 600.0 * f[m - 1] + 600.0 * f[m - 2]
    fpi += -400.0 * f[m - 3] + 150.0 * f[m - 4] - 24.0 * f[m - 5]
    fpi /= 120.0 * dr

    fppi = 225.0 * f[m] - 770.0 * f[m - 1] + 1070.0 * f[m - 2]
    fppi += -780.0 * f[m - 3] + 305.0 * f[m - 4] - 50.0 * f[m - 5]
    fppi /= 60.0 * dr^2

    fi = fpi * (diam - r[m]) + f[m]
    delf = fd - fi
    delfp = fpd - fpi
    delfpp = fppd - fppi
    a0 = delf - diam * delfp + diam^2 * delfpp / 2.0
    a1 = delfp - diam * delfpp
    a2 = delfpp / 2.0

    for i in 1:m
        f[i] += a0 + a1 * r[i] + a2 * r[i]^2
    end

    return a0, a1, a2
end

function remove_continuous!(f, q, diam, params)
    (a0, a1, a2) = params

    for i in eachindex(f)
        if q[i] > 0.0
            part1 = a0 + 2.0 * diam * a1 + 3.0 * diam^2 * a2
            part1 -= 6.0 * a2 / q[i]^2
            part1 *= sin(q[i] * diam) / q[i]^3

            part2 = -diam * a0 + 2.0 * a1 / q[i]^2 - diam^2 * a1
            part2 += 6.0 * diam * a2 / q[i]^2
            part2 = (part2 - diam^3 * a2) * cos(q[i] * diam) / q[i]^2

            part3 = -2.0 * a1 / q[i]^4

            f[i] -= 4.0 * π * (part1 + part2 + part3)
        else
            f[i] = 0.0
        end
    end

    return nothing
end

@inline function extrapolation(g, f, rho, drho)
    a = (f .- g) ./ drho
    b = @. f - a * rho
    result = @. a * (rho + drho) + b

    return result
end

function inner_product(f, df)
    (integral, _) = romberg(df, f)

    return integral
end

function pairwise_inner_product(x::AbstractMatrix, dx::AbstractRange)
    n = size(x, 1)
    result_matrix = zeros(eltype(x), (n, n))

    # First, deal with the diagonal
    for i in 1:n
        temp_prod = x[i, :] .* x[i, :]
        result_matrix[i, i] = inner_product(temp_prod, dx)
    end

    for i in 1:n
        for j in (i + 1):n
            temp_prod = x[i, :] .* x[j, :]
            result_matrix[i, j] = inner_product(temp_prod, dx)
            result_matrix[j, i] = result_matrix[i, j] # Fill in the symmetric part
        end
    end

    return result_matrix
end
