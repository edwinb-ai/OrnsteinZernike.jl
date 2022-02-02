function make_fft_plan(x)
    xc = zeros(2 * length(x))
    p = plan_rfft(xc)

    return (xc, p)
end

function _fft_with_plan(x, r, p)
    n = length(x)
    r[1:n] = x
    result = imag(conj!(p * r))

    return result[1:(end - 1)]
end

function fft_oz(f, r, rmax, nr, p)
    for i in 1:nr
        f[i] *= (i - 1)
    end

    result = _fft_with_plan(f, r, p)

    normalization = (4.0 * rmax^3) / nr^2
    for i in 2:nr
        result[i] *= normalization / (i - 1)
    end

    return result
end

function ifft_oz(f, r, rmax, nr, p)
    for i in 1:nr
        f[i] *= (i - 1)
    end

    result = _fft_with_plan(f, r, p)

    normalization = nr * (1.0 / 2.0 / rmax^3)
    for i in 2:nr
        result[i] *= normalization / (i - 1)
    end

    return result
end
