function oz_solve(st::Structure, γ, br, params::Parameters)
    ur = apply_potential(st.pot, st.r) ./ params.ktemp
    broadcast!(x -> isnan(x) ? 0.0 : x, ur, ur)
    cr = closure_relation(γ, br, ur, st.r)
    broadcast!(x -> isnan(x) ? 0.0 : x, cr, cr)

    fits = continuous!(cr, st.r, params.diam)

    r, p = make_fft_plan(cr)
    ck = fft_oz(cr, r, params.rmax, params.mr, p)
    remove_continuous!(ck, st.q, params.diam, fits)

    δ = @. 1.0 - params.δρ * ck
    γk = @. params.δρ * ck^2 / δ

    r, p = make_fft_plan(γk)
    new_gamma = ifft_oz(γk, r, params.rmax, params.mr, p)

    cr = closure_relation(new_gamma, br, ur, st.r)
    broadcast!(x -> isnan(x) ? 0.0 : x, cr, cr)

    return cr, new_gamma
end

function build_solve_pairwise(x::AbstractMatrix, y::AbstractVector, r)
    pairwise_product_matrix = pairwise_inner_product(x, r)

    n = size(pairwise_product_matrix, 1)
    difference_product = zeros(eltype(x), n)

    for i in 1:n
        temp_product = y .* x[i, :]
        difference_product[i] = inner_product(temp_product, r)
    end

    # Solve the linear system
    return pairwise_product_matrix \ difference_product
end

function create_solution(coef::AbstractVector, gmat::AbstractMatrix)
    n = size(gmat, 1)
    result = (1.0 - sum(coef)) .* gmat[end, :]

    for i in 1:(n - 1)
        @. result += coef[i] * gmat[n - i, :]
    end

    return result
end

function check_precision(t::Interaction, y::AbstractVector, r)
    cr, g = oz_solve(t.structure, y, t.bridge, t.params)
    diff_vector = g .- y
    prec = inner_product(diff_vector .^ 2, r)

    return (√prec, diff_vector, g, cr)
end

function differences_coefs!(cmat, dmat, y)
    for i in 1:(size(cmat, 1))
        cmat[i, :] = y .- dmat[end - i, :]
    end
end

function solve_to_precision(t::Interaction, cmat, dmat, gmat; prec=1e-3)
    new_prec = Inf
    view_gmat = gmat[2:end, :]
    last_diff = view_gmat[end, :] .- view_gmat[end - 1, :]
    differences_coefs!(cmat, dmat, last_diff)
    cr = zeros(eltype(gmat), size(gmat[1, :]))

    while (new_prec - prec) > 0
        coefs = build_solve_pairwise(cmat, last_diff, t.structure.r)
        fnew = create_solution(coefs, view_gmat)
        new_prec, diff_vector, g, cr = check_precision(t, fnew, t.structure.r)
        dmat[end, :] = diff_vector
        view_gmat[end, :] = g
        if (new_prec - prec) < 0 # early stopping
            break
        end
        view_gmat = circshift(view_gmat, (-1, 0))
        dmat = circshift(dmat, (-1, 0))
        cr, view_gmat[end, :] = oz_solve(
            t.structure, view_gmat[end - 1, :], t.bridge, t.params
        )
        last_diff = view_gmat[end, :] .- view_gmat[end - 1, :]
        differences_coefs!(cmat, dmat, last_diff)
    end

    gmat[2:end, :] = view_gmat

    return cr
end

function ng_method_first!(t::Interaction, γ_matrix)
    for i in 2:size(γ_matrix, 1)
        _, γ_matrix[i, :] = oz_solve(t.structure, γ_matrix[i - 1, :], t.bridge, t.params)
    end

    γ_matrix[1, :] = γ_matrix[end, :]

    return nothing
end

function ng_method!(t::Interaction, γ_matrix)
    for i in 2:size(γ_matrix, 1)
        _, γ_matrix[i, :] = oz_solve(t.structure, γ_matrix[i - 1, :], t.bridge, t.params)
    end
    dif_matrix = diff(γ_matrix; dims=1)

    n, m = size(dif_matrix)
    coeff_matrix = zeros(n - 1, m)

    cr = solve_to_precision(t, coeff_matrix, dif_matrix, γ_matrix)

    return cr
end

function solve(t::Interaction; npoints=5)
    γ_matrix = zeros(npoints + 2, t.params.mr) # +2 for initial and one last computation
    temp_gamma = zeros(t.params.mr)

    δρ = t.params.ρ / t.params.nρ
    rho_range = range(δρ, t.params.nρ * δρ; length=t.params.nρ)

    # First iteration
    t.params.δρ = rho_range[1]
    ng_method_first!(t, γ_matrix)
    temp_gamma = γ_matrix[end, :]

    r = Result()

    for j in rho_range[2:end]
        t.params.δρ = j

        cr = ng_method!(t, γ_matrix)
        r = Result(cr, γ_matrix[end, :], t.structure, t.params)

        extrapolated = extrapolation(temp_gamma, γ_matrix[end, :], t.params.ρ, j)
        γ_matrix[1, :] = extrapolated
        temp_gamma = γ_matrix[end, :]
    end

    return r
end
