abstract type Closure end

struct HypernettedChain <: Closure end

struct PercusYevick <: Closure end

struct ModifiedVerlet <: Closure end

struct MeanSpherical <: Closure end

struct SoftMeanSpherical <: Closure end

struct HMSA <: Closure
    α::Float64
end

struct ConstantClosure <: Closure
    a::Real
end

struct Parameters{T}
    diam::T
    rmax::T
    ρ::T
    mr::Integer
    nρ::Integer
    δρ::T
    ktemp::T
end

function Parameters(diam, rmax, ρ, mr::Integer, nr, ktemp)
    d, rm, r, n = promote(diam, rmax, ρ, nr)

    return Parameters{typeof(d)}(d, rm, r, mr, n, zero(typeof(d)), ktemp)
end

struct Structure{T<:AbstractArray}
    r::T
    q::T
    pot::Potential
end

function Structure(p::Parameters, pot::Potential)
    dr = p.rmax / p.mr
    dq = π / p.rmax
    r = collect(range(0.0, (p.mr - 1) * dr; length=p.mr))
    q = collect(range(0.0, (p.mr - 1) * dq; length=p.mr))

    return Structure{typeof(r)}(r, q, pot)
end

struct Result{T<:AbstractArray}
    r::T
    q::T
    gr::T
    cr::T
    sq::T
end

function Result()
    vecs = [zeros(2) for i in 1:5]
    return Result{Vector{Float64}}(vecs...)
end

function Result(cr, gam, st::Structure, params::Parameters)
    # First, we deal with the RDF
    gr = @. gam[2:end] + cr[2:end] + 1.0
    result_cr = copy(cr)

    # Next, we deal with the structure factor
    # fits = continuous!(cr, st.r, params.diam)

    r, p = make_fft_plan(cr)
    ck = fft_oz(cr, r, params.rmax, params.mr, p)
    # remove_continuous!(ck, st.q, params.diam, fits)

    sq_inv = @. 1.0 - params.ρ * ck[2:end]

    # Also save the new structure information
    new_r = collect(st.r[2:end])
    new_q = collect(st.q[2:end])

    return Result{typeof(gr)}(new_r, new_q, gr, result_cr, 1.0 ./ sq_inv)
end

struct Interaction
    params::Parameters
    structure::Structure
    bridge::Closure

    Interaction(p, s, b) = new(p, s, b)
end
