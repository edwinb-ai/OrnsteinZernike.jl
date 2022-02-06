abstract type Closure end

# ! FIXME: Add a way to use external functions as a closure relation.

@inline closure_relation(γ, br, u, r) = exp.(-u .+ γ .+ br) .- γ .- 1.0
struct HypernettedChain <: Closure end

struct PercusYevick <: Closure end

struct ModifiedVerlet <: Closure end

struct MeanSpherical <: Closure end

struct SoftMeanSpherical <: Closure end

struct HMSA{T<:Real} <: Closure
    α::T
end

struct ConstantClosure{T<:Real} <: Closure
    a::T
end

@inline closure_relation(γ, ::HypernettedChain, u, r) = closure_relation(γ, 0.0, u, r)
@inline closure_relation(γ, c::ConstantClosure, u, r) = closure_relation(γ, c.a, u, r)
function closure_relation(γ, ::MeanSpherical, u, r)
    cr = similar(r)
    less_one = r .< 1.0
    larger_one = r .> 1.0

    @. cr[less_one] = -1.0 - γ[less_one]
    @. cr[larger_one] = -u[larger_one]

    return cr
end

function _softmeanspherical!(cr, u1, u2, gamma, idxs)
    @. cr[idxs] = exp(-u1[idxs]) * (1.0 + gamma[idxs] - u2[idxs])
    @. cr[idxs] += -gamma[idxs] - 1.0

    return nothing
end

function closure_relation(γ, ::SoftMeanSpherical, u, r)
    cr = similar(r)
    u₁ = similar(u)
    u₂ = similar(u)

    less_one = r .< 1.0
    larger_one = r .>= 1.0

    u₁[less_one] .= u[less_one]
    u₁[larger_one] .= 0.0
    u₂[less_one] .= 0.0
    u₂[larger_one] .= u[larger_one]

    _softmeanspherical!(cr, u₁, u₂, γ, less_one)
    _softmeanspherical!(cr, u₁, u₂, γ, larger_one)

    return cr
end

function _fmix!(fm, alpha, r, idxs)
    @. fm[idxs] = 1.0 - exp(-alpha * r[idxs])

    return nothing
end

function _hmsa(cr, u1, u2, gamma, fm, idxs)
    @. cr[idxs] = exp(-u1[idxs])
    @. cr[idxs] *= 1.0 + ((exp((gamma[idxs] - u2[idxs]) * fm[idxs]) - 1.0) / fm[idxs])
    @. cr[idxs] += -gamma[idxs] - 1.0

    return nothing
end

function closure_relation(γ, cls::HMSA, u, r)
    cr = similar(r)
    f_mix = similar(r)
    u₁ = similar(u)
    u₂ = similar(u)

    less_one = r .< 1.0
    larger_one = r .>= 1.0

    u₁[less_one] .= u[less_one]
    u₁[larger_one] .= 0.0
    u₂[less_one] .= 0.0
    u₂[larger_one] .= u[larger_one]

    _fmix!(f_mix, cls.α, r, less_one)
    _fmix!(f_mix, cls.α, r, larger_one)

    _hmsa(cr, u₁, u₂, γ, f_mix, less_one)
    _hmsa(cr, u₁, u₂, γ, f_mix, larger_one)

    return cr
end

@inline function closure_relation(γ, ::PercusYevick, u, r)
    return closure_relation(γ, -γ .+ log.(1.0 .+ γ), u, r)
end

@inline function closure_relation(γ, ::ModifiedVerlet, u, r)
    new_br = similar(γ)

    @inbounds for i in eachindex(γ)
        if γ[i] < 0.0
            new_br[i] = -0.5 * γ[i]^2
            new_br[i] /= (1.0 - 0.8 * γ[i])
        else
            new_br[i] = -0.5 * γ[i]^2
            new_br[i] /= (1.0 + 0.8 * γ[i])
        end
    end

    return closure_relation(γ, new_br, u, r)
end
