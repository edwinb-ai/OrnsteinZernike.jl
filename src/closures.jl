abstract type Closure end

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
@inline closure_relation(γ, ::MeanSpherical, u, r) = -u

function closure_relation(γ, ::SoftMeanSpherical, u, r)
    cr = similar(r)
    u₁ = 0.0
    u₂ = 0.0

    for i in eachindex(r)
        if r[i] < 1.0
            u₁ = u[i]
            u₂ = 0.0
            cr[i] = exp(-u₁) * (1.0 + γ[i] - u₂)
            cr[i] += -γ[i] - 1.0
        else
            u₁ = 0.0
            u₂ = u[i]
            cr[i] = exp(-u₁) * (1.0 + γ[i] - u₂)
            cr[i] += -γ[i] - 1.0
        end
    end

    return cr
end

function closure_relation(γ, cls::HMSA, u, r)
    u₁ = 0.0
    u₂ = 0.0
    f_mix = 0.0
    cr = similar(r)

    for i in eachindex(r)
        if r[i] < 1.0
            f_mix = 1.0 - exp(-cls.α * r[i])
            u₁ = u[i]
            u₂ = 0.0
            cr[i] = exp(-u₁) * (1.0 + ((exp((γ[i] - u₂) * f_mix) - 1.0) / f_mix))
            cr[i] += -γ[i] - 1.0
        else
            f_mix = 1.0 - exp(-cls.α * r[i])
            u₁ = 0.0
            u₂ = u[i]
            cr[i] = exp(-u₁) * (1.0 + ((exp((γ[i] - u₂) * f_mix) - 1.0) / f_mix))
            cr[i] += -γ[i] - 1.0
        end
    end

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
