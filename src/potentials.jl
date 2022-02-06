abstract type Potential end

apply_potential(p::Potential, r) = p.f.(r)

struct PseudoHS <: Potential
    f::Function
end

PseudoHS() = PseudoHS(pseudohs)
PseudoHS(λ::Real) = PseudoHS((x) -> pseudohs(x; λ=λ))

function pseudohs(r; λ=50.0)
    # FIXME: This formula only applies for a 3D system
    dT = 1.4543 + 1.199 / (λ^1.0545)
    A = λ * (λ / (λ - 1.0))^(λ - 1.0)
    B = λ / (λ - 1.0)
    eps_energy = 1.0 / dT

    ur = 0.0

    if r < B
        xr = 1.0 / (r^λ)
        xa = 1.0 / (r^(λ - 1.0))
        ur = A * eps_energy * (xr - xa) + eps_energy
    else
        ur = 0.0
    end

    return ur
end

struct HardSphere <: Potential
    f::Function

    HardSphere() = new(hard_sphere)
end

function hard_sphere(r)
    ur = 0.0

    if r < 1.0
        ur = Inf
    else
        ur = 0.0
    end

    return ur
end

struct SquareWell <: Potential
    f::Function
end

SquareWell(λ::Real) = SquareWell((x) -> square_well(x; λ=λ))

function square_well(r; λ=1.5)
    uij = 0.0

    if r < 1.0
        uij = Inf
    elseif 1.0 <= r < λ
        uij = -1.0
    else
        uij = 0.0
    end

    return uij
end

struct LennardJones <: Potential
    f::Function

    LennardJones() = new(lj)
end

function lj(r)
    r6 = 1.0 / r^6
    r12 = r6 * r6

    return 4.0 * (r12 - r6)
end
