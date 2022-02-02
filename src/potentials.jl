abstract type Potential end

apply_potential(p::Potential, r) = p.f.(r)

struct PseudoHS <: Potential
    f::Function
end

PseudoHS() = PseudoHS(pseudohs)
PseudoHS(λ=50.0) = PseudoHS((x) -> pseudohs(x, λ))

function pseudohs(r, λ)
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

SquareWell() = SquareWell(square_well)
SquareWell(λ=1.5) = SquareWell((x) -> square_well(x, λ))

function square_well(r, λ)
    uij = 0.0

    if r < 1.0
        ur = Inf
    elseif 1.0 <= r < λ
        ur = -1.0
    else
        ur = 0.0
    end

    return uij
end
