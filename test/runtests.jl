using OrnsteinZernike
using Test

function oz_sw(clsr)
    rho = 0.8
    nrho = 30
    ktemp = 1.5
    p = OrnsteinZernike.Parameters(1.0, 6.0, rho, 2^11, nrho, ktemp)
    pot = OrnsteinZernike.SquareWell(1.5)
    st = OrnsteinZernike.Structure(p, pot)

    inter = OrnsteinZernike.Interaction(p, st, clsr)

    return OrnsteinZernike.solve(inter)
end

function oz_hs(br)
    rho = 0.8
    nrho = 30
    p = OrnsteinZernike.Parameters(1.0, 8.0, rho, 2^11, nrho, 1.0)
    pot = OrnsteinZernike.HardSphere()
    st = OrnsteinZernike.Structure(p, pot)
    inter = OrnsteinZernike.Interaction(p, st, br)
    result = OrnsteinZernike.solve(inter)

    return result
end

@testset "OrnsteinZernike.jl" begin
    include("physicstests.jl")
    include("potentials.jl")
    include("closures.jl")
    include("typetests.jl")
end
