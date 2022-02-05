@testset "Square Well" begin
    rho = 0.8
    nrho = 50
    p = OrnsteinZernike.Parameters(1.0, 8.0, rho, 2^8, nrho, 1.0)
    pot = OrnsteinZernike.SquareWell()
    st = OrnsteinZernike.Structure(p, pot)
    cls = OrnsteinZernike.HypernettedChain()

    inter = OrnsteinZernike.Interaction(p, st, cls)

    result = OrnsteinZernike.solve(inter)
    if 1.0 ∈ result.r
        first_contact = result.gr[result.r .== 1.0][1]
    else
        first_contact = maximum(result.gr)
    end
    @show first_contact
    @test first_contact isa Real
end
