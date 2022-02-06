@testset "Lennard-Jones" begin
    # The parameters in this test come from the NIST computed values
    # in https://mmlapps.nist.gov/srs/LJ_PURE/md.htm
    rho = 0.84
    nrho = 200
    p = OrnsteinZernike.Parameters(1.0, 2.5, rho, 2^12, nrho, 0.85)
    pot = OrnsteinZernike.LennardJones()
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

@testset "Square Well" begin
    rho = 0.8
    nrho = 500
    p = OrnsteinZernike.Parameters(1.0, 2.5, rho, 2^13, nrho, 1.0)
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
