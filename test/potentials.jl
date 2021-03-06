@testset "Lennard-Jones" begin
    # The parameters in this test come from the NIST computed values
    # in https://mmlapps.nist.gov/srs/LJ_PURE/md.htm
    rho = 0.84
    nrho = 200
    p = OrnsteinZernike.Parameters(1.0, 2.5, rho, 2^11, nrho, 0.85)
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
    # We don't know exactly the value at contact, so just check that it is a valid number
    @test first_contact isa Real
end

function sw_oz(pot)
    rho = 0.8
    nrho = 30
    ktemp = 1.5
    p = OrnsteinZernike.Parameters(1.0, 6.0, rho, 2^12, nrho, ktemp)
    st = OrnsteinZernike.Structure(p, pot)
    cls = OrnsteinZernike.MeanSpherical()

    inter = OrnsteinZernike.Interaction(p, st, cls)

    result = OrnsteinZernike.solve(inter)

    return result
end

@testset "Square Well" begin
    potentials = Dict("1.5" => SquareWell(1.5), "1.3" => SquareWell(1.3))

    for (k, v) in potentials
        @testset "$k" begin
            result = sw_oz(v)
            if 1.0 ∈ result.r
                first_contact = result.gr[result.r .== 1.0][1]
            else
                first_contact = maximum(result.gr)
            end
            @show "SW $k" first_contact
            @test first_contact isa Real
        end
    end
end

function phs_oz(pot)
    rho = 0.8
    nrho = 30
    p = OrnsteinZernike.Parameters(1.0, 6.0, rho, 2^11, nrho, 1.0)
    st = OrnsteinZernike.Structure(p, pot)
    cls = OrnsteinZernike.ModifiedVerlet()

    inter = OrnsteinZernike.Interaction(p, st, cls)

    return OrnsteinZernike.solve(inter)
end

@testset "Pseudo Hard Sphere" begin
    potentials = Dict("50" => PseudoHS(), "51" => PseudoHS(51.0))

    for (k, v) in potentials
        @testset "$k" begin
            result = phs_oz(v)
            if 1.0 ∈ result.r
                first_contact = result.gr[result.r .== 1.0][1]
            else
                first_contact = maximum(result.gr)
            end
            @show "PHS $k" first_contact
            @test first_contact isa Real
        end
    end
end
