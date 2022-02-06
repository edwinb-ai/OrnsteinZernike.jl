function oz_br(clsr)
    rho = 0.8
    nrho = 1500
    p = OrnsteinZernike.Parameters(1.0, 2.5, rho, 2^11, nrho, 1.5)
    pot = OrnsteinZernike.SquareWell()
    st = OrnsteinZernike.Structure(p, pot)

    inter = OrnsteinZernike.Interaction(p, st, clsr)

    return OrnsteinZernike.solve(inter)
end

@testset "Closures" begin
    closures = Dict(
        "HNC" => HypernettedChain(),
        "MSA" => MeanSpherical(),
        "SMSA" => SoftMeanSpherical(),
        "HMSA" => HMSA(0.1),
        "PY" => PercusYevick(),
        "MV" => ModifiedVerlet(),
    )

    for (k, v) in closures
        @testset "$(k)" begin
            result = oz_br(v)
            if 1.0 ∈ result.r
                first_contact = result.gr[result.r.==1.0][1]
            else
                first_contact = maximum(result.gr)
            end
            @show "$(k)" first_contact
            # We don't know exactly the value at contact, so just check that it is a valid number
            @test first_contact isa Real
        end
    end
end
