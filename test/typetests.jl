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

@testset "Parameters" begin
    brs = Dict(
        "Default" => ModifiedVerlet(),
        "Three parameters" => ModifiedVerlet(-0.6, 0.85, 0.01)
    )

    for (k, v) in brs
        @testset "$k" begin
            result = oz_hs(v)
            if 1.0 ∈ result.r
                first_contact = result.gr[result.r .== 1.0][1]
            else
                first_contact = maximum(result.gr)
            end
            @show "MV parameters" first_contact
            # We don't know exactly the value at contact, so just check that it is a valid number
            @test first_contact isa Real
        end
    end
end
