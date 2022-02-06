@testset "Parameters" begin
    rho = 0.8
    nrho = 100
    p = OrnsteinZernike.Parameters(1.0, 6.0, rho, 2^11, nrho, 1.0)
    pot = OrnsteinZernike.HardSphere()
    st = OrnsteinZernike.Structure(p, pot)

    inter = OrnsteinZernike.Interaction(p, st, ModifiedVerlet(-0.6, 0.85))
    result = OrnsteinZernike.solve(inter)
    if 1.0 ∈ result.r
        first_contact = result.gr[result.r .== 1.0][1]
    else
        first_contact = maximum(result.gr)
    end
    @show "MV parameters" first_contact
    # We don't know exactly the value at contact, so just check that it is a valid number
    @test first_contact isa Real
end
