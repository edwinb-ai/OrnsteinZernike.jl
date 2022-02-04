@testset "Results" begin
    rho = 0.8
    nrho = 50
    p = OrnsteinZernike.Parameters(1.0, 6.0, rho, 2^9, nrho, 1.0)
    pot = OrnsteinZernike.HardSphere()
    st = OrnsteinZernike.Structure(p, pot)
    cls = OrnsteinZernike.PercusYevick()

    inter = OrnsteinZernike.Interaction(p, st, cls)

    result = OrnsteinZernike.solve(inter)
    if 1.0 ∈ result.r
        first_contact = result.gr[result.r .== 1.0][1]
    else
        first_contact = maximum(result.gr)
    end
    ground_truth = 3.47 # For Percus-Yevick, this is the best it can get
    @show first_contact
    @test isapprox(first_contact, ground_truth; atol=0.01)
end
