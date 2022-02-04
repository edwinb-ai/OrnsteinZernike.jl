@testset "Results" begin
    rho = 0.8
    nrho = 500
    p = OrnsteinZernike.Parameters(1.0, 6.0, rho, 2^9, nrho, 1.0)
    pot = OrnsteinZernike.HardSphere()
    st = OrnsteinZernike.Structure(p, pot)
    cls = OrnsteinZernike.PercusYevick()

    inter = OrnsteinZernike.Interaction(p, st, cls)

    result = OrnsteinZernike.solve(inter)
    first_contact = maximum(result.gr)
    ground_truth = 3.47 # For Percus-Yevick, this is the best it can get
    @show first_contact
    @test isapprox(first_contact, ground_truth; atol=0.01)
end
