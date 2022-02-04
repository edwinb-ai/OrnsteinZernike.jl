@testset "Results" begin
    phi = 0.8
    nrho = 5000
    p = OrnsteinZernike.Parameters(1.0, 5.0, rho, 2^12, nrho)
    pot = OrnsteinZernike.HardSphere()
    st = OrnsteinZernike.Structure(p, pot)
    cls = OrnsteinZernike.PercusYevick()

    inter = OrnsteinZernike.Interaction(p, st, cls)

    result = OrnsteinZernike.solve(inter; npoints=npoints)
    first_contact = max(result.gr)
    @show first_contact
    @assert first_contact isa Real
end
