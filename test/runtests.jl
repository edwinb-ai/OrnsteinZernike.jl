using OrnsteinZernike
using Test

@testset "OrnsteinZernike.jl" begin
    include("physicstests.jl")
    include("potentials.jl")
    include("closures.jl")
end
