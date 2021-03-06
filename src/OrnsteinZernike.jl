module OrnsteinZernike

using LinearAlgebra
using FFTW
using Romberg

include("potentials.jl")
export Potential, PseudoHS, HardSphere, SquareWell, LennardJones
include("closures.jl")
export HypernettedChain,
    PercusYevick,
    ModifiedVerlet,
    MeanSpherical,
    SoftMeanSpherical,
    HMSA,
    Parameters,
    Structure,
    Result,
    Interaction
include("types.jl")
include("fft.jl")
include("utils.jl")
include("solver.jl")
export solve

end
