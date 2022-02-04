module OrnsteinZernike

using LinearAlgebra
using FFTW

include("potentials.jl")
export Potential, PseudoHS, HardSphere, SmoothSW
include("types.jl")
export HypernettedChain,
    PercusYevick,
    ModifiedVerlet,
    MeanSpherical,
    SoftMeanSpherical,
    HMSA,
    ConstantClosure,
    Parameters,
    Structure,
    Result,
    Interaction

include("fft.jl")
include("utils.jl")
include("closures.jl")
include("solver.jl")
export solve

end
