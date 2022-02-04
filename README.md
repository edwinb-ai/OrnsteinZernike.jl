# OrnsteinZernike.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://edwinb-ai.github.io/OrnsteinZernike.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://edwinb-ai.github.io/OrnsteinZernike.jl/dev)
[![Build Status](https://github.com/edwinb-ai/OrnsteinZernike.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/edwinb-ai/OrnsteinZernike.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/edwinb-ai/OrnsteinZernike.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/edwinb-ai/OrnsteinZernike.jl)

This is an attempt to build a general purpose code for solving the 
[Ornstein-Zernike equation](http://www.sklogwiki.org/SklogWiki/index.php/Ornstein-Zernike_relation),
which is an exact formalism to describe the microstructure of simple and complex
fluids.

## Features
This code provides several features:
- It is _fast_. It can solve the Percus-Yevick equation for the hard sphere fluid in milliseconds (benchmarks pending).
- It is _composable_. Through public APIs, new closures and potentials can be added, and the same code can be reused.
- More to come...

## Installation

As with any `Julia` package, you add this to your environment
```julia
julia> using Pkg
julia> pkg"add https://github.com/edwinb-ai/OrnsteinZernike.jl.git
```

## How to use it

For now, you can look at the tests inside the `physicstests.jl` file within the `tests` directory from this repository.
Documentation is still pending.