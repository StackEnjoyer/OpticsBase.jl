<p align="center">
  <img src="docs/src/assets/logo.svg" width="300" alt="OpticsBase.jl logo">
</p>

# OpticsBase.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://StackEnjoyer.github.io/OpticsBase.jl/dev/)
[![CI](https://github.com/StackEnjoyer/OpticsBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/StackEnjoyer/OpticsBase.jl/actions/workflows/CI.yml?query=branch%3Amain)

A minimal interface package for coupling optical solvers in Julia. Every solver of the
ecosystem, whether it traces beamlets, runs a BPM, sums an angular spectrum or solves
Maxwell's equations on a mesh, is a solution of Maxwell's equations. They exchange one
thing: a `PlaneField`, the tangential electric and magnetic field sampled on a plane in
3D space.

```julia
struct PlaneField{T, A}
    E::A          # nx × ny × 2: (Eu, Ev) in V/m
    H::A          # nx × ny × 2: (Hu, Hv) in A/m
    spacing       # (Δu, Δv) in m
    origin        # plane center, global frame
    axes          # columns u, v, n (right-handed, n = direction of travel)
    λ             # vacuum wavelength
    n             # refractive index at the plane
    R             # radius of the reference sphere (Inf = none)
end
```

By the surface equivalence theorem, tangential E and H on a plane determine the field on
both sides of it, so the format covers one-way and bidirectional solvers alike. Five
functions come with it: `coordinates`, `reference_phase`, `power` (Poynting flux),
`forward` and `backward`.

Solvers couple by function composition. A package offers methods that return a
`PlaneField` and methods that take one:

```julia
f1 = PlaneField(detector; size = (256, 256), spacing = (0.25e-6, 0.25e-6))  # BeamletOptics
f2, stats = propagate(f1, fiber, FDBPM(dz = 1e-6); grid)                   # BeamletFibers
src = WavefrontBeamletDecomposition(f2)                                     # BeamletOptics
```

There is no solver interface, no trait system and no converter registry. The glue lives
in the solver packages. The binding conventions (SI units, exp(−iωt), peak amplitudes,
frames, sampling, phase reference) are on the
[Conventions](https://StackEnjoyer.github.io/OpticsBase.jl/dev/conventions/) page.

**Status:** early development, the API is not stable yet.
