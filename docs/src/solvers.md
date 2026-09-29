```@meta
CurrentModule = OpticsBase
```

# Solvers and extensions

`OpticsBase` contains no solver numerics. Solver packages connect through package
extensions that are loaded automatically when both packages are loaded.

## Free-space propagation

Coaxial propagation of a planar [`SampledField`](@ref) with the angular-spectrum method of
[WaveOpticsPropagation.jl](https://github.com/JuliaPhysics/WaveOpticsPropagation.jl):

```julia
using OpticsBase, WaveOpticsPropagation

port_out = PlanarPort(origin + z * n, n, u)          # parallel to the input port
prob = PropagationProblem(field, FreeSpace(), port_out)
sol = solve(prob, AngularSpectrumMethod())            # sol.field isa SampledField at port_out
```

```@docs
FreeSpace
AngularSpectrumMethod
```

## BeamletOptics

Rays and Gaussian beamlets traced with
[BeamletOptics.jl](https://github.com/JuliaPhysics/BeamletOptics.jl) are handed over at a
BeamletOptics `Detector`: after `solve_system!`, the hits of the detector become a
[`RayBundle`](@ref) at a port in the detector plane.

```julia
using OpticsBase, BeamletOptics

solve_system!(system, beam)
bundle = RayBundle(detector)                                        # beamlets → N = 1 or 3 with Q
field = convert_field(GaussianBeamletSummation(grid), bundle)       # → SampledField
```

```@docs
RayBundle(::Main.BeamletOptics.Detector)
```

Two complete, runnable chains:

  - [`examples/bmo_to_fourier.jl`](https://github.com/StackEnjoyer/OpticsBase.jl/blob/main/examples/bmo_to_fourier.jl):
    BeamletOptics beamlet → detector → `RayBundle` → beamlet summation → angular-spectrum
    propagation, checked against the analytic Gaussian beam.
  - [`examples/bmo_oap_focus.jl`](https://github.com/StackEnjoyer/OpticsBase.jl/blob/main/examples/bmo_oap_focus.jl):
    polarized BeamletOptics rays over an off-axis parabolic mirror → detector →
    `RayBundle` → [`DebyeWolf`](@ref) → [`PlaneWaveSummation`](@ref): the vectorial focal
    field, checked against analytic mirror rays.

## Numerical backends

Two converters use light optional packages for their numerics:

  - [`PlaneWaveDecomposition`](@ref) needs FFTW.jl (`using FFTW`, extension
    `OpticsBaseFFTWExt`).
  - [`DebyeWolf`](@ref) computes solid angles per ray with DelaunayTriangulation.jl
    (`using DelaunayTriangulation`, extension `OpticsBaseDelaunayTriangulationExt`), unless
    they are given explicitly.

Without the package, the conversion throws a `MethodError` whose message names it.
