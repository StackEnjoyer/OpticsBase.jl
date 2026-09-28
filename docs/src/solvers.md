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

A complete, runnable chain (BeamletOptics beamlet → detector → `RayBundle` → beamlet
summation → angular-spectrum propagation, checked against the analytic Gaussian beam) is
in [`examples/bmo_to_fourier.jl`](https://github.com/StackEnjoyer/OpticsBase.jl/blob/main/examples/bmo_to_fourier.jl).
