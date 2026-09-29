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

## Solver-side glue that lives here temporarily

`OpticsBase` knows no solver. BeamletOptics.jl is unchanged and OpticSim.jl is third
party, so their glue lives in this repository as the extensions
`OpticsBaseBeamletOpticsExt` and `OpticsBaseOpticSimExt`. It is solver-side code: it may
use the solver's internals, adds no numerics to `src/`, and is meant to move into the
solver packages once `OpticsBase` is registered.

Ray tracers hand rays over as [`RayBundle`](@ref) (geometry only) to other ray tracers,
and fields as [`SampledField`](@ref) to field solvers. Rays are only physical at a port
away from caustics, see [Physical validity at the port](@ref).

## BeamletOptics

Rays and Gaussian beamlets traced with
[BeamletOptics.jl](https://github.com/JuliaPhysics/BeamletOptics.jl) are handed over at a
BeamletOptics `Detector`. After `solve_system!`, the hits of the detector become

  - a [`SampledField`](@ref) on a grid: BeamletOptics sums its beamlets coherently with its
    own field functions. Stigmatic beamlets give a scalar `SampledField{1}`, astigmatic
    beamlets (which carry the polarization) a vectorial `SampledField{3}`;
  - a [`RayBundle`](@ref) (positions and directions of the hits) for other ray tracers;
  - a [`PolarizedRayBundle`](@ref) for polarized rays, e.g. as input to
    [`DebyeWolf`](@ref); BeamletOptics rays carry no power, so the total power is passed.

The reverse direction builds a BeamletOptics `Beam` from a `RayBundle`, one ray per bundle entry.

```julia
using OpticsBase, BeamletOptics

solve_system!(system, beam)
field = SampledField(detector, grid)                       # beamlets → SampledField{1} or {3}
bundle = RayBundle(detector)                               # rays for another ray tracer
pbundle = PolarizedRayBundle(detector; power = P0)         # polarized rays with power
spectrum = convert_field(DebyeWolf(port_focus), pbundle)   # → PlaneWaveSpectrum{3}
```

```@docs
SampledField(::Main.BeamletOptics.Detector, ::AbstractGrid{2})
RayBundle(::Main.BeamletOptics.Detector)
PolarizedRayBundle(::Main.BeamletOptics.Detector)
Main.BeamletOptics.Beam(::RayBundle, ::Integer)
```

## OpticSim

Rays traced with [OpticSim.jl](https://github.com/brianguenter/OpticSim.jl) are converted
to a [`RayBundle`](@ref) at a port and back (OpticSim works in mm, `OpticsBase` in m; the
glue converts explicitly). Together with the BeamletOptics glue this hands a ray fan from
BeamletOptics to OpticSim and back:

```julia
using OpticsBase, BeamletOptics, OpticSim

bundle = RayBundle(detector)   # BeamletOptics rays at plane A → RayBundle
# … OpticSim rays from `bundle`, traced through a lens to plane B, and back into a RayBundle
```

The exact function names are documented in the extension `OpticsBaseOpticSimExt`.

Complete, runnable chains:

  - [`examples/bmo_to_fourier.jl`](https://github.com/StackEnjoyer/OpticsBase.jl/blob/main/examples/bmo_to_fourier.jl):
    BeamletOptics beamlet → detector → `SampledField` → angular-spectrum propagation,
    checked against the analytic Gaussian beam.
  - [`examples/bmo_oap_focus.jl`](https://github.com/StackEnjoyer/OpticsBase.jl/blob/main/examples/bmo_oap_focus.jl):
    polarized BeamletOptics rays over an off-axis parabolic mirror → detector →
    `PolarizedRayBundle` → [`DebyeWolf`](@ref) → [`PlaneWaveSummation`](@ref): the vectorial
    focal field, checked against analytic mirror rays.

## Numerical backends

Two converters use light optional packages for their numerics:

  - [`PlaneWaveDecomposition`](@ref) needs FFTW.jl (`using FFTW`, extension
    `OpticsBaseFFTWExt`).
  - [`DebyeWolf`](@ref) computes solid angles per ray with DelaunayTriangulation.jl
    (`using DelaunayTriangulation`, extension `OpticsBaseDelaunayTriangulationExt`), unless
    they are given explicitly.

Without the package, the conversion throws a `MethodError` whose message names it.
