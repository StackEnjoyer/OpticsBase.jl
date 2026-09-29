<p align="center">
  <img src="docs/src/assets/logo.svg" width="300" alt="OpticsBase.jl logo">
</p>

# OpticsBase.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://StackEnjoyer.github.io/OpticsBase.jl/dev/)
[![CI](https://github.com/StackEnjoyer/OpticsBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/StackEnjoyer/OpticsBase.jl/actions/workflows/CI.yml?query=branch%3Amain)

A lightweight base package for an ecosystem of interoperable optical solvers in Julia.
Solvers do not know each other, only `OpticsBase`: optical propagation data is handed
from one solver to the next through well-defined exchange formats (ray bundles, sampled
fields, plane-wave spectra) at ports with explicit frames and conventions.

`OpticsBase` contains no solver code. It provides abstract types, exchange formats,
traits, ports, conventions and generic converters, and uses the
[CommonSolve.jl](https://github.com/SciML/CommonSolve.jl) interface for solvers.

**Status:** early development, the API is not stable yet.

## Exchange formats

Everything a solver takes as input or returns as output is an `AbstractOpticalData`: one
monochromatic, coherent component given at a port. The hierarchy is ordered by physical
model, Maxwell fields on one side and rays (geometrical optics) on the other:

```text
AbstractOpticalData                  port, wavelength
├── AbstractOpticalField{N}          Maxwell field; N = 3 vectorial, N = 1 scalar
│   │                                (aliases AbstractVectorField, AbstractScalarField)
│   ├── SampledField{N}              complex E on a grid at the port (plane or volume)
│   └── PlaneWaveSpectrum{N}         plane waves: directions, spectral density, solid angles
└── AbstractRayBundle                rays: positions, directions
    ├── RayBundle                    geometry only; what ray tracers hand to ray tracers
    └── PolarizedRayBundle           + optical path length, power, polarization per ray
```

Fields go to field solvers, rays to other ray tracers or into a converter to the wave
world. A solver declares what it needs with `input_representation`, at any level of the
tree: `SampledField{3}` for exactly that encoding, `AbstractVectorField` for any vectorial
field, `AbstractRayBundle` for any rays. Conversions are explicit converter objects:
`DebyeWolf` (`PolarizedRayBundle` → `PlaneWaveSpectrum`), `PlaneWaveSummation`
(`PlaneWaveSpectrum` → `SampledField`) and `PlaneWaveDecomposition` (`SampledField` →
`PlaneWaveSpectrum`).

## Connecting a solver package

A solver package connects itself to `OpticsBase`; `OpticsBase` does not know it. Ideally
the glue lives in a package extension of the solver (`MySolverOpticsBaseExt`), so the
solver does not depend on `OpticsBase` for users who don't chain it.

What a package implements depends on its role:

- **Source or sink only** (e.g. a ray tracer that hands its result on): build or read an
  exchange format at a port, e.g. a `PlanarPort`. Ray tracers hand rays
  (`AbstractRayBundle`, e.g. `RayBundle`) to other ray tracers and fields
  (`AbstractOpticalField`, e.g. `SampledField` or `PlaneWaveSpectrum`) to field solvers.
  Nothing else is needed.
- **Propagation stage** (field in at one port, field out at another), following the
  [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl) pattern
  `init(prob, alg) -> integrator`, `solve!(integrator) -> solution`:
  - an algorithm type `MyAlg <: AbstractPropagationAlgorithm`,
  - the traits `input_representation(::MyAlg)` and `output_representation(::MyAlg)`
    (the field types it accepts and returns),
  - `OpticsBase.__init(prob::PropagationProblem, alg::MyAlg; kwargs...)` returning the
    solver's own integrator type, and `CommonSolve.solve!(integrator)` returning a
    `PropagationSolution` whose field is given at `prob.port_out`. Optionally
    `CommonSolve.step!(integrator)`.
  - Shortcut for solvers without state worth keeping: only
    `OpticsBase.__solve(prob, alg::MyAlg; kwargs...)` returning the `PropagationSolution`
    (then `init` is not supported).

  `prob.system` is the solver's own description of the optical system.
- **Converter** between representations: a type `<: AbstractFieldConverter` with the same
  two traits and `OpticsBase.__convert_field(conv, field)`.
- **New exchange format** (rarely needed): a field type `<: AbstractOpticalField{N}` with
  `port`, `wavelength`, `total_power` and `is_coherent` (`is_vectorial` follows from `N`),
  or a ray type `<: AbstractRayBundle` with `port`, `wavelength`, `positions`,
  `directions` and `length`.

```julia
struct MyAlg <: AbstractPropagationAlgorithm end

OpticsBase.input_representation(::MyAlg) = SampledField
OpticsBase.output_representation(::MyAlg) = SampledField

struct MyIntegrator{P, A, C}
    prob::P
    alg::A
    cache::C   # precomputed data, reusable across solve! calls
end

function OpticsBase.__init(prob::PropagationProblem, alg::MyAlg; kwargs...)
    return MyIntegrator(prob, alg, ...)   # set up from prob.field, prob.system, prob.port_out
end

function CommonSolve.solve!(integ::MyIntegrator)
    field_out = ...   # propagate integ.prob.field to integ.prob.port_out
    return PropagationSolution(field_out, integ.prob, integ.alg)
end

prob = PropagationProblem(field_in, system, port_out)
sol = solve(prob, MyAlg())          # = solve!(init(prob, MyAlg())), plus output checks
```

Why `__init` instead of `init`: `OpticsBase` owns `init` and `solve` for
`PropagationProblem` (as SciMLBase.jl does for its problems). They first check the input
against `input_representation` (a mismatch raises a `MissingConverterError` before any
solver code runs, instead of a wrong result), then call `__init`/`__solve`; `solve` also
checks the returned field against `output_representation` and `prob.port_out`. Solver
packages never add methods to `solve` or `init`. `convert_field` works the same way on top
of `__convert_field`. All inputs and outputs follow the binding
[conventions](https://StackEnjoyer.github.io/OpticsBase.jl/dev/conventions/) (SI units,
exp(−iωt), power normalization, polarization basis at the port).

The glue for BeamletOptics.jl, OpticSim.jl and WaveOpticsPropagation.jl currently lives in this
repository (`ext/`) as a temporary exception, until it can move into those packages.
