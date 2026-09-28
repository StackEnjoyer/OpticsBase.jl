# OpticsBase.jl

`OpticsBase.jl` is a lightweight base package that enables an ecosystem of interoperable
optical solvers. Solvers do not know each other, only `OpticsBase`: optical propagation
data is passed from one solver to the next through well-defined exchange formats at
ports.

An example chain: a ray/beamlet tracer (e.g.
[BeamletOptics.jl](https://github.com/JuliaPhysics/BeamletOptics.jl)) hands its result
to a fiber solver, which in turn hands its output to a Fourier-optics solver, with
`OpticsBase` in between each stage.

## Design principles

- **No solver code.** `OpticsBase` contains abstract types, exchange formats, traits,
  ports, conventions and generic converters only.
- **Minimal dependencies.** Only `StaticArrays`, `CommonSolve` and `LinearAlgebra` are
  hard dependencies. Everything else (FFT, NUFFT, GPU, units) is loaded via package
  extensions.
- **Solvers connect via extensions.** Pairwise specialized converters live as package
  extensions in the respective solver package, not in `OpticsBase`.
- **Explicit converters.** Conversions between representations are almost always
  approximations. Converters are therefore explicit objects that carry their parameters
  (grid, sampling, number of modes, coherence assumption), never implicit
  `Base.convert` methods.
- **CommonSolve interface.** Solvers implement `solve`/`init`/`solve!`/`step!` from
  [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl) on `OpticsBase` problem types.
- **Traits for compatibility checks.** Pipelines use traits to check whether stages fit
  together and report a missing converter instead of computing something wrong.

The binding physical conventions are listed on the [Conventions](@ref) page.

## Package overview

```@docs
OpticsBase
```
