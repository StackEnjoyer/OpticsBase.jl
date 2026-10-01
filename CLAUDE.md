# OpticsBase.jl

A minimal interface package for coupling optical solvers, i.e. numerical solutions of
Maxwell's equations. It has one exchange format, `PlaneField`: the tangential E and H
field sampled on a plane in 3D space. Example chain: BeamletOptics.jl → `PlaneField` →
BeamletFibers.jl → `PlaneField` → BeamletOptics.jl.

Current repository: `StackEnjoyer/OpticsBase.jl`. Planned home:
`JuliaPhysics/OpticsBase.jl`; on transfer, update the repo URLs in `docs/make.jl` and
`README.md`. The binding physical conventions are in
[docs/src/conventions.md](docs/src/conventions.md). The design and its rationale are in
`plans/huygens-core.md` (plans are local and not checked in), which replaced an earlier,
much larger design (ports, ray bundles, plane-wave spectra, converters, CommonSolve
interface; see git history before the rewrite).

## Settled architecture decisions

1. **One exchange format.** `PlaneField` with tangential `(Eu, Ev, Hu, Hv)` in the local
   frame of the plane. By the surface equivalence theorem it determines the field on both
   sides, so it serves one-way and bidirectional solvers. Normal components are not
   stored (they follow from the curl equations). No scalar type: scalar solvers use `Eu`.
2. **No solver interface.** No problem/solution types, no traits, no converters, no
   CommonSolve. Solvers couple by function composition: a package offers methods that
   return a `PlaneField` and methods that take one.
3. **No solver code.** Numerics that need an FFT, a mesh or a ray tracer belong in solver
   packages. The package only has `coordinates`, `reference_phase`, `power`, `forward`,
   `backward` and the constructors.
4. **Glue lives in the solver packages** (e.g. `BeamletOpticsOpticsBaseExt` in
   BeamletOptics, a hard dependency in BeamletFibers). OpticsBase knows no solver.
5. **Minimal dependencies:** `LinearAlgebra` and `StaticArrays` only. Adding one needs a
   very good reason.
6. **Beyond one plane by convention, not by types:** closed Huygens surfaces,
   polychromatic and incoherent light are `Vector{PlaneField}`.
7. **One non-minimal field, the reference sphere `R`,** keeps converging or diverging
   beams cheap to sample. Tilt is not factored out; turn the plane instead.

Changing any of these needs a plan in `plans/<topic>.md` and approval first.

## Conventions

Binding; the canonical text is [docs/src/conventions.md](docs/src/conventions.md) — keep
it and this summary in sync.

- **Time convention:** exp(−iωt). Plane wave: exp(i(k·r − ωt)).
- **Units:** SI (m, s, W, V/m, A/m); Z₀ = 376.730313668 Ω (same literal as BMO).
- **Amplitudes:** complex peak values; the real field is Re(E e^{−iωt}).
- **Power:** Poynting flux along n, ½ Re ∫(Eu Hv* − Ev Hu*) dA; exact and signed.
- **Reference direction:** n for R = ∞, else the ray of the reference sphere through the
  sample. The E-only constructor and `forward`/`backward` treat each sample as a local
  plane wave along it: exact there, error of order 1 − cos θ at angle θ to it.
- **Medium:** homogeneous, isotropic, lossless, non-magnetic at the plane; no planes
  through inhomogeneous structures (e.g. fiber cross-sections).
- **Frame:** `axes` = (u, v, n), orthonormal and right-handed; n is the reference
  direction (positive power, `forward`); u fixes the polarization basis; no global
  optical axis.
- **Sampling:** sample i at (i − 1 − N÷2)·Δ from the origin (fftshift center); no grid
  offsets.
- **Phase:** physical field = stored field × `reference_phase`; full spatial phase, OPL
  absolute from the reference point of the chain.
- **Polarization:** Jones vectors in (u, v); right-circular = positive helicity =
  (1, i)/√2.

When a convention is unclear: do not guess. Ask and record the decision in
`docs/src/conventions.md`.

## Running Julia

Always use the newest Julia version installed on the machine, not just whatever a
pinned Manifest specifies as a minimum. Where `juliaup` is installed, `julia` is fine;
otherwise pick the executable explicitly, e.g. on Windows:

```
"/c/Users/uitt_hu/AppData/Local/Programs/Julia-1.13.0/bin/julia.exe" --project=. -e '...'
```

The package itself supports Julia ≥ 1.10 (LTS); do not use language features newer than
that in `src/`.

## Tests

Tests are modular: each `test/Test<Name>.jl` is one module, registered in the
`TEST_MODULES` list in [test/runtests.jl](test/runtests.jl).

```
julia --project=. -e 'using Pkg; Pkg.test()'                               # all
julia --project=. -e 'using Pkg; Pkg.test(test_args=["TestPlaneField"])'   # one module
```

`TestAqua` runs Aqua.jl. Chain tests with real solvers live in the solver packages, not
here.

## Julia code conventions

- **Type-stable, generic code.** Concrete parametric fields; generic in the element type
  (`Float32`, `Float64`) and the array type (GPU arrays, tested with JLArrays): use `similar`,
  broadcasting and views, no scalar indexing in `src/`.
- **`StaticArrays` for 3-vectors and 3×3 matrices.**
- **No type piracy**, no methods on `Base` functions for foreign types.
- **Validation lives in the constructor** that owns the invariant; functions taking a
  `PlaneField` trust it.
- Style baseline is the [SciML Style Guide](https://github.com/SciML/SciMLStyle);
  formatting via JuliaFormatter with `style = "sciml"` (see `.JuliaFormatter.toml`).

## Documentation philosophy & docstring style

- **Docs live off docstrings.** Pages embed docstrings via `@docs` blocks rather than
  re-explaining them; if a page needs to explain something a docstring is missing,
  extend the docstring.
- **Every public function and type has a docstring with a convention note:** units,
  frame (global vs. local) and normalization of every physical quantity.
- **Standard structure:** indented signature-synopsis line, blank line, prose with
  `` [`Type`](@ref) `` cross-references, then an optional `# Fields` or `# Arguments`
  bullet list.
- **Constructor docstrings must be self-sufficient for the REPL**, since solver authors
  implement against `?PlaneField` alone.

## Documentation build

```
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

Pitfalls:

- **Windows-only link bug:** Documenter treats a bare `[text]` followed by a
  parenthesized aside as a malformed link. Escape literal square brackets as `\[m\]`
  (in docstrings `\\[m\\]`).
- **Section titles must not equal a type name** (e.g. no `# PlaneField` heading).
- Every exported name must appear in a `@docs` block.

## Contributing expectations

New or changed functionality comes with tests, docstrings, and documentation where
relevant. Every new convention or deviation from one is recorded in
`docs/src/conventions.md` in the same change.
