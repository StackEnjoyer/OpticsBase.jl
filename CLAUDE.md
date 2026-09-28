# OpticsBase.jl

A lightweight Julia base package that enables an ecosystem of interoperable optical
solvers. Solvers do not know each other, only `OpticsBase`. Optical propagation data is
handed between solvers via defined exchange formats at ports. Example chain:
BeamletOptics.jl → OpticsBase → fiber solver → OpticsBase → Fourier solver.

Current repository: `StackEnjoyer/OpticsBase.jl`. Planned home:
`JuliaPhysics/OpticsBase.jl` (next to BeamletOptics.jl and WaveOpticsPropagation.jl);
on transfer, update the repo URLs in `docs/make.jl` and `README.md`. Binding physical conventions are in
[docs/src/conventions.md](docs/src/conventions.md).

Architecture role models:
- **SciMLBase.jl / CommonSolve.jl**: shared problem/solution types and generic
  `solve`/`init`/`solve!`/`step!`.
- **Tables.jl**: every source can talk to every sink without N² converters.
- **VirtualLab Fusion ("field tracing")**: source of ideas for converters and ports.

## Settled architecture decisions

1. **No solver code in OpticsBase.** The package contains only abstract types, concrete
   exchange formats, traits, ports, conventions and generic converters.
2. **Minimal dependencies.** Allowed hard deps: `StaticArrays`, `CommonSolve`,
   `LinearAlgebra`. Everything else (FFTW, NFFT, CUDA, Unitful) only as a package
   extension (`[weakdeps]` + `ext/`).
3. **Solvers connect via extensions.** Pairwise specialized converters live as extensions
   in the respective solver package (e.g. `BeamletOpticsOpticsBaseExt`), not here.
4. **Converters are explicit objects**, never implicit `Base.convert`. They carry their
   parameters (grid, sampling, number of modes, coherence assumption), because
   conversions are almost always approximations.
5. **Solver interface via CommonSolve:**
   ```julia
   prob = PropagationProblem(field_in, system, port_out)
   sol  = solve(prob, alg)   # sol.field isa AbstractOpticalField
   ```
   No type piracy: the first argument of `solve`/`init`/`solve!`/`step!` is always a
   type owned by OpticsBase or the solver package.
6. **Traits for compatibility checks:** `input_representation(alg)`,
   `output_representation(alg)`, `is_vectorial(field)`, `is_coherent(field)`. Pipelines
   use them to check whether stages fit together and report missing converters instead
   of computing something wrong.

## Exchange formats (initial scope)

```julia
abstract type AbstractOpticalField end

RayBundle        # rays/beamlets: pos, dir, opl, power, polarization, λ, coherence_id
SampledField{N}  # complex field on plane/volume; N = 1 scalar, N = 3 vectorial; with grid + port
AngularSpectrum  # samples on the k-sphere: k vectors, E vectors, weights (Jacobian/apodization), λ
ModalField       # reference to a mode basis + complex coefficients (later)
```

`AngularSpectrum` is the central node between the ray and the wave world: a ray is a
k vector with amplitude and phase. Format names describe the representation, not the
dimensionality — not every format is a 2D/3D grid.

## Ports

Every handover happens at a port: a surface in global coordinates with origin, normal
and local axes. The port also fixes the polarization basis. No field handover without a
port. OpticsBase has no global optical axis.

## Conventions

Binding; the canonical text is [docs/src/conventions.md](docs/src/conventions.md) — keep
it and this summary in sync.

- **Time convention:** exp(−iωt). Plane wave: exp(i(k·r − ωt)).
- **Units:** SI internally (m, s, W). Unitful only optionally via extension.
- **Phase reference:** carrier wave not included; OPL referenced to the port origin. If a
  format deviates, document it explicitly.
- **Power normalization:** ∫|E|² dA = P with one documented factor for refractive index
  and impedance; one factor, used everywhere.
- **Polarization:** Jones vectors in the local port basis; 3D E vectors in the global
  frame; one handedness convention for circular polarization, applied everywhere.
- **Coherence:** polychromatic fields are lists of monochromatic components. Each
  component carries a coherence group; add coherently (fields) only within a group,
  incoherently (intensities) between groups.
- **Vector fields:** transversality k·E = 0 must hold for all `AngularSpectrum` samples.

When a convention is unclear: do not guess. Ask me and record the decision in
`docs/src/conventions.md`.

## Converters (initial scope)

- `RayBundle` → `SampledField`: coherent summation of Gaussian beamlets onto a grid at
  the port.
- `SampledField` ↔ `AngularSpectrum`: FFT or NUFFT (via extension).
- `AngularSpectrum` → `SampledField` (focus/far field): Debye–Wolf via chirp-z for single
  planes; 3D/4D gridding (NUFFT, cf. Lorbeer et al., Opt. Express 23, 3341 (2015)) only
  for volumes/time.
- Later: `SampledField` → `RayBundle` (Gaussian beam decomposition or phase gradient).
  Prototype early — this return path shows whether the abstractions hold.

## Shared test suite

A common test suite that every implementation must pass (as a submodule or a separate
package `OpticsBaseTests`):
- Energy conservation across every converter.
- Round trips, e.g. `SampledField` → `AngularSpectrum` → `SampledField`.
- Analytical references: Gaussian beam (waist, Rayleigh length, Gouy phase).
- Vectorial: transversality; x-polarized pupil at NA 0.9 (elongated focus, E_z lobes);
  radially polarized pupil (strong on-axis E_z).

## Approach

- **Start small:** first `RayBundle`, `SampledField` (scalar + vectorial), ports,
  conventions and one real chain BeamletOptics → Fourier solver. `ModalField` and
  coherence groups only once the conventions have proven themselves.
- Wrap existing packages that work on plain arrays (e.g. WaveOpticsPropagation.jl) with
  thin wrappers instead of re-implementing their functionality.
- Before larger design decisions, present a plan in `plans/<topic>.md` and wait for
  approval.

## Running Julia

Always use the newest Julia version installed on the machine, not just whatever a
pinned Manifest specifies as a minimum. There is no `juliaup` on this machine, so pick
the right `julia.exe` explicitly, e.g.:

```
"/c/Users/uitt_hu/AppData/Local/Programs/Julia-1.13.0/bin/julia.exe" --project=. -e '...'
```

Re-check installed versions each session (`AppData\Local\Programs\Julia-<version>\bin`)
rather than assuming a specific version stays the newest. The package itself supports
Julia ≥ 1.10 (LTS); do not use language features newer than that in `src/`.

## Tests

Tests are modular: each `test/Test<Name>.jl` is one module, registered in the
`TEST_MODULES` list in [test/runtests.jl](test/runtests.jl). New test files must be added
there.

- Single module (or several): run from the repo root
  ```
  julia --project=. -e 'using Pkg; Pkg.test(test_args=["TestAqua"])'
  ```
- Full suite:
  ```
  julia --project=. -e 'using Pkg; Pkg.test()'
  ```

`TestAqua` runs Aqua.jl (ambiguities, piracy, stale deps, compat bounds); it is slow, so
skip it while iterating on unrelated code.

## Julia code conventions

- **Multiple dispatch, not `if` cascades over types.** Behavior that depends on the
  representation belongs in methods on the respective type or trait.
- **Type-stable functions.** Concrete, parametric struct fields; no abstractly typed
  fields in hot data structures. Check with `@code_warntype` / `JET` when in doubt.
- **`StaticArrays` for 3-vectors and small matrices** (`SVector{3,T}`, `SMatrix{3,3,T}`),
  generic in the element type `T <: Real` so `Float32`, `Float64` and dual numbers work.
- **No type piracy**, no methods on `Base` functions for foreign types, and no implicit
  `Base.convert` between exchange formats (see decision 4).
- **Traits over deep type hierarchies.** Capabilities (vectorial, coherent,
  representation) are expressed as trait functions, so solver packages can opt in
  without subtyping constraints.
- **DRY: precondition/invariant-establishing calls live in the function whose contract
  needs them, not at every call site that happens to reach it.** If a function's
  contract requires some setup or validation (e.g. checking transversality or that a
  field and a port are consistent), that call belongs inside the contract-owning
  function itself. A caller one layer up should trust the callee to uphold its own
  contract rather than re-establishing the precondition; doing so is dead work, not
  useful defensive redundancy. Conversely, do not strip such a call from a public
  function that is also called directly, since that breaks its standalone contract.
- Style baseline is the [SciML Style Guide](https://github.com/SciML/SciMLStyle);
  formatting via JuliaFormatter with `style = "sciml"` (see `.JuliaFormatter.toml`).

## Documentation philosophy & docstring style

- **Docs live off docstrings, not hand-duplicated prose.** A `.md` page embeds
  docstrings via `@docs` blocks (and figures where useful) rather than re-explaining in
  free text what a docstring should already cover. If a page needs to explain something
  a docstring is missing, extend the docstring — don't compensate for it in the page.
- **Every public function and type has a docstring with a convention note:** units,
  frame (global vs. local port frame), and normalization of every physical quantity it
  takes or returns.
- **Standard structure:** indented signature-synopsis line, blank line, prose
  description using `` [`Type`](@ref) ``/`` [`func`](@ref) `` cross-references, then an
  optional `# Fields` (structs) or `# Arguments` (functions) bullet list
  (`` - `name`: description ``). Extra `# <Section>` headers (e.g. `# Conventions`) are
  fine if an existing one doesn't already fit; don't repeat prose and bullets.
- **Constructor docstrings must be self-sufficient for the REPL.** `?SampledField` is
  often a user's only source of information — the constructor docstring must cover
  units, frames, sign/phase conventions, physical assumptions and limitations, not just
  parameter types.
- **Type/abstract type docstrings must state the interface they define or fulfill** —
  which functions a subtype needs to implement, and which traits a concrete type
  participates in (use an `# Interface` section). This matters doubly here: solver
  authors implement OpticsBase interfaces from the docstrings alone.
- **Converter docstrings state their approximation:** what is assumed (sampling,
  paraxiality, coherence, number of modes) and what is conserved (e.g. power).

## Documentation build

Documenter.jl, pages in `docs/src/`. The `docs` environment references the local
package via `[sources]`, so from the repo root:

```
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

The result is written to `docs/build/`. Pitfalls:

- **Windows-only link bug:** Documenter treats a bare `[text]` followed by a
  parenthesized aside — e.g. `` `λ` in [m] (vacuum wavelength) `` — as a malformed
  markdown link, which fails the build on Windows with "colons not allowed in paths".
  Escape literal square brackets as `\[m\]` when they aren't meant to be a link. This
  hits unit annotations constantly here.
- **Section titles must not equal a type name** (e.g. no `# SampledField` heading):
  `@ref` then confuses the heading with the embedded docstring.
- Example-generating code in `@example` blocks must end with a suppressing statement
  (`nothing # hide` or a trailing `;`), otherwise the last expression's `show` output
  leaks into the rendered page.
- Every exported name must appear in some `@docs`/`@autodocs` block; the build fails on
  missing docstrings otherwise.

## Contributing expectations

New or changed functionality comes with tests, docstrings, and documentation/examples
where relevant. Every new convention or deviation from one is recorded in
`docs/src/conventions.md` in the same change.

## Relation to BeamletOptics.jl (secondary)

BeamletOptics (BMO, `~/.julia/dev/BeamletOptics`) is the first client and the source of
the first real chain. Keep in mind when designing the interface, but do not adopt
BMO-specific conventions in OpticsBase:

- BMO uses the global +y axis as its optical axis and right-handed frames with CCW
  rotations. OpticsBase has no global optical axis; the mapping onto ports happens in
  `BeamletOpticsOpticsBaseExt` inside BMO.
- BMO-side code (the extension) follows BMO's own `CLAUDE.md`, not this file.
