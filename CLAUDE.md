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
   Exception (plans/chain-bmo-fourier.md, D3/D4): BeamletOptics must stay unchanged and
   WaveOpticsPropagation is third party, so their glue lives here as
   `OpticsBaseBeamletOpticsExt` and `OpticsBaseWaveOpticsPropagationExt`; the BMO glue may
   move into BMO once OpticsBase is registered. Further temporary exception
   (plans/field-hierarchy.md): OpticSim.jl (third party) glue lives here as
   `OpticsBaseOpticSimExt` (rays in/out); it moves into the solver package later, too.
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
abstract type AbstractOpticalData end                       # root; interface port, wavelength
abstract type AbstractOpticalField{N} <: AbstractOpticalData end   # Maxwell fields, N = 1 scalar, N = 3 vectorial
  # aliases AbstractScalarField = AbstractOpticalField{1}, AbstractVectorField = AbstractOpticalField{3}
abstract type AbstractRayBundle <: AbstractOpticalData end  # geometrical optics; positions, directions

RayBundle          # <: AbstractRayBundle: per ray position and direction; λ + port per bundle (what ray tracers exchange)
PolarizedRayBundle # <: AbstractRayBundle: adds per ray opl, power, unit transverse 3D polarization vector (input of DebyeWolf)
SampledField{N}    # <: AbstractOpticalField{N}: complex field on plane/volume; N = 1 scalar, N = 3 vectorial; with grid + port
PlaneWaveSpectrum  # <: AbstractOpticalField{N}: samples on the k-sphere: unit directions, spectral density ℰ per solid angle (V/m/sr), solid angles (sr); λ + port (origin = phase reference)
ModalField         # reference to a mode basis + complex coefficients (later)
```

Fields go to field solvers, rays go to other ray tracers or into a ray-to-field converter;
OpticsBase knows no solver, so beamlets, detectors etc. never appear in a format. Rays must
not be handed over at caustics (see "Physical validity at the port" in `conventions.md`).

`PlaneWaveSpectrum` is the central node between the ray and the wave world: a ray is a
k vector with amplitude and phase. Format names describe the representation, not the
dimensionality — not every format is a 2D/3D grid. The name is deliberately not
`AngularSpectrum`: WaveOpticsPropagation exports its propagator under that name, and
exported OpticsBase names must not collide with those of solver packages.

## Ports

Every handover happens at a port: a surface in global coordinates with origin, normal
and local axes. The port also fixes the polarization basis and the medium (refractive
index). No field handover without a port. OpticsBase has no global optical axis. The
local frame `(u, v, n)` is right-handed, `n` points downstream (`n·k > 0`), and `u` is
always passed explicitly. OpticsBase checks that representations fit together, not that
they are physically valid at the port: the user chooses ports where the representation
holds (no rays at caustics, no undersampled phase in sampled fields).

## Conventions

Binding; the canonical text is [docs/src/conventions.md](docs/src/conventions.md) — keep
it and this summary in sync.

- **Time convention:** exp(−iωt). Plane wave: exp(i(k·r − ωt)).
- **Units:** SI internally (m, s, W). Unitful only optionally via extension.
- **Phase reference:** only the carrier exp(−iωt) is excluded; fields store the full
  spatial phase. OPL is absolute, counted from the reference point of the coherent
  component (set by the first solver of the chain); the port origin is r₀ only for
  position-dependent phase terms. Ray phase = k₀·OPL + arg(phasor). If a format deviates,
  document it explicitly.
- **Power normalization:** E is the physical peak amplitude in V/m (real field
  Re(E·e^{−iωt})); P = n/(2Z₀) ∫|E|² dA with n the port index and
  Z₀ = 376.730313668 Ω (same literal as BMO's `Z_vacuum`). Exact along n, paraxial
  otherwise; no obliquity weighting in E. One factor, used everywhere
  (`power_normalization`).
- **Polarization:** 3D E vectors in the global frame. Jones vectors in a right-handed
  basis: `(u, v)` at the port, for oblique rays the ray basis `R(n→d)·(u, v)` (minimal
  rotation, Richards–Wolf). Circular handedness by helicity: right = positive helicity =
  `(1, i)/√2`, rotating u → v under exp(−iωt).
- **Coherence:** each field is one monochromatic, coherent component; polychromatic or
  incoherent light is a collection of fields. Components add coherently (fields) only
  within a coherence group, incoherently (intensities) between groups; coherence groups
  are not yet in the API.
- **Grids:** regular grids in port-local coordinates, sample i at (i − (N÷2 + 1))·Δ, i.e.
  the port origin is the fftshift center; no grid offsets (move the port instead).
- **Vector fields:** transversality k·E = 0 must hold for all `PlaneWaveSpectrum` samples.
- **Plane-wave spectra:** E(r) = Σ w ℰ exp(i k s·(r − r₀)), ℰ per solid angle, w in sr,
  r₀ the port origin; `total_power` is the exact flux κ λₘ² Σ w ‖ℰ‖² (differs from the
  paraxial `SampledField` power by O(θ²)). Converging rays become plane waves (Debye) with
  the factor −i.

When a convention is unclear: do not guess. Ask me and record the decision in
`docs/src/conventions.md`.

## Converters (initial scope)

- No beamlet summation in OpticsBase: BMO sums its own beamlets onto a `SampledField`
  (`SampledField(detector, grid)` in the BMO glue extension); rays are exchanged as
  `RayBundle`, never as beamlets (plans/field-hierarchy.md, D3/D5).
- `SampledField` → `PlaneWaveSpectrum`: FFT (`PlaneWaveDecomposition`, FFTW extension);
  NUFFT later.
- `PolarizedRayBundle` → `PlaneWaveSpectrum`: Debye approximation for converging bundles
  (`DebyeWolf`), solid angles per ray from Voronoi cells (DelaunayTriangulation
  extension) or given explicitly.
- `PlaneWaveSpectrum` → `SampledField` (focus/far field, planes and volumes): direct
  summation (`PlaneWaveSummation`, exact, O(M·N_points)); chirp-z for single planes and
  3D/4D gridding (NUFFT, cf. Lorbeer et al., Opt. Express 23, 3341 (2015)) later as fast
  paths (plans/plane-wave-spectrum.md, D4).
- Later: `SampledField` → `PolarizedRayBundle` (Gaussian beam decomposition or phase gradient).
  Prototype early — this return path shows whether the abstractions hold.

## Shared test suite

A common test suite that every implementation must pass (as a submodule or a separate
package `OpticsBaseTests`):
- Energy conservation across every converter.
- Round trips, e.g. `SampledField` → `PlaneWaveSpectrum` → `SampledField`.
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

Integration tests with the heavy weak dependencies (BeamletOptics needs Julia ≥ 1.12,
WaveOpticsPropagation pulls CUDA.jl and Zygote) live in their own environment
`test/integration/` with its own `TEST_MODULES` list and a separate CI job on Julia 1.
Never add these packages to the main test targets. The light weak dependencies FFTW and
DelaunayTriangulation are in the main test targets; test modules load them themselves
(test missing-extension hints before the `using`).

- Single integration module (or several), from the repo root:
  ```
  julia --project=test/integration -e 'using Pkg; Pkg.instantiate(); include("test/integration/runtests.jl")' TestChain
  ```
- All integration modules: the same command without a module name.

The "full suite" at the end of a plan means both: the core suite and all integration
modules.

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
  the extension `ext/OpticsBaseBeamletOpticsExt.jl` (here, see settled decision 3).
- BeamletOptics itself is not changed for OpticsBase. The extension adapts to BMO's
  conventions (detector frame, amplitude normalization); the mapping is listed in the
  Rationale of `plans/core-types.md` ("BMO compatibility").
