# Core types: ports, RayBundle, SampledField, problem interface

**Status:** implemented · **Tier:** full

## Goal

Phase 1 of the roadmap in `CLAUDE.md`: ports, `RayBundle`, `SampledField` (scalar +
vectorial), traits, `PropagationProblem`/`solve`, and all open conventions fixed.
Out of scope (next plan): converters, `AngularSpectrum`, BMO extension, Fourier-solver
wrapper, `OpticsBaseTests`, coherence groups, `ModalField`, curved ports, Unitful.

## Change map

| Path | Change |
|------|--------|
| `src/AbstractTypes.jl` | new: abstract types, trait/interface function stubs with docstrings |
| `src/Constants.jl` | new: vacuum impedance, `power_normalization` (κ) |
| `src/Ports.jl` | new: `PlanarPort`, frame transforms, `ray_basis`, Jones helpers |
| `src/Grids.jl` | new: `RegularGrid` |
| `src/SampledField.jl` | new: `SampledField` |
| `src/RayBundle.jl` | new: `RayBundle` |
| `src/Problem.jl` | new: `PropagationProblem`, `PropagationSolution`, `solve`/`init` entry, compatibility check |
| `src/OpticsBase.jl` | edit: includes, exports |
| `test/TestPorts.jl`, `TestGrids.jl`, `TestSampledField.jl`, `TestRayBundle.jl`, `TestProblem.jl` | new |
| `test/runtests.jl` | edit: `TEST_MODULES` |
| `docs/src/conventions.md` | edit: resolve all "Not yet fixed" boxes per D1–D4 |
| `docs/src/ports.md`, `formats.md`, `interface.md`, `docs/make.jl` | new/edit: `@docs` pages |
| `docs/src/reference.md` | edit: drop the catch-all `@autodocs` (duplicates the `@docs` pages and fails the build), keep an `@index` |

## API delta

Everything is new (the package currently only re-exports CommonSolve). Signatures:

```diff
+ abstract type AbstractOpticalField end
+ abstract type AbstractPort end
+ abstract type AbstractGrid{D} end
+ abstract type AbstractPropagationAlgorithm end

+ is_vectorial(field::AbstractOpticalField)::Bool
+ is_coherent(field::AbstractOpticalField)::Bool
+ port(field::AbstractOpticalField)::AbstractPort
+ wavelength(field::AbstractOpticalField)::Real                 # vacuum, m
+ total_power(field::AbstractOpticalField)::Real                # W
+ input_representation(alg::AbstractPropagationAlgorithm)::Type
+ output_representation(alg::AbstractPropagationAlgorithm)::Type

+ const VACUUM_IMPEDANCE                                        # Ω
+ power_normalization(port::AbstractPort)                        # κ, D3

+ PlanarPort(origin, normal, u; refractive_index = 1)
+ origin(port), normal(port), local_axes(port), refractive_index(port)
+ to_local(port, r), to_global(port, ξ)
+ ray_basis(port, dir) -> (x̂, ŷ)
+ jones_to_global(port, dir, jones) -> SVector{3,Complex}
+ circular_jones(helicity::Integer) -> SVector{2,Complex}       # D2

+ RegularGrid(dims::NTuple{D,Integer}, spacing::NTuple{D,Real})
+ size(grid), spacing(grid), coordinates(grid, dim)

+ SampledField(E::AbstractArray{<:Complex}, grid, port, wavelength)
+ grid(field), field_array(field)

+ RayBundle(port, wavelength, position, direction, opl, power, phasor; beamlet = nothing)
+ length(bundle), has_beamlets(bundle)

+ PropagationProblem(field, system, port_out)
+ PropagationSolution(field, prob, alg, stats = nothing)
+ solve(prob::PropagationProblem, alg::AbstractPropagationAlgorithm; kwargs...)
+ init(prob::PropagationProblem, alg::AbstractPropagationAlgorithm; kwargs...)
+ __solve(prob, alg; kwargs...), __init(prob, alg; kwargs...)  # implemented by solvers (D10)
+ is_compatible(field, alg)::Bool, check_compatibility(field, alg)
+ struct MissingConverterError <: Exception
```

Type hierarchy (owner in parentheses):

```text
AbstractOpticalField                         (W1)
├── RayBundle{N,T,P,B}                       (W2)
└── SampledField{N,D,T,A,G,P}                (W3)
AbstractPort                                 (W1)
└── PlanarPort{T}                            (W1)
AbstractGrid{D}                              (W1)
└── RegularGrid{D,T}                         (W3)
AbstractPropagationAlgorithm                 (W1; subtypes in solver packages)
PropagationProblem{F,S,P}                    (W4)
PropagationSolution{F,Pr,A,St}               (W4)
```

| Symbol group | Status | Breaking | Call sites |
|---|---|---|---|
| all of the above | new | no | none (new package) |

Export policy: export types, `solve`/`init`/`solve!`/`step!`, trait functions,
`total_power`, `is_compatible`, `check_compatibility`. Accessors and geometry helpers
(`port`, `wavelength`, `grid`, `origin`, `normal`, `to_local`, `ray_basis`, …) are
documented but not exported, to avoid clashes when users load a solver package (e.g.
BMO) alongside. On Julia ≥ 1.11 they are marked `public` (via a version-guarded
`eval`, so 1.10 still loads).

## Workstreams

| W  | Owns | Needs | Executor | Summary |
|----|------|-------|----------|---------|
| W1 | `src/AbstractTypes.jl`, `src/Constants.jl`, `src/Ports.jl`, `test/TestPorts.jl` | – | lead | abstract types, trait stubs, κ, `PlanarPort`, Jones helpers |
| W2 | `src/RayBundle.jl`, `test/TestRayBundle.jl` | W1 | opus | `RayBundle` with validation |
| W3 | `src/Grids.jl`, `src/SampledField.jl`, `test/TestGrids.jl`, `test/TestSampledField.jl` | W1 | opus | `RegularGrid`, `SampledField`, `total_power` |
| W4 | `src/Problem.jl`, `test/TestProblem.jl` | W1 | opus | problem/solution types, `solve`/`init` entry, compatibility |
| –  | `src/OpticsBase.jl`, `test/runtests.jl`, `docs/**`, `CLAUDE.md` | all | lead | integration: wiring, docs pages, conventions |

Before spawning W2–W4 the lead wires all includes, exports and `TEST_MODULES` entries,
creates empty stub files and commits W1 plus the wiring, so implementers touch only
their own files. W2–W4 run in separate git worktrees: every `Pkg.test` compiles the
whole package, so a half-written file of one workstream would break the test runs of
the others. The file sets are disjoint, so merging back is conflict-free.

## Acceptance

**W1**
- [ ] `TestPorts` passes.
- [ ] `local_axes` is orthonormal and right-handed: `‖AᵀA − I‖ < 10⁻¹²`, `u × v = n`.
- [ ] `to_local(p, to_global(p, ξ)) ≈ ξ` (`rtol 10⁻¹²`); `ray_basis(p, normal(p)) == (u, v)`; `x̂ × ŷ = dir` for oblique `dir`.
- [ ] D12: for an oblique `dir` at azimuth φ = 45°, θ = 64° the basis maps `ê_ρ` to `ê_θ` and leaves `ê_φ` unchanged (Richards–Wolf), `atol 10⁻¹²`.
- [ ] `circular_jones(+1)` realizes the D2 convention, checked via `Re(J·e^{−iωt})` at `ωt = 0, π/2`.
- [ ] Constructor throws `ArgumentError` for zero normal and for `u ∥ normal`.

**W2**
- [ ] `TestRayBundle` passes: scalar (N = 1) and vectorial (N = 3) bundles, with and without beamlets.
- [ ] Each invariant in Details/W2 has one test that triggers its `ArgumentError`.
- [ ] `is_vectorial`, `is_coherent`, `total_power`, `length` covered.

**W3**
- [ ] `TestGrids` and `TestSampledField` pass.
- [ ] `coordinates(grid, d)[n ÷ 2 + 1] == 0` for even and odd `n`.
- [ ] Scalar field from a 2D array shares memory with the input (no copy).
- [ ] `total_power` of a sampled Gaussian (w = 1 mm, P = 1 W, n = 1.5) satisfies `|P − 1 W| < 10⁻⁹ W`.

**W4**
- [ ] `TestProblem` passes (uses a test-local mock field and mock algorithm only).
- [ ] Incompatible input throws `MissingConverterError` **before** `__solve` runs.
- [ ] Solution whose field has the wrong representation or a port `≠ port_out` throws.
- [ ] `init` + `solve!` path and direct `__solve` path both work.

**Integration**
- [ ] Full suite passes, including `TestAqua` (no piracy, no stale deps, no undefined exports).
- [ ] `docs/make.jl` builds without missing-docstring errors; `conventions.md` has no "Not yet fixed" boxes.
- [ ] `Project.toml` `[deps]` unchanged (StaticArrays, CommonSolve, LinearAlgebra only).

## Decisions

Guiding constraint (user): no decision may require a change to BeamletOptics (BMO).
All mapping happens in the future BMO extension; see "BMO compatibility" in Rationale.

**D1 — Port frame:** a) `(u, v, n)` right-handed, `n` points downstream (fields at a port propagate into `n·k > 0`), `u` must be passed explicitly (chosen)

**D2 — Circular polarization handedness:** a) by helicity (IEEE): right-circular = positive helicity = Jones `(1, i)/√2` in `(u, v)`; with exp(−iωt) the field rotates u → v, i.e. counter-clockwise for an observer facing the source (chosen)

**D3 — Power normalization κ:** a) E is the physical complex field amplitude in V/m (peak, not RMS; real field = Re(E·e^{−iωt})), `P = n/(2Z₀) ∫|E|² dA` with n the port medium index; exact for fields along n, paraxial otherwise; no obliquity weighting in E. `Z₀ = 376.730313668 Ω` (CODATA 2018, the same literal as BMO's `Z_vacuum`) (chosen)

**D4 — Phase reference and carrier:** a) carrier = exp(−iωt) only; fields store the full spatial phase; `opl` is absolute from the reference point of the coherent component (set by the first solver of the chain); the port origin is the reference point r₀ for position-dependent phase terms such as exp(i k·(r − r₀)) (chosen)

**D5 — RayBundle: wavelength and coherence granularity:** a) per bundle: one bundle = one monochromatic, mutually coherent component; polychromatic or incoherent sets are collections of bundles; coherence groups deferred (chosen)

**D6 — Ray amplitude representation:** a) `RayBundle{N}` parallel to `SampledField{N}`: per-ray `power` (W) + unit complex `phasor` ∈ ℂᴺ (N = 1: phase factor; N = 3: E direction in the global frame, `e·dir = 0`); the phasor carries all phase not in `opl` (Fresnel/coating, Gouy) (chosen)

**D7 — Gaussian beamlet data:** a) optional per-ray symmetric complex 2×2 matrix Q (Im Q positive definite), field ∝ exp(i k xᵀQx/2) in `ray_basis`, k = 2πn/λ; `nothing` for pure rays (chosen)

**D8 — SampledField storage:** a) one complex array of size `(grid..., N)`, also for N = 1; components in the global frame (chosen)

**D9 — Grid convention:** a) regular Cartesian grid in port-local coordinates `(ξ, η[, ζ])`, sample `i` at `(i − (n÷2 + 1))·Δ` on every axis, i.e. the port origin sits at the fftshift center; no separate offset (move the port instead) (chosen)

**D10 — Solve entry point:** a) SciML pattern: OpticsBase owns `solve`/`init` on `(PropagationProblem, AbstractPropagationAlgorithm)`, validates input and output, and dispatches to `__solve`/`__init`, which solver packages implement (chosen)

**D11 — Representation traits:** a) `input_representation`/`output_representation` return field types (e.g. `SampledField{3}`, `Union{RayBundle, SampledField}`); compatible ⇔ `field isa input_representation(alg)` (chosen)

**D12 — Local ray basis for oblique rays:** b) minimal rotation: `(x̂, ŷ) = R(n→d)·(u, v)`, R the rotation about `n × d`; defines what a Jones vector means for a ray with `dir ≠ n` and the frame of the beamlet Q in D7 (chosen)

**D13 — Output validation on the `init` path:** a) only `solve` validates the solution; `init` returns the solver's integrator unchanged, so `solve!(init(prob, alg))` is not checked; stated in the `init` docstring (chosen)

## Details

Common rules for all workstreams: follow `CLAUDE.md` (docstring style with
convention notes, type stability, `T <: Real` generic, escape `\[m\]` in docstrings);
invariants are checked in inner constructors and throw `ArgumentError` with a message
naming the violated invariant; default tolerance `rtol = √eps(T)`. Run only your own
test modules (`Pkg.test(test_args=[...])`, see `CLAUDE.md`). Do not edit
`src/OpticsBase.jl`, `test/runtests.jl`, `Project.toml` or docs; they are pre-wired:
every name in API delta is already exported or marked `public` there, and the interface
functions `port`, `wavelength`, `total_power`, `is_vectorial`, `is_coherent`,
`input_representation`, `output_representation` already exist (add methods, do not
redeclare them). `runtests.jl` evaluates each test file in its own module with
`using OpticsBase, Test` already loaded; add further `using` lines yourself if needed.
You work in your own git worktree; the stub files you own exist and are empty.

<details>
<summary>W1 — Abstract types, constants, ports (lead)</summary>

Written down for a lead in a fresh session without the planning conversation.

- `AbstractTypes.jl`: the four abstract types plus function stubs (`function is_vectorial end`
  etc.) for the trait/interface functions in API delta. The docstrings define the
  interface: `AbstractOpticalField` has an `# Interface` section listing `port`,
  `wavelength`, `total_power`, `is_vectorial`, `is_coherent`; `AbstractPropagationAlgorithm`
  lists `input_representation`, `output_representation` and points to `__solve`/`__init` (W4).
  No default methods for the representation traits (a missing method is an interface error).
- `Constants.jl`: `VACUUM_IMPEDANCE = 376.730313668` Ω, the literal from BMO's
  `Z_vacuum` (CODATA 2018, μ₀c). `power_normalization(port) = refractive_index(port) / (2 VACUUM_IMPEDANCE)`
  in the port's element type (D3).
- `PlanarPort{T} <: AbstractPort` with fields `origin::SVector{3,T}`,
  `axes::SMatrix{3,3,T,9}` (columns u, v, n) and `refractive_index::T`. The constructor
  normalizes `normal`, orthogonalizes `u` against it (Gram–Schmidt) and normalizes it, and sets
  `v = n × u` (so `u × v = n`). Throw `ArgumentError` if `‖normal‖` is zero or
  `‖u − (u·n)n‖ < √eps(T)·‖u‖`, or if `refractive_index ≤ 0`. Promote inputs to a common `T`.
- `to_local(port, r) = axesᵀ (r − origin)` returns `(ξ, η, ζ)`; `to_global(port, ξ)` accepts
  length 2 (on the plane, ζ = 0) or length 3.
- `ray_basis(port, dir)` per D12 b): `(x̂, ŷ) = (R u, R v)` with
  ```text
  c = n·d,   k = n × d
  R x = c x + k × x + (k·x) k / (1 + c)
  ```
  Exactly `(u, v)` for `d = n` (k = 0); right-handed with `d`. Requires `d·n > 0` (D1),
  so `1 + c > 1` and R is never singular; throw `ArgumentError` otherwise.
- `jones_to_global(port, dir, J) = J₁x̂ + J₂ŷ` with `(x̂, ŷ) = ray_basis(port, dir)`.
- `circular_jones(h)`: `h = +1` gives right-circular / positive helicity, `(1, i)/√2` under D2 a),
  `h = −1` the conjugate; any other `h` throws `ArgumentError`.
- After W1: create empty stubs for the W2–W4 files, wire includes (order: AbstractTypes,
  Constants, Ports, Grids, SampledField, RayBundle, Problem), exports and `TEST_MODULES`,
  then spawn W2–W4.

</details>

<details>
<summary>W2 — RayBundle</summary>

Struct of arrays, `Vector` storage:

```text
RayBundle{N, T<:Real, P<:AbstractPort, B} <: AbstractOpticalField
  # B = Nothing or Vector{SMatrix{2,2,Complex{T},4}}
  port::P
  wavelength::T                          # vacuum, m
  position::Vector{SVector{3,T}}         # global, m, on the port plane
  direction::Vector{SVector{3,T}}        # global, unit
  opl::Vector{T}                         # m, see conventions (D4)
  power::Vector{T}                       # W per ray
  phasor::Vector{SVector{N,Complex{T}}}  # unit norm (D6)
  beamlet::B                             # Q matrices (D7) or nothing
```

Invariants (one `ArgumentError` test each):
1. `N ∈ (1, 3)`; all per-ray vectors (and `beamlet` if present) have equal length.
2. `wavelength > 0`; `power .≥ 0`.
3. `‖direction‖ ≈ 1`; `direction · normal(port) > 0` (D1).
4. Position on the port plane: `|(p − origin)·n| ≤ √eps(T) · max(‖p − origin‖, λ)`.
5. `‖phasor‖ ≈ 1`; for N = 3 transversality `|Σ eᵢ dᵢ| ≤ √eps(T)` (no conjugation).
6. Q symmetric within `√eps(T)·‖Q‖`; Im Q positive definite (trace > 0 and det > 0).

- Outer convenience: accept `phasor::AbstractVector{<:Number}` for N = 1 (wrap into `SVector{1}`); accept `beamlet` as scalar `1/q` values and build `Q = (1/q)·I`.
- Promote all real inputs to a common `T`.
- Implement: `length`, `port`, `wavelength`, `is_vectorial` (`N == 3`), `is_coherent` (always `true` in phase 1, D5), `total_power` (sum of `power`), `has_beamlets`.
- Tests: build vectorial phasors with `OpticsBase.jones_to_global` for oblique rays and check they pass invariant 5.
- Do not normalize inputs silently; do not add per-ray wavelength or coherence ids.

</details>

<details>
<summary>W3 — RegularGrid and SampledField</summary>

`RegularGrid{D,T} <: AbstractGrid{D}` with fields `dims::NTuple{D,Int}`,
`spacing::NTuple{D,T}` (m, > 0). `coordinates(grid, d)` returns a range with sample `i`
at `(i − (n÷2 + 1))·Δ` (D9); D = 2 means `(ξ, η)` along `(u, v)`, D = 3 adds `ζ` along
`n` with the same centering rule. `size(grid) = dims`.

```text
SampledField{N, D, T<:Real, A, G, P} <: AbstractOpticalField
  # A<:AbstractArray{Complex{T}}, G<:AbstractGrid{D}, P<:AbstractPort
  E::A              # size (size(grid)..., N), V/m, global frame (D8)
  grid::G
  port::P
  wavelength::T     # vacuum, m
```

- Constructor infers N: `ndims(E) == D` → N = 1 via `reshape` (no copy); `ndims(E) == D + 1` → N = `size(E, D+1)`. Error unless N ∈ (1, 3) and `size(E)[1:D] == size(grid)`.
- Keep `A` generic (`AbstractArray`) so GPU arrays work later; do not `collect`.
- `field_array(f)` returns `E`; `grid(f)`, `port(f)`, `wavelength(f)`; `is_vectorial` (`N == 3`); `is_coherent` (`true`).
- `total_power(f)` for D = 2: `power_normalization(port) · Σ|E|² · Δξ·Δη` (all components). For D = 3 throw `ArgumentError` (power through a volume is undefined).
- Test the Gaussian power case from Acceptance: `|E₀|² = 4 Z₀ P/(n π w²)`, grid 512 × 512, Δ = 20 µm.
- Do not add FFT code or any dependency.

</details>

<details>
<summary>W4 — Problem, solution, solve entry point</summary>

```text
PropagationProblem{F<:AbstractOpticalField, S, P<:AbstractPort}
    (field, system, port_out)
PropagationSolution{F<:AbstractOpticalField, Pr<:PropagationProblem,
                    A<:AbstractPropagationAlgorithm, St}
    (field, prob, alg, stats)   # stats = nothing by default
```

Flow of `solve(prob, alg; kw...)`:
`check_compatibility(prob.field, alg) → sol = __solve(prob, alg; kw...) → _check_solution(sol, prob, alg) → sol`.

1. `is_compatible(field, alg) = field isa input_representation(alg)` (D11).
2. `check_compatibility` throws `MissingConverterError(typeof(field), input_representation(alg))`; its `showerror` says which converter is missing.
3. `_check_solution`: `sol.field isa output_representation(alg)` and `port(sol.field) == prob.port_out`, otherwise `ArgumentError`.
4. `init(prob, alg; kw...)` checks compatibility, then returns `__init(prob, alg; kw...)`.
5. Fallback `__solve(prob, alg; kw...) = solve!(__init(prob, alg; kw...))`, so a solver implements either `__solve` or `__init` + `CommonSolve.solve!` on its own integrator type. Check the solution in `solve` only (not in the fallback) so it runs exactly once.

- Methods are on `CommonSolve.solve`/`CommonSolve.init` with OpticsBase-owned first argument (no piracy).
- `__solve`/`__init` have no methods except the fallback; missing implementations surface as `MethodError`.
- Tests use a test-local `MockField <: AbstractOpticalField` (implements `port`) and mock algorithms (one via `__solve`, one via `__init` + `solve!`, one returning a wrong port) — do not depend on `RayBundle`/`SampledField`.
- Docstrings of `AbstractPropagationAlgorithm`, `__solve`, `__init` form the solver-author interface: list exactly what a solver package must implement.

</details>

<details>
<summary>Deviations</summary>

- W1: `ray_basis` normalizes `dir` internally and returns exactly `(u, v)` when `dir ∥ n`
  (checked on the unnormalized input). Directions with `d·n > 0` only by rounding
  (grazing) are accepted, consistent with `RayBundle` invariant 3.
- W2: the outer `RayBundle` constructor copies its inputs into new `Vector`s; N is
  inferred from the first phasor unless `SVector`s are passed; `has_beamlets` is a runtime
  check.
- W3: extra invariants (`dims ≥ 1`, finite spacings, `wavelength > 0`), `size(grid, d)`,
  and a type-stable `SampledField{N}(E, grid, port, λ)` constructor; `coordinates` is a
  `StepRangeLen` anchored at 0 in the center sample.
- W4: constructor synopses live in the type docstrings; `_check_solution` also rejects a
  `__solve` that does not return a `PropagationSolution`.
- Process: the implementer worktrees were cut from the scaffold commit instead of the W1
  commit; the implementers tested against an exported copy of the W1 commit and the lead
  copied the files over.
- Integration: `runtests.jl` evaluates each test file in its own module; the module
  docstring moved to the Home page after dropping the catch-all `@autodocs`.

</details>

<details>
<summary>Rationale</summary>

- Struct-of-arrays `RayBundle`: cache/SIMD friendly, cheap column access, maps onto
  Tables.jl later. `Vector` storage for now; GPU-generic storage would add type
  parameters and is not needed for the first chain.
- Invariants in inner constructors follow the DRY rule in `CLAUDE.md`: the type owns its
  contract, consumers never re-check.
- Accessors not exported: BMO and other solvers export similar names (`refractive_index`
  etc.); qualified access avoids ambiguity warnings for users loading both.
- W4 tests use mocks so W2–W4 run in parallel after W1.
- `total_power` for `SampledField` uses κ from W1, so D3 is implemented in one place.

Rejected options (rejected, do not implement):
- D1 b) `u` optional via a fixed rule — the polarization basis depends on `u`, a hidden
  default would be a silent convention.
- D2 b) traditional right = clockwise facing the source — depends on viewing direction.
- D3 b) normalized amplitudes with κ = 1; D3 with Z₀ from CODATA 2022 (376.730313412 Ω,
  relative difference 7·10⁻¹⁰) — would make BMO and OpticsBase powers disagree unless BMO
  changes its constant.
- D4 b) phase relative to the port origin — loses the relative phase of fields reaching
  one port via different paths; D4 c) removed spatial carrier — can be added later as an
  optional field.
- D5 b) per-ray λ and `coherence_id`.
- D6 b) Jones vectors in the port basis — ill-defined for oblique rays; D6 c) one complex
  vector with `|a|² = P`.
- D7 b) stigmatic scalar q only — cannot hold BMO's astigmatic beamlets; D7 c) defer.
- D8 b) array of `SVector{N,Complex}`.
- D9 b) grid with its own offset relative to the port origin.
- D10 b) solvers implement `CommonSolve.solve` directly.
- D11 b) singleton trait types.
- D12 a) projection basis — deviates from the Richards–Wolf mapping by ≈ 21° in
  polarization azimuth at NA 0.9, φ = 45°.
- D13 b) OpticsBase-owned integrator wrapper.

BMO compatibility (checked against BMO `origin/master` 0135efec; every item is handled
in the future BMO extension, BMO itself stays unchanged):
- D1: BMO `Detector` local `(x, z)` is left-handed with respect to propagation and its
  normal points against the beam → port `n = −normal`, `(u, v) = (x, −z)`.
- D3: BMO amplitudes are `√(2 Z₀ I)` in every medium (no √n) → `E = E_BMO/√n`; BMO
  detector maps multiply beamlet fields by `√|cos θ|` → the extension evaluates the
  beamlet fields itself without that factor. Z₀ is identical by construction.
- D4, D5, D6: BMO's OPL is absolute from the source and inherited by child beams, its
  `PolarizedRay.E0` is a global transverse 3-vector with Fresnel phases, detectors add all
  hits coherently — these map 1:1. BMO rays carry no power → the converter takes the
  per-ray power from the source definition.
- D7: stigmatic beamlet → `Q = (1/q)·I` from `gauss_parameters`; astigmatic beamlet →
  `Q = U·H⁻¹` from the geometric parabasal rays, with `k = 2πn/λ` (BMO's own field
  evaluation uses k₀, identical at n = 1). M² ≠ 1 is not representable by one Q.
- D8: BMO vector detector output is `Matrix{Point3{Complex}}` → one `permutedims` copy.
- D9: BMO detector grids are `LinRange(min, max, n)` around the hit centroid → the
  extension places the port origin at sample `n÷2 + 1`.
- D2, D12: BMO defines neither. Note for the extension: BMO's `XZBasis`/`YZBasis` Jones
  bases are left-handed with respect to +y/+x propagation, so `[1, 0, im]` along +y has
  negative helicity (left-circular under D2).

</details>
