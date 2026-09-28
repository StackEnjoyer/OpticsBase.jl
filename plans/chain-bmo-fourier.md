# First chain: BeamletOptics → OpticsBase → Fourier solver

**Status:** implemented · **Tier:** full

## Goal

One real chain: BMO beamlets at a detector → `RayBundle` → coherent beamlet summation →
`SampledField` → angular-spectrum propagation (WaveOpticsPropagation.jl) → `SampledField`,
validated against the analytic Gaussian beam.

## Change map

| Path | Change |
|------|--------|
| `src/Converters.jl` | new: `AbstractFieldConverter`, `convert_field`, compatibility methods |
| `src/BeamletSummation.jl` | new: `GaussianBeamletSummation` (`RayBundle` → `SampledField`) |
| `src/FreeSpace.jl` | new: `FreeSpace`, `AngularSpectrumMethod`, geometry checks, backend stub |
| `ext/OpticsBaseWaveOpticsPropagationExt.jl` | new: backend call into WaveOpticsPropagation |
| `ext/OpticsBaseBeamletOpticsExt.jl` | new: `RayBundle(::BeamletOptics.Detector)` |
| `src/OpticsBase.jl` | edit: includes, exports, `__init__` (error hint) |
| `Project.toml` | edit: `[weakdeps]`, `[extensions]`, `[compat]` for BMO and WOP; `[deps]` unchanged |
| `test/TestConverters.jl`, `TestBeamletSummation.jl`, `TestFreeSpace.jl` | new (core suite) |
| `test/integration/` | new environment: `Project.toml`, `runtests.jl`, `TestWaveOpticsPropagation.jl`, `TestBeamletOptics.jl`, `TestChain.jl` |
| `test/runtests.jl` | edit: `TEST_MODULES` |
| `.github/workflows/CI.yml` | edit: integration job on Julia 1 |
| `docs/src/converters.md`, `docs/src/solvers.md`, `docs/make.jl`, `docs/Project.toml` | new/edit |
| `CLAUDE.md` | edit: integration test commands; settled decision 3 amended (D3) |

## API delta

```diff
+ abstract type AbstractFieldConverter end
+ convert_field(conv::AbstractFieldConverter, field::AbstractOpticalField) -> AbstractOpticalField
+ input_representation(conv::AbstractFieldConverter), output_representation(conv::AbstractFieldConverter)
+ is_compatible(field::AbstractOpticalField, conv::AbstractFieldConverter)::Bool
+ check_compatibility(field::AbstractOpticalField, conv::AbstractFieldConverter)

+ GaussianBeamletSummation(grid::AbstractGrid{2})  <: AbstractFieldConverter

+ struct FreeSpace end
+ AngularSpectrumMethod(; padding = true, pad_factor = 2, bandlimit = true)  <: AbstractPropagationAlgorithm
+ __solve(prob::PropagationProblem{<:SampledField, FreeSpace}, alg::AngularSpectrumMethod)

+ RayBundle(detector::BeamletOptics.Detector; port = nothing, power = nothing)   # extension
```

| Symbol group | Status | Breaking | Call sites |
|---|---|---|---|
| converter API, `GaussianBeamletSummation` | new | no | none |
| `FreeSpace`, `AngularSpectrumMethod` | new | no | none |
| `is_compatible`/`check_compatibility` for converters | new methods | no | existing algorithm methods unchanged |
| `RayBundle(::Detector)` | new method (extension) | no | none |

Existing signatures (`is_compatible(field, alg::AbstractPropagationAlgorithm)`,
`input_representation(alg)`, `RayBundle(port, λ, …)`) are unchanged; their docstrings are
widened to mention converters.

## Logic

Coherent summation of Gaussian beamlets on the port plane (W1). Per ray j: position p,
unit direction d, ray basis `(x̂, ŷ) = ray_basis(port, d)`, beamlet matrix Q (at p),
k = 2πn/λ, k₀ = 2π/λ, κ = `power_normalization(port)`.

```text
a_j   = √(P_j · k · √det(Im Q) / (κ π)) · exp(i k₀ opl_j) · phasor_j      (peak amplitude)
E(r)  = Σ_j a_j · det(I + sQ)^(−1/2) · exp(i k s) · exp(i k xᵀ Q(s) x / 2)
s     = (r − p)·d,   x = ((r − p)·x̂, (r − p)·ŷ),   Q(s) = Q (I + sQ)⁻¹
det(I + sQ)^(1/2) = √(1 + sμ₁) · √(1 + sμ₂)   (μ eigenvalues of Q, principal roots)
```

1. check `field isa RayBundle` with beamlets (D2), grid 2D
2. allocate `E` of size `(size(grid)..., N)` with the bundle's element type
3. for each ray j: precompute `a_j`, `(x̂, ŷ)`, μ₁, μ₂
4. for each grid point r = `to_global(port, (ξ, η))`, for each ray j:
5.   s, x as above; D = I + sQ; Q(s) = Q·D⁻¹
6.   add `a_j · exp(i k s) · exp(i k xᵀQ(s)x/2) / (√(1 + sμ₁)√(1 + sμ₂))`
7. return `SampledField(E, grid, port, λ)`

Trace (awkward case: oblique stigmatic beamlet, point off its transverse plane): port
normal z, u = x; ray at origin, d = (sin 30°, 0, cos 30°), w = 1 mm, λ = 1 µm, n = 1,
Q = i b I with b = λ/(π w²) = 0.3183 m⁻¹, point r = (w, 0, 0).

| step | value |
|---|---|
| s = r·d | 0.5 mm |
| x̂ = R u, x = (r·x̂, r·ŷ) | x̂ = (cos 30°, 0, −sin 30°), x = (0.866 mm, 0) |
| √(1 + sμ₁)√(1 + sμ₂) = 1 + i b s | 1 + 1.59·10⁻⁴ i |
| exp(i k xᵀQ(s)x/2) | magnitude ≈ exp(−0.75) = 0.472 |
| exp(i k s) | phase 3141.6 rad |
| κ ∫\|E\|² dA over the port plane | P / cos 30° (no obliquity weighting, D3 of core-types) |

## Workstreams

| W  | Owns | Needs | Executor | Summary |
|----|------|-------|----------|---------|
| W1 | `src/Converters.jl`, `src/BeamletSummation.jl`, `test/TestConverters.jl`, `test/TestBeamletSummation.jl` | – | lead | converter API, beamlet summation |
| W2 | `src/FreeSpace.jl`, `ext/OpticsBaseWaveOpticsPropagationExt.jl`, `test/TestFreeSpace.jl`, `test/integration/TestWaveOpticsPropagation.jl` | – | opus | angular-spectrum wrapper |
| W3 | `ext/OpticsBaseBeamletOpticsExt.jl`, `test/integration/TestBeamletOptics.jl` | – | opus | BMO detector → `RayBundle` |
| –  | `src/OpticsBase.jl`, `Project.toml`, `test/runtests.jl`, `test/integration/{Project.toml, runtests.jl, TestChain.jl}`, CI, docs, `CLAUDE.md` | all | lead | wiring, chain test, docs |

Before spawning W2/W3 the lead wires includes, exports, `__init__`, weakdeps, the
integration environment and empty stubs, runs the full suite, commits and pushes, so the
implementer worktrees contain the wiring.

## Acceptance

**W1**
- [ ] `TestConverters`, `TestBeamletSummation` pass.
- [ ] Single beamlet at normal incidence (w = 1 mm, grid ±8w, Δ = w/16): `|P − P₀|/P₀ < 10⁻⁹`, profile equals the analytic Gaussian `max|ΔE|/max|E| < 10⁻¹²`.
- [ ] Oblique beamlet (θ = 20°): `κ∫|E|²dA = P₀/cos θ` within `rtol 10⁻⁴`; matches an independent scalar tilted-Gaussian formula `rtol 10⁻¹⁰`.
- [ ] Two beamlets with phasors 1 and −1 at the same position: `max|E| < 10⁻¹²·max|E₁|`; two beamlets at ±θ: fringe period λ/(2 sin θ) within `rtol 10⁻⁶`.
- [ ] Vectorial bundle yields `SampledField{3}`; pure rays and non-`RayBundle` input throw `MissingConverterError`.

**W2**
- [ ] `TestFreeSpace` passes without WOP loaded: every invalid geometry throws before the backend runs; the missing-backend `MethodError` names WaveOpticsPropagation.
- [ ] Integration: analytic Gaussian (w₀ = 0.25 mm, λ = 1.064 µm) propagated by ±z_R/2 matches `w(z)` within `rtol 10⁻³`, on-axis phase `kz − atan(z/z_R)` within 10⁻³ rad, power within `rtol 10⁻⁶`; N = 3 components match the scalar result.

**W3**
- [ ] Integration: astigmatic BMO beamlet (w₀ = 0.5 mm, λ = 1.064 µm, P₀ = 1 mW) at a detector 0.2 m after the waist: one ray, N = 3, `total_power` = BMO `optical_power` within `rtol 10⁻⁹`, `Q = (1/q)·I` with `q = z − i z_R` within `rtol 10⁻⁶`, Gouy phase `−atan(z/z_R)` within 10⁻⁶ rad.
- [ ] Astigmatic waists (w₀ₓ ≠ w₀ᵧ) give the two q values on the diagonal of Q in the support-axis basis; stigmatic `GaussianBeamlet` gives N = 1 with the same Q; polarized and plain rays give pure rays with the given total power.
- [ ] Port from the detector is right-handed with `d·n > 0` for all hits, also for a tilted detector.

**Integration**
- [ ] `TestChain`: BMO beamlet → detector → `RayBundle` → summation → ASM (z = z_R/2) agrees with the analytic Gaussian (width `rtol 10⁻³`, power `rtol 10⁻⁶`) and with BMO's own detector field at that plane (`max|ΔE|/max|E| < 10⁻³`).
- [ ] Core suite passes on Julia 1.10 without BMO/WOP; integration suite passes on Julia 1; Aqua clean.
- [ ] Docs build; `[deps]` unchanged; BMO unchanged.

## Decisions

Guiding constraint: no change to BeamletOptics; WaveOpticsPropagation is third party.

**D1 — Converter API:** a) `convert_field(conv, field)` with `abstract type AbstractFieldConverter`; converters declare `input_representation`/`output_representation`, and `convert_field` runs `check_compatibility` first (chosen)

**D2 — Pure rays in the beamlet summation:** a) rejected with `MissingConverterError`: the summation needs Q per ray (chosen)

**D3 — Home of the BMO glue:** b) extension in OpticsBase (`OpticsBaseBeamletOpticsExt`, BMO as weakdep, compat `0.13.10`); amends settled decision 3 in `CLAUDE.md` (chosen)

**D4 — Home of the Fourier wrapper:** a) `FreeSpace`, `AngularSpectrumMethod` and the port-geometry checks in core; only the array call to WOP in `OpticsBaseWaveOpticsPropagationExt`; without WOP a `MethodError` with a hint to load it (chosen)

**D5 — Geometry supported by `AngularSpectrumMethod`:** a) coaxial only: output port parallel to the input port, same `u`, origin shifted by `z·n` (z of either sign), same refractive index; output grid = input grid; anything else throws `ArgumentError` (chosen)

**D6 — Tests for heavy extensions:** a) separate environment `test/integration/` with its own runner and a CI job on Julia 1; the core suite stays light and keeps running on 1.10 (chosen)

**D7 — BMO handover mechanism:** a) detector-based: after `solve_system!`, `RayBundle(detector)` reads the hits of a BMO `Detector`; the port is derived from the detector frame (chosen)

## Details

Common rules: follow `CLAUDE.md` and the conventions page; use only public OpticsBase API
plus the functions named here; do not edit files you do not own (wiring, `Project.toml`,
runtests and docs are pre-wired). Core tests: `Pkg.test(test_args=[...])`. Integration
tests: `julia --project=test/integration -e 'using Pkg; Pkg.instantiate(); include("test/integration/runtests.jl")' <TestName>`.

<details>
<summary>W2 — Angular-spectrum wrapper</summary>

- `src/FreeSpace.jl` (core, no WOP): `struct FreeSpace end` (homogeneous medium between two
  ports; the medium is the ports' refractive index). `AngularSpectrumMethod` with fields
  `padding::Bool`, `pad_factor::Int` (≥ 1), `bandlimit::Bool`, keyword constructor.
  Traits: input and output representation `SampledField{<:Any, 2}`.
- `__solve(prob::PropagationProblem{<:SampledField, FreeSpace}, alg::AngularSpectrumMethod)`:
  1. both ports `PlanarPort`; normals equal and `u` equal within `√eps`; else `ArgumentError`
  2. Δ = origin_out − origin_in, z = Δ·n; lateral part `‖Δ − z n‖ ≤ √eps·max(‖Δ‖, λ)`
  3. refractive indices equal; λₘ = λ/n
  4. L = `size(grid) .* spacing(grid)` (dim 1 ↔ L[1]; WOP's "y, x" names are labels only)
  5. per component c: `_angular_spectrum(E[:, :, c], z, λₘ, L, alg)` (copy the slice)
  6. return `PropagationSolution(SampledField(E_out, grid, port_out, λ), prob, alg)`
- `function _angular_spectrum end` in core without methods; the extension adds
  `OpticsBase._angular_spectrum(E::AbstractMatrix{<:Complex}, z, λₘ, L::NTuple{2}, alg)`
  using `WaveOpticsPropagation.AngularSpectrum(E, z, λₘ, L; padding, pad_factor, bandlimit)(E)`
  and returning a new array (the WOP result is a view into its buffer). Do not use
  `angular_spectrum` (it requires square fields).
- `_register_error_hints()` in `FreeSpace.jl`: registers a `MethodError` hint for
  `_angular_spectrum` telling the user to `using WaveOpticsPropagation`. The lead calls it
  from `OpticsBase.__init__` (pre-wired).
- WOP's kernel `exp(i k z √(1 − (λf)²))` with `k = 2π/λ` matches exp(−iωt); pass the medium
  wavelength λₘ. Its grids are centered at sample N÷2 + 1 (`CenterFT`), like `RegularGrid`.
- `TestFreeSpace` (core): all geometry errors (thrown before the backend is called), D = 3
  grid rejected by the representation check (`MissingConverterError`), and the hint text of
  the missing-backend `MethodError` via `sprint(showerror, e)`. Do not add backend methods
  in the core suite and never load WOP there; the successful path is tested in the
  integration environment.
- `TestWaveOpticsPropagation` (integration): see Acceptance W2; build the analytic Gaussian
  `SampledField` directly (no RayBundle).

</details>

<details>
<summary>W3 — BMO detector → RayBundle</summary>

Extension module `OpticsBaseBeamletOpticsExt` (uses `OpticsBase`, `BeamletOptics`,
`StaticArrays`, `LinearAlgebra`). Target the registered BMO 0.13.10; BMO internals
(`hits`, hit fields, `parabasal_ray_parameters`, `gauss_parameters`, `optical_power`,
`polarization`, `optical_path_length`) are allowed, BMO itself is not changed.

- `OpticsBase.RayBundle(det::BeamletOptics.Detector; port = nothing, power = nothing)`.
- Port: default `PlanarPort(position(det), n, u; refractive_index = n_ray)` with
  `n = −orientation(det)[:, 2]` (BMO's detector normal points against the beam) and
  `u = −orientation(det)[:, 1]` (BMO's local x). Verify `d·n > 0` for all hits. A given
  `port` must be a `PlanarPort` in the detector plane with the same normal. All hits must
  share λ and refractive index, else `ArgumentError`.
- `AstigmaticGaussianBeamletHit` → N = 3 beamlet:
  - position: chief hit point; direction: `hit.d0`; opl: optical path length of the chief
    beam up to the detector (parent included)
  - parabasal vectors at the hit plane: `h(l) = h + l·u` with l the chief segment length,
    from the cached `h1, u1, h2, u2` at the segment start
  - `(x̂, ŷ) = OpticsBase.ray_basis(port, d)`; H = [x̂·h1 x̂·h2; ŷ·h1 ŷ·h2], U likewise
    (non-conjugating products); `Q = U·H⁻¹`, symmetrized. BMO's u are geometric slopes, so
    Q fits the core convention with k = 2πn/λ; check Im Q > 0 at a waist.
  - phasor: `E_seg/‖E_seg‖` (polarization of the chief ray on the hit segment) times the
    unit Gouy factor of `√(area_ref/area_hit)`
  - power: `κ · A² · π/(k·√det(Im Q))` with `A = ‖E_seg‖·|√(area_ref/area_hit)|`; in free
    space this equals BMO's `optical_power(agb)`
- `GaussianBeamletHit` → N = 1 beamlet: `1/q = R + i λ/(π n w²)` from `gauss_parameters`
  (BMO's R is the curvature 1/r), `Q = (1/q)·I`; phasor `exp(i(arg E0 + ψ))` with BMO's
  Gouy phase ψ; power: BMO `optical_power(gauss)`; opl: `optical_path_length(gauss)`.
- `PolarizedRayHit` → N = 3 pure rays (phasor `E0/‖E0‖`); `RayHit` → N = 1 pure rays
  (phasor 1); opl from `hit.opl`; power: keyword `power` (total W) split equally, required
  for pure rays (`ArgumentError` otherwise).
- Never use BMO's detector field functions or its √cos projection factor.
- Tests: see Acceptance W3, plus a tilted detector (rotated about BMO's z axis).

</details>

<details>
<summary>Deviations</summary>

- W1: `AbstractFieldConverter` lives in `src/AbstractTypes.jl`, and the existing
  `is_compatible`/`check_compatibility` methods were widened to
  `Union{AbstractPropagationAlgorithm, AbstractFieldConverter}` instead of adding separate
  converter methods; `convert_field` also checks the output representation. The converter
  hook is `__convert_field` (public), mirroring `__solve`.
- W2: refractive indices are compared within `√eps`; the missing-backend hint is printed
  only while the extension is not loaded; `__solve` accepts any `AbstractGrid{2}`.
- W3: the beamlet power uses the exact Lagrange invariant, `κ·‖E‖²·π·|area_ref|/2`, so it
  equals BMO's `optical_power` to machine precision; BMO's parabasal rays fulfil the
  invariant only to O(θ²), so the amplitude derived from power and Q differs from BMO's
  field by ≈ 4·10⁻⁸. Extra checks: `power` must not be given for beamlets, the given port
  must match the hits' refractive index, empty detectors are rejected.
- Integration: the chain agrees with BMO's own field to 1.2·10⁻⁷, so `TestChain` checks
  `< 10⁻⁶` instead of the planned `10⁻³`. Added `examples/bmo_to_fourier.jl` (run in the
  integration CI job); the docs environment loads BeamletOptics so the extension
  docstring of `RayBundle(::Detector)` is embedded.
- Process: implementer worktrees were again cut from the last pushed commit; the
  implementers worked on an export of the local base commit.

</details>

<details>
<summary>Rationale</summary>

- Chain 1 uses free space only, so every stage has an analytic reference (Gaussian beam).
- The summation is the core converter from `CLAUDE.md` ("coherent summation of Gaussian
  beamlets onto a grid at the port"); it needs no dependency.
- D3 b) keeps BMO untouched (user constraint) at the price of OpticsBase tracking BMO
  internals; the compat bound `0.13.10` and the integration job detect breaking BMO
  releases.
- Rejected (rejected, do not implement): D1 b) converters as algorithms via `solve`,
  D1 c) callable converters; D2 b) default beamlet width for pure rays; D3 a) glue in BMO,
  D3 c) separate glue package; D4 b) types only in the extension, D4 c) upstream PR;
  D5 b) lateral shifts; D6 b) BMO/WOP in the main test targets; D7 b) full BMO solver
  wrapper with reverse conversion.
- Name clash to resolve before the `AngularSpectrum` exchange format is added:
  WaveOpticsPropagation exports `AngularSpectrum` (its propagator), so `using OpticsBase,
  WaveOpticsPropagation` would make both unusable unqualified.

</details>
