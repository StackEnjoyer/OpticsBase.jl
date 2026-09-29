# Exchange formats as a hierarchy of fields and rays

**Status:** implemented · **Tier:** full

## Goal

Exchange formats form one hierarchy (Maxwell fields, rays); OpticsBase knows no solver
(no beamlets, no BMO detector); ray exchange is tested between BMO and OpticSim.jl.

## Change map

| Path | Change |
|------|--------|
| `src/AbstractTypes.jl` | edit: root type (D8), `AbstractOpticalField{N}` + aliases, `AbstractRayBundle` + interface |
| `src/RayBundle.jl` | edit: minimal `RayBundle` (port, λ, position, direction); `PolarizedRayBundle` (D6) |
| `src/BeamletSummation.jl`, `test/TestBeamletSummation.jl` | delete |
| `src/DebyeWolf.jl`, `test/TestDebyeWolf.jl` | edit: input `PolarizedRayBundle` (D6) |
| `src/SampledField.jl`, `src/PlaneWaveSpectrum.jl` | edit: supertype `AbstractOpticalField{N}`, generic `is_vectorial` |
| `src/Problem.jl`, `src/Converters.jl` | edit: accept the root type (fields and rays) |
| `src/OpticsBase.jl`, `Project.toml`, `test/runtests.jl`, `test/integration/runtests.jl` | edit: exports, weak deps, module lists (integration) |
| `ext/OpticsBaseBeamletOpticsExt.jl` | edit: BMO-side glue: rays in/out, `SampledField(detector, grid)` via `BMO.electric_field` |
| `ext/OpticsBaseOpticSimExt.jl` | new: OpticSim-side glue: rays in/out (temporary, like BMO) |
| `test/integration/Project.toml` | edit: OpticSim from Git, pinned commit |
| `test/integration/TestOpticSim.jl`, `TestRayExchange.jl` | new |
| `test/integration/TestBeamletOptics.jl`, `TestChain.jl`, `TestFocus.jl`, both examples | edit |
| `docs/src/*.md`, `README.md`, `CLAUDE.md` | edit |

## API delta

```diff
+ abstract type AbstractOpticalData end                  # name: D8; interface port, wavelength
- abstract type AbstractOpticalField end
+ abstract type AbstractOpticalField{N} <: AbstractOpticalData end
+ const AbstractVectorField = AbstractOpticalField{3}
+ const AbstractScalarField = AbstractOpticalField{1}
+ is_vectorial(::AbstractOpticalField{N}) where {N} = N == 3
+ abstract type AbstractRayBundle <: AbstractOpticalData end   # + position, direction, length

- struct SampledField{N, D, T, A, G, P} <: AbstractOpticalField
+ struct SampledField{N, D, T, A, G, P} <: AbstractOpticalField{N}
- struct PlaneWaveSpectrum{N, T, P} <: AbstractOpticalField
+ struct PlaneWaveSpectrum{N, T, P} <: AbstractOpticalField{N}

- struct RayBundle{N, T, P, B} <: AbstractOpticalField
- RayBundle(port, wavelength, position, direction, opl, power, phasor; beamlet = nothing)
+ struct RayBundle{T, P} <: AbstractRayBundle
+ RayBundle(port::AbstractPort, wavelength::Real, position::AbstractVector, direction::AbstractVector)
+ struct PolarizedRayBundle{T, P} <: AbstractRayBundle
+ PolarizedRayBundle(port::AbstractPort, wavelength::Real, position::AbstractVector,
+     direction::AbstractVector, opl::AbstractVector, power::AbstractVector, polarization::AbstractVector)
- has_beamlets(::RayBundle); is_coherent(::RayBundle); is_vectorial(::RayBundle); total_power(::RayBundle)

- struct GaussianBeamletSummation{G <: AbstractGrid{2}} <: AbstractFieldConverter   (+ 3 methods)
  DebyeWolf(port; solid_angle = nothing)                    # input now PolarizedRayBundle, output PlaneWaveSpectrum{3}

- PropagationProblem{F <: AbstractOpticalField, S, P}
+ PropagationProblem{F <: AbstractOpticalData, S, P}
- convert_field(conv::AbstractFieldConverter, field::AbstractOpticalField)
+ convert_field(conv::AbstractFieldConverter, data::AbstractOpticalData)

  # ext, BMO-side glue (BeamletOptics loaded)
- RayBundle(detector::BMO.Detector; port = nothing, power = nothing)     # N, beamlets
+ RayBundle(detector::BMO.Detector; port = nothing)                      # geometry only
+ PolarizedRayBundle(detector::BMO.Detector; port = nothing, power)      # D6a
+ SampledField(detector::BMO.Detector, grid::AbstractGrid{2}; port = nothing)   # BMO field functions; N = 1 stigmatic, N = 3 astigmatic
+ BMO rays from a RayBundle (D7; name fixed by W2, e.g. `BMO.Ray` vector)
  # ext, OpticSim-side glue (OpticSim loaded)
+ RayBundle(rays, port, wavelength)  and  OpticSim rays from a RayBundle   (names fixed by W3)
```

| Symbol | Status | Breaking | Call sites |
|---|---|---|---|
| `AbstractOpticalData` | new, exported | no | `Problem.jl`, `Converters.jl` |
| `AbstractOpticalField{N}`, `AbstractVectorField`, `AbstractScalarField` | changed / new | yes for external subtypes (none known) | – |
| `AbstractRayBundle` | new, exported | no | – |
| `RayBundle` (7-arg constructor, `N`, `B`) | narrowed to 4 args | yes | tests, ext, examples |
| `has_beamlets`, `is_coherent/is_vectorial/total_power(::RayBundle)` | removed | yes | `TestRayBundle` |
| `GaussianBeamletSummation` | removed, export removed | yes | `TestBeamletSummation`, `TestChain`, example, docs |
| `PolarizedRayBundle` | new, exported | no | – |
| `PropagationProblem`, `convert_field` | widened to root type | no | – |
| `RayBundle(::BMO.Detector)` | geometry only, `power` kwarg removed | yes | integration tests, examples |
| `SampledField(::BMO.Detector, grid)` | new (ext) | no | – |

Type hierarchy afterwards (recommended options):

```text
┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐
│ SampledField    │ │ PlaneWave-      │ │ RayBundle       │ │ PolarizedRay-   │
│   {N, D, …}     │ │   Spectrum{N,…} │ │   pos, dir, λ   │ │   Bundle (D6)   │
└────────┬────────┘ └────────┬────────┘ └────────┬────────┘ └────────┬────────┘
         │                   │                   │                   │
         └─────────┬─────────┘                   └─────────┬─────────┘
                   ▼                                       ▼
     ┌───────────────────────────┐           ┌───────────────────────────┐
     │ AbstractOpticalField{N}   │           │ AbstractRayBundle         │
     │   Maxwell, N = 1 or 3     │           │   geometrical optics      │
     └─────────────┬─────────────┘           └─────────────┬─────────────┘
                   │                                       │
                   └───────────────────┬───────────────────┘
                                       ▼
                        ┌─────────────────────────────┐
                        │ AbstractOpticalData (D8)    │
                        │   port, wavelength          │
                        └─────────────────────────────┘
```

## Workstreams

| W  | Owns | Needs | Executor | Summary |
|----|------|-------|----------|---------|
| W1 | `src/AbstractTypes.jl`, `src/RayBundle.jl`, `src/BeamletSummation.jl`, `src/DebyeWolf.jl`, `src/SampledField.jl`, `src/PlaneWaveSpectrum.jl`, `src/Problem.jl`, `src/Converters.jl`, core tests in the change map | – | lead | hierarchy, minimal rays, Debye per D6 |
| W2 | `ext/OpticsBaseBeamletOpticsExt.jl`, `test/integration/TestBeamletOptics.jl`, `TestChain.jl`, `TestFocus.jl`, both examples | W1 | sonnet | BMO-side glue |
| W3 | `ext/OpticsBaseOpticSimExt.jl`, `test/integration/Project.toml`, `TestOpticSim.jl`, `TestRayExchange.jl` | W1; W2 for `TestRayExchange` | opus | OpticSim glue and BMO ↔ OpticSim test |
| W4 | `docs/src/*.md`, `README.md`, `CLAUDE.md` | W1 | sonnet | docs, README roles, settled decisions |
| –  | `src/OpticsBase.jl`, `Project.toml`, `test/runtests.jl`, `test/integration/runtests.jl`, `.github/workflows/CI.yml` | all | lead | exports, weak deps, module lists, full suite |

## Acceptance

**W1**
- [ ] `TestRayBundle`, `TestSampledField`, `TestPlaneWaveSpectrum`, `TestConverters`, `TestProblem` and `TestDebyeWolf` pass.
- [ ] Tests: `SampledField{3} <: AbstractVectorField`, `PlaneWaveSpectrum{1} <: AbstractScalarField`, `RayBundle <: AbstractRayBundle <: AbstractOpticalData`, `!(RayBundle <: AbstractOpticalField)`; `PropagationProblem(bundle, sys, port)` constructs.
- [ ] `RayBundle` keeps the geometric checks of today (unit directions, `n·d > 0`, positions on the port plane) with the same tolerances.
- [ ] `TestDebyeWolf` numerical assertions keep their tolerances; only the input type changes.
- [ ] No identifier `beamlet`, `Q` or `BeamletOptics` in `src/` (grep).

**W2**
- [ ] Integration modules `TestBeamletOptics`, `TestChain`, `TestFocus` pass.
- [ ] `TestChain` (astigmatic beamlet): `SampledField(detector, grid) isa SampledField{3}` with `total_power ≈ P0` (rtol 10⁻⁶) and beam radius vs. analytic Gaussian (rtol 10⁻³) at both ports, as today.
- [ ] Single astigmatic beamlet: `‖E‖` equals `|BMO.electric_field(detector)|` within 10⁻¹² relative; a stigmatic beamlet gives `SampledField{1}`.
- [ ] BMO → `RayBundle` → BMO rays → `RayBundle` reproduces positions and directions within 10⁻¹².
- [ ] Both examples run.

**W3**
- [ ] `test/integration` instantiates on Julia 1.12+ with BMO and OpticSim together.
- [ ] `TestOpticSim`: OpticSim → `RayBundle` → OpticSim reproduces positions and directions within 10⁻¹².
- [ ] `TestRayExchange`: a ray fan traced by BMO to plane A, handed over as `RayBundle`, traced by OpticSim through a spherical lens to plane B, agrees with BMO tracing the same lens from A to B: positions within 10⁻⁹ m, directions within 10⁻¹². Plus the reverse direction (OpticSim → BMO).

**W4**
- [ ] `grep -rn "GaussianBeamletSummation\|has_beamlets\|sum_beamlets" docs README.md CLAUDE.md` finds nothing.
- [ ] Docs build passes.

**Integration**
- [ ] Full core suite incl. `TestAqua`, all integration modules and the docs build pass.

## Decisions

**D1 — Scalar/vector axis in the type:** a) `AbstractOpticalField{N}` with `AbstractVectorField`, `AbstractScalarField` (chosen)

**D2 — Rays:** `AbstractRayBundle` with the concrete `RayBundle` (port, position, direction, wavelength) and a type for polarized rays; ray exchange tested between BMO and OpticSim.jl (chosen)

**D3 — Beamlet summation:** removed from OpticsBase; OpticsBase does not know beamlets (chosen)

**D4 — BMO handover:** implemented BMO-side, for now in `OpticsBaseBeamletOpticsExt` (chosen)

**D5 — Summation of BMO's beamlets:** done with BMO's own field functions; stigmatic `GaussianBeamlet` hits give a scalar `SampledField{1}`, `AstigmaticGaussianBeamlet` hits a vectorial `SampledField{3}` (BMO carries the polarization there) (chosen)

**D6 — Debye–Wolf:** a) `PolarizedRayBundle` as a full type now (port, λ, position, direction, OPL, power per ray, unit transverse 3D polarization vector per ray; today's vectorial `RayBundle` without beamlets); `DebyeWolf` stays a generic converter `PolarizedRayBundle → PlaneWaveSpectrum{3}`; scalar Debye tests are dropped (chosen)

**D7 — Direction of the ray exchange test:** a) both ways, BMO → OpticSim and OpticSim → BMO (chosen)

**D8 — Name of the common supertype:** a) `AbstractOpticalData` (chosen)

**D9 — WaveOpticsPropagation glue:** a) not part of this plan; follow-up plan (chosen)

## Details

<details>
<summary>W2 — BMO-side glue</summary>

- The extension is BMO's code that lives here temporarily. It may use BMO internals; it
  must not add solver numerics to `src/`.
- `RayBundle(detector; port)`: positions and directions of the hits, wavelength (common to
  all hits, as today). Drop OPL, power, phasor and beamlet mapping from this constructor.
- `PolarizedRayBundle(detector; port, power)` (D6): today's `_raybundle` for
  `PolarizedRayHit` without beamlets; BMO rays carry no power, so `power` (total, W) is
  required, as in today's `RayBundle(det; power)`.
- `SampledField(detector, grid; port)`, dispatch on the hit type:
  - `GaussianBeamletHit` (stigmatic): `BMO.electric_field(detector; n, x_min, x_max, z_min,
    z_max)` → `SampledField{1}`.
  - `AstigmaticGaussianBeamletHit`: `SampledField{3}`. BMO's detector method sums only the
    scalar part (`E_ref_amp = norm(polarization)`); the polarization lives on the beamlet
    (`BMO.polarized_field(agb, r, z)`, `polarization(chief)`). Sum per hit: scalar detector
    field of that single hit (`BMO.electric_field(detector, [hit]; …)`, keeps BMO's
    obliquity and phase terms) times the hit's normalized polarization
    `polarization(chief) / hit.E_ref_amp` in the global frame. Verify on a single beamlet
    that `‖E‖` equals `|BMO.electric_field(detector)|` within 10⁻¹² relative and that the
    direction equals `BMO.polarized_field` at the grid points.
  - Pure-ray hits (unscaled field in BMO): `ArgumentError`.
  - Sampling: BMO's samples must coincide with the grid samples. Mapping as in
    `examples/bmo_to_fourier.jl`: `ξ = x`, `η = −z` (reverse the second axis). Require a
    square `n × n` grid with equal spacing (BMO takes one `n`), otherwise `ArgumentError`.
    Our sample `i` is at `(i − (n÷2 + 1))·Δ`, so `x_min = −(n÷2)·Δ`,
    `x_max = (n − 1 − n÷2)·Δ`.
  - Check BMO's amplitude normalization against `power_normalization(port)`; if it
    differs, rescale and document the factor in the docstring.
- BMO rays from a `RayBundle` (D7): construct BMO ray objects at the bundle positions and
  directions with its wavelength; pick the BMO type that `solve_system!` accepts.
- `TestChain`: replace the converter call by `SampledField(detector, grid)`; keep the
  analytic checks, drop the self-comparison with `BMO.electric_field`.
  `TestFocus`/OAP example: `convert_field(DebyeWolf(portF), PolarizedRayBundle(detector; power = P0))`.
- Docstrings: self-sufficient (units, frames, normalization, approximations), and state
  that the method is BMO-side glue that will move into BeamletOptics.
- Do not touch `src/`, `Project.toml` or the module lists. Run only the integration
  modules above.

</details>

<details>
<summary>W3 — OpticSim glue and ray exchange test</summary>

- OpticSim 0.7.1 (registered) pins JET 0.9 and therefore Julia 1.11; BMO needs Julia
  ≥ 1.12. OpticSim's `main` branch supports Julia 1.12 (same version number). Add it to
  `test/integration/Project.toml` via `[sources]` with the Git URL
  `https://github.com/brianguenter/OpticSim.jl` and a pinned `rev` (commit hash, not
  `main`).
- `ext/OpticsBaseOpticSimExt.jl` (trigger `OpticSim`; the lead adds the weak dep): map
  OpticSim rays at a plane to `RayBundle` and back. Read OpticSim's ray and tracing API
  (`OpticalRay`, `trace`, sources/emitters) first; units: OpticSim works in mm, OpticsBase
  in m — convert explicitly and test it.
- `TestOpticSim`: round trip OpticSim → `RayBundle` → OpticSim, positions and directions
  within 10⁻¹².
- `TestRayExchange` (needs W2's BMO mapping): a fan of about 7 collimated rays; BMO traces
  to plane A, hand over, OpticSim traces through a spherical singlet to plane B; reference:
  BMO traces the same singlet from A to B. Use the same constant refractive index in both
  tracers (if OpticSim needs a glass, evaluate its index at λ and pass that constant to
  BMO). Reverse direction as well (D7).
- Do not touch `src/`, other extensions, `Project.toml` or the module lists.

</details>

<details>
<summary>W4 — Docs, README, CLAUDE.md</summary>

- `docs/src/formats.md`: hierarchy with root (D8), fields (`AbstractOpticalField{N}`,
  aliases, `SampledField`, `PlaneWaveSpectrum`) and rays (`AbstractRayBundle`,
  `RayBundle`, `PolarizedRayBundle`); every exported name in an `@docs` block.
- `docs/src/converters.md`: drop `GaussianBeamletSummation`; `DebyeWolf` per D6.
- `docs/src/solvers.md`: BMO and OpticSim glue as solver-side code living here
  temporarily; show `SampledField(detector, grid)` and the ray handover.
- `docs/src/conventions.md`: "Physical validity at the port" — rays are handed over only
  between ray tracers or into a ray-to-field converter, never at caustics.
- `README.md`, "Connecting a solver package": ray tracers hand over rays
  (`AbstractRayBundle`) to other ray tracers and fields to field solvers; OpticsBase knows
  no solver.
- `CLAUDE.md`: "Exchange formats" (hierarchy), "Converters (initial scope)" (no beamlet
  summation), settled decision 3 (OpticSim glue as further temporary exception). Keep the
  conventions summary in sync with `conventions.md`.
- Escape unit brackets as `\[m\]`. Do not touch `src/`, `ext/`, `test/`.

</details>

<details>
<summary>Deviations</summary>

- W1: the ray accessors are `positions(bundle)` and `directions(bundle)` (plural, public,
  not exported), because BMO exports `position`. `src/Constants.jl` still names BMO's
  `Z_vacuum` as the source of the Z₀ literal (a convention statement, not code).
- W2: `PolarizedRayBundle(detector; power)` splits the power equally over the rays (BMO
  rays carry none), as the old `RayBundle(det; power)` did. The D7 method is
  `BeamletOptics.Beam(bundle::RayBundle, i)` (one BMO beam per ray). `SampledField(detector,
  grid; port)` requires a port with the detector axes. BMO's amplitude normalization equals
  `power_normalization`; no rescale.
- W3: `TestRayExchange` checks directions within 10⁻¹⁰ instead of 10⁻¹² (measured
  6·10⁻¹¹): BMO finds lens surfaces by ray marching to ~10⁻¹⁰ m, while OpticSim agrees with
  an exact analytic trace to 10⁻¹⁶. Positions agree within 9·10⁻¹³ m. OpticSim is pinned to
  commit `9a03f0b`; resolving it downgraded ForwardDiff, IntervalArithmetic,
  LogExpFunctions and OrderedCollections in the (untracked) integration Manifest.
- W4: no `@docs` entries for the OpticSim methods; the docs environment does not load
  OpticSim.

</details>

<details>
<summary>Rationale</summary>

- Exchanged is what is physical and solver-independent: Maxwell fields on a screen, and
  rays as the geometrical-optics limit. Beamlets are BMO's internal decomposition; their
  coherent sum is BMO's job (D3, D5).
- The hierarchy has two branches because ray-to-ray exchange is wanted (BMO ↔ OpticSim);
  `PropagationProblem` and `convert_field` accept the common root so a ray tracer can be a
  propagation stage with rays in and out.
- Scalar/vector is the parameter `N` of the field branch (D1). Julia subtyping means "is
  a", so a poorer representation is not modeled as a subtype of a richer one.
- D5 consequence: vectorial handover from BMO needs astigmatic beamlets; stigmatic
  `GaussianBeamlet`s give a scalar field. The vectorial beamlet summation of OpticsBase
  is dropped; the glue combines BMO's scalar per-hit field with BMO's polarization.
- D6 rejected (do not implement): `PolarizedRayBundle` as stub with Debye–Wolf moved into
  the BMO extension; stub with Debye–Wolf and the OAP chain removed.
- D7 rejected (do not implement): exchange test only BMO → OpticSim.
- D8 rejected (do not implement): `AbstractExchangeFormat`; `AbstractOpticalField` as root
  with `AbstractWaveField{N}` for fields.
- D9 rejected (do not implement): moving `FreeSpace`/`AngularSpectrumMethod` in this plan.
  They stay in core as a known exception until the follow-up plan.
- Rejected earlier in this plan (do not implement): `sum_beamlets`/`debye_wolf` as helper
  functions in OpticsBase; `BeamletField` as exchange format; `AbstractWaveField{N}` layer
  under a field-only root.

</details>
