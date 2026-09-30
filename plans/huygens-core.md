# Radical rewrite: one exchange format, the Huygens plane

**Status:** approved 2026-09-30; phases 1–4 implemented (uncommitted, local branches) · **Tier:** full · replaces all previous plans

## Goal

OpticsBase shrinks to a single exchange format that every Maxwell solver can produce and
consume: the tangential electric and magnetic field sampled on a plane in 3D space
(`PlaneField`). No solver interface, no traits, no converters, no ray formats. Solvers
couple by function composition: each package offers methods that return a `PlaneField`
and methods that take one.

```julia
f1 = PlaneField(detector; size = (256, 256), spacing = (0.25e-6, 0.25e-6))  # BeamletOptics
f2, stats = propagate(f1, fiber, FDBPM(dz = 1e-6); grid)                   # BeamletFibers
src = WavefrontBeamletDecomposition(f2)                                     # BeamletOptics
```

Only the logo (`docs/src/assets/logo.svg`) and the conventions that still apply
(time convention, SI, peak amplitudes, Z₀, helicity naming) survive.

## Why tangential E and H (D1)

The package should cover every solver of the ecosystem, not only one-way beam
propagators: beamlets (BeamletOptics), BPM (BeamletFibers), angular spectrum/Fresnel
(WaveOpticsPropagation), and later FDTD, FEM, RCWA, T-matrix/Mie. Candidates were
tangential E only `(Eu, Ev)`, global 3D E, and tangential E + H. The choice is
**tangential E and H in the local frame of the plane**, `(Eu, Ev, Hu, Hv)`:

1. **Exact and complete (surface equivalence theorem).** Tangential E and H on a surface
   determine the field on either side, for any medium and both directions of travel. E
   alone does so only for a field known to travel one way in a known homogeneous
   medium. Reflections, cavities, standing waves and bidirectional solvers
   (FDTD/FEM/RCWA) need H.
2. **No information lost.** The normal components follow locally from Maxwell's curl
   equations, without an FFT or a paraxial assumption:
   `En = i Z₀/(k₀ n²) (∂u Hv − ∂v Hu)`, `Hn = −i/(k₀ Z₀) (∂u Ev − ∂v Eu)`.
   A global 3D E array is therefore not more expressive, only redundant and
   one-way.
3. **Exact power.** The power through the plane is the Poynting flux
   `P = ½ Re ∫ (Eu Hv* − Ev Hu*) dA`, signed along `n`. It needs no refractive index
   and no paraxial `n/(2Z₀) ∫|E|²` factor.
4. **The native currency of other solvers.** FDTD TF/SF sources and DFT flux monitors,
   near-to-far-field transforms, FEM port boundaries and RCWA layer interfaces all work
   with tangential E and H.
5. **Local frame.** Moving or rotating the plane changes only `origin` and `axes`, never
   the arrays. BPM and Jones-type solvers already work in `(u, v)`.

The cost is falling on one-way producers that only know E. The constructor
`PlaneField(E, …)` handles that case by filling in H for a local plane wave along `n`,
and `forward(f)`/`backward(f)` let one-way consumers take the part that travels their way.

## API (complete)

```julia
const VACUUM_IMPEDANCE = 376.730313668   # Ω, same literal as BMO's Z_vacuum

struct PlaneField{T <: AbstractFloat, A <: AbstractArray{Complex{T}, 3}}
    E::A                       # nx × ny × 2: (Eu, Ev) in V/m, peak, full phase
    H::A                       # nx × ny × 2: (Hu, Hv) in A/m, peak, full phase
    spacing::NTuple{2, T}      # (Δu, Δv) in m
    origin::SVector{3, T}      # plane center, global frame, m
    axes::SMatrix{3, 3, T, 9}  # columns u, v, n; orthonormal, right-handed
    λ::T                       # vacuum wavelength, m
    n::T                       # refractive index of the medium at the plane
    R::T                       # radius of the reference sphere, m; Inf = none (D3)
end

PlaneField(E, H, spacing, origin, axes, λ; n = 1, R = Inf)  # validates
PlaneField(E, spacing, origin, axes, λ; n = 1, R = Inf)     # H of a local plane wave along n
                                                           # E: nx×ny×2, or nx×ny (scalar = u, D2)
coordinates(f, d)       # sample positions along u (d = 1) or v (d = 2), m
reference_phase(f)      # nx × ny factor exp(i k₀ n s (√(ρ² + R²) − |R|)), s = sign(R)
power(f)                # net Poynting flux along n, W
forward(f), backward(f) # parts travelling along ±n (local plane wave split, D4)
```

`A` is generic, so GPU arrays and `Float32` work. The only dependencies are
`StaticArrays` and `LinearAlgebra`; `CommonSolve` is dropped.

## Logic

- **Sample positions:** sample `i` sits at `(i − 1 − N÷2)·Δ` along its axis (fftshift
  center at the origin), as before. No grid offsets: move the origin instead.
- **Stored vs. physical field:** the physical tangential fields are
  `E .* reference_phase(f)` and `H .* reference_phase(f)`. With `R = Inf` the factor is 1.
  `R > 0`: diverging from the point `origin − R n`; `R < 0`: converging to
  `origin + |R| n`. The phase factors cancel in `power`, so `power` works on the stored arrays.
- **Local plane wave H:** `H = (n/Z₀) n̂ × E`, i.e. `Hu = −(n/Z₀) Ev`,
  `Hv = (n/Z₀) Eu`.
- **Split:** `E± = ½ (Et ∓ (Z₀/n) n̂ × Ht)`, i.e.
  `E+u = ½ (Eu + (Z₀/n) Hv)`, `E+v = ½ (Ev − (Z₀/n) Hu)`; H of each part from the
  local plane wave rule with sign ±. This is exact for waves along `n`, and a paraxial
  approximation otherwise. The exact split (per plane wave in the spectrum) needs an FFT
  and belongs in a solver.
- **Validation (constructor only):** equal sizes of E and H, third dimension 2,
  `axes` orthonormal with `det = +1`, `λ > 0`, `n > 0`, `R ≠ 0`.

## Decisions

**D1 — Field content:** tangential E and H in the local frame (see above). *(decided)*

**D2 — Scalar fields:** there is no scalar type. A scalar solver hands over `Eu = ψ`,
`Ev = 0`; a scalar consumer takes `Eu` or projects onto its polarization. This avoids a
second type and N² scalar/vector cases. *(recommended)*

**D3 — Reference sphere `R`:** the one non-minimal field. A converging beam of 25 mm
diameter at NA 0.1 and λ = 1 µm needs Δ < 5 µm, i.e. about 5000² samples × 4 components ×
16 B ≈ 1.6 GB, if the full phase is sampled. With the sphere removed, only the residual
aberration remains, typically a few hundred samples. This replaces the old
`PlaneWaveSpectrum` and `RayBundle` routes into a focus: a pupil field with `R = −f` is
exactly the input of a Debye–Wolf solver. A tilt term was left out on purpose: turn
the plane instead. *(recommended)*

**D4 — No solver interface:** `PropagationProblem`, `PropagationSolution`, CommonSolve,
`input_representation`/`output_representation` and `check_compatibility` are removed.
With one format there is nothing to check or convert. Each package names its own
verbs (`propagate`, `WavefrontBeamletDecomposition`, …). *(recommended)*

**D5 — Removed formats and numerics:** `RayBundle`, `PolarizedRayBundle`,
`PlaneWaveSpectrum`, `DebyeWolf`, `PlaneWaveSummation`, `PlaneWaveDecomposition`,
`FreeSpace`/`AngularSpectrumMethod` and all five extensions are removed. They are solver
numerics. Debye–Wolf and plane wave summation can return as a small solver package
that maps `PlaneField` → `PlaneField`. Consequence: direct ray exchange between ray
tracers (BMO ↔ OpticSim) is gone and goes through fields instead. *(recommended)*

**D6 — Glue location:** in the solver packages. BeamletFibers takes OpticsBase as a hard
dependency (as now). BeamletOptics gets `ext/BeamletOpticsOpticsBaseExt.jl` (weak
dependency); BMO may now be changed for this. The glue for third-party packages
(WaveOpticsPropagation, OpticSim) is dropped for now. *(recommended)*

**D7 — Beyond a single plane:** by convention, not by types. A closed Huygens box is a
`Vector{PlaneField}` (e.g. 6 faces around a Mie scatterer). Polychromatic light and pulses
(FDTD) are a `Vector{PlaneField}` over λ with a common time origin `t = 0`. Incoherent
light is a `Vector{PlaneField}` whose powers add. Curved ports are not supported.
*(recommended)*

**D8 — Kept conventions:** `exp(−iωt)`, SI units, peak amplitudes, Z₀ literal,
right-handed `(u, v, n)` with explicit `u`, absolute optical path in the phase, helicity
naming of circular polarization. *(decided)*

## Change map

| Path | Change |
|---|---|
| `src/*` | delete everything; new `src/OpticsBase.jl` (module + `PlaneField` + 5 functions, ~150 lines with docstrings) |
| `ext/*`, `examples/*`, `plans/*` (except this one), `test/integration/*` | delete |
| `Project.toml` | deps `StaticArrays`, `LinearAlgebra`; drop `CommonSolve`, all weakdeps/extensions; version 0.2.0 |
| `test/` | `TestPlaneField.jl`, `TestAqua.jl`, `runtests.jl` |
| `docs/src/` | keep `assets/logo.svg`; `index.md` (pitch + chain example), `conventions.md` (rewritten per D1–D8), `reference.md` (`@docs`); delete the rest |
| `CLAUDE.md`, `README.md` | rewrite to the new scope |
| BeamletFibers `src/OpticsBaseSolver.jl` | `propagate(f::PlaneField, system, alg; grid, facets)` → `(PlaneField, PropagationResult)`; entry via `forward(f)`, exit via the E-only constructor |
| BeamletOptics `ext/BeamletOpticsOpticsBaseExt.jl` | new: `PlaneField(detector; size, spacing)` (E from `electric_field`, H as Σ over hits of `(n/Z₀) dᵢ × Eᵢ` with the beamlet direction `dᵢ`), `WavefrontBeamletDecomposition(::PlaneField)` via `forward(f)` |

## Phases

1. OpticsBase core (this repo), tests, docs. *Done.*
2. BeamletFibers glue on the new core; the existing core tests move over. *Done*
   (branch `opticsbase-planefield`): `propagate(f::PlaneField, segments, alg; grid, …)
   -> (PlaneField, PropagationResult)`, exit plane coaxial at the total length.
3. BeamletOptics extension. *Done* (branch `feature/opticsbase-ext`):
   `PlaneField(detector; size, spacing)` and `WavefrontBeamletDecomposition(::PlaneField)`.
4. Chain test BMO → fiber → BMO. *Done:* `test/integration/TestChain.jl` in
   BeamletFibers.

### Findings from phases 3 and 4

- BMO's detector sum multiplies each beamlet by `√|cos θ|` so that `|E|²` integrates
  to the power per detector area. That factor is not part of the physical field; the
  `PlaneField` glue leaves it out, and the Poynting flux of an oblique beam is then the
  full beam power.
- Two phase bugs in the first glue draft, both invisible to power and amplitude tests:
  the projection `dot(E, u)` conjugated the field (Julia's `dot` conjugates its first
  argument), and the complex beamlet amplitude entered twice (once in the scalar
  field, once in the polarization vector). Tests now compare complex fields, including
  the absolute phase, against an independent solver.
- The chain uses a large-core fiber (200 µm core, 60 µm beam). BMO's wavefront
  decomposition places one beamlet per sample with a waist of 1.2 samples, which needs
  beamlets much larger than λ; a single-mode fiber mode (w ≈ 5 µm) is too small for
  that return path. The decomposition smooths the field with that Gaussian, which keeps
  the power fraction `w²/(w² + w_b²)` of a Gaussian beam.
- Result: BMO beamlets and BeamletFibers' angular spectrum solver, both started from the
  fiber exit field, agree after 20 mm of air with a complex correlation of 0.996 and a
  phase offset of 0.006 rad.

## Acceptance (phase 1)

- [x] Gaussian beam built with the E-only constructor: `power(f)` equals the analytic
      power (`rtol 1e-9`), `forward(f) ≈ f`, `backward(f) ≈ 0`.
- [x] Two equal counter-propagating plane waves: `power(f) ≈ 0`, and `forward`/`backward`
      recover both parts.
- [x] `power(backward-travelling wave) < 0`.
- [x] Rigid motion: a new `origin`/`axes` leaves `E`, `H` and `power` untouched.
- [x] `reference_phase`: value at the center 1, and at radius ρ it equals the analytic
      spherical phase for both signs of `R`; `power` is independent of `R`.
- [x] `PlaneField` with a scalar matrix puts it into `Eu`.
- [x] Constructor errors: size mismatch, wrong component count, non-orthonormal or
      left-handed axes, `λ ≤ 0`, `n ≤ 0`, `R = 0`.
- [x] `Float32` and a non-`Array` `AbstractArray` work.
- [x] Aqua passes; docs build.
