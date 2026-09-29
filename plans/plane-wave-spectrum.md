# PlaneWaveSpectrum and a high-NA focus chain

**Status:** implemented · **Tier:** full

## Goal

Add the exchange format `PlaneWaveSpectrum` with three converters (rays → plane waves by
Debye, sampled field → plane waves by FFT, plane waves → sampled field by summation) and a
second real chain: BeamletOptics rays → vectorial high-NA focus.

## Change map

| Path | Change |
|------|--------|
| `src/PlaneWaveSpectrum.jl` | new: format `PlaneWaveSpectrum` (type, invariants, traits, `total_power`) |
| `src/PlaneWaveSummation.jl` | new: converter `PlaneWaveSummation` (`PlaneWaveSpectrum` → `SampledField`, planar or volume grid, any port) |
| `src/PlaneWaveDecomposition.jl` | new: converter `PlaneWaveDecomposition` (`SampledField` → `PlaneWaveSpectrum`), FFT backend stub + error hint |
| `ext/OpticsBaseFFTWExt.jl` | new: centered FFT via FFTW |
| `src/DebyeWolf.jl` | new: converter `DebyeWolf` (`RayBundle` → `PlaneWaveSpectrum`), triangulation stub + error hint |
| `ext/OpticsBaseDelaunayTriangulationExt.jl` | new: dual cell areas via DelaunayTriangulation.jl (D5) |
| `test/TestPlaneWaveSpectrum.jl`, `TestPlaneWaveSummation.jl`, `TestPlaneWaveDecomposition.jl`, `TestDebyeWolf.jl` | new (core suite) |
| `test/integration/TestFocus.jl`, `examples/bmo_oap_focus.jl` | new: BMO focus chain (D7) |
| `src/OpticsBase.jl`, `Project.toml`, `test/runtests.jl`, `test/integration/{Project.toml, runtests.jl}`, CI | edit: wiring, weakdeps FFTW + DelaunayTriangulation, test targets, example step |
| `docs/src/{formats,converters,conventions,solvers}.md`, `CLAUDE.md` | edit: new format and converters, conventions for the spectrum |

## API delta

```diff
+ PlaneWaveSpectrum(port, wavelength, direction, amplitude, weight)  <: AbstractOpticalField
+ struct PlaneWaveSpectrum{N, T <: Real, P <: AbstractPort}
+     port::P; wavelength::T
+     direction::Vector{SVector{3,T}}           # unit s_j, global frame
+     amplitude::Vector{SVector{N,Complex{T}}}  # spectral density ℰ_j [V/m per sr], global frame
+     weight::Vector{T}                         # solid angle w_j [sr]
+ end
+ PlaneWaveSummation(grid::AbstractGrid, port = nothing)   <: AbstractFieldConverter
+ PlaneWaveDecomposition(; pad_factor = 1)                 <: AbstractFieldConverter
+ DebyeWolf(port::AbstractPort; solid_angle = nothing)     <: AbstractFieldConverter
```

| Symbol | Status | Breaking | Call sites |
|---|---|---|---|
| `PlaneWaveSpectrum` (+ `port`, `wavelength`, `total_power`, `is_vectorial`, `is_coherent`, `length`) | new | no | none |
| `PlaneWaveSummation`, `PlaneWaveDecomposition`, `DebyeWolf` (+ representation traits, `__convert_field`) | new | no | none |

Existing signatures (`convert_field`, `__convert_field`, `RayBundle(…)`,
`SampledField(…)`, `RayBundle(::Detector)`) are unchanged. `[deps]` is unchanged; FFTW and
DelaunayTriangulation become weakdeps.

## Logic

```text
┌──────────────────────────────┐        ┌──────────────────────────────┐
│ RayBundle                    │        │ SampledField                 │
│   rays at a pupil port       │        │   planar grid at a port      │
└──────────────┬───────────────┘        └──────────────┬───────────────┘
               │                                       │
               ▼                                       ▼
        DebyeWolf(port_F)                 PlaneWaveDecomposition()
               │                                       │
               └───────────────────┬───────────────────┘
                                   ▼
                    ┌──────────────────────────────┐
                    │ PlaneWaveSpectrum            │
                    │   samples at port F          │
                    └──────────────┬───────────────┘
                                   │
                                   ▼
                     PlaneWaveSummation(grid, port)
                                   │
                                   ▼
                    ┌──────────────────────────────┐
                    │ SampledField                 │
                    │   2D or 3D grid, any port    │
                    └──────────────────────────────┘
```

Model and formulas (port F: origin r₀, axes u, v, n, index n_F; λₘ = λ/n_F, k = 2π/λₘ,
k₀ = 2π/λ, κ = `power_normalization(port)`):

```text
Spectrum (D1, D2):
  E(r) = Σ_j w_j · ℰ_j · exp(i k s_j·(r − r₀))       s_j unit, s_j·n > 0
  P    = κ · λₘ² · Σ_j w_j · ‖ℰ_j‖²                   exact flux through a plane
  N = 3: s_jᵀ ℰ_j = 0                                 (no conjugation)

DebyeWolf, ray j (position p, direction d, opl, power P, phasor e) → sample j:
  s_j = d_j
  ℰ_j = −(i/λₘ) · √(P_j / (κ w_j)) · e_j · exp(i [k₀ opl_j − k d_j·(p_j − r₀)])
  w_j = A_j / (d_j·n),  A_j = dual cell area of (d_j·u, d_j·v)   (D5, or given)
  ⇒ Σ over the spectrum of P equals Σ_j P_j exactly, for any w
```

DebyeWolf (W3, lead):

1. `refractive_index(port_F) ≈ refractive_index(port(bundle))` (rtol √eps), else `ArgumentError`
2. per ray: `(r₀ − p_j)·d_j > 0` (converging, D6) and `d_j·n > 0`, else `ArgumentError`
3. `w = conv.solid_angle`, or: `a_j = d_j·u`, `b_j = d_j·v`; duplicates → `ArgumentError`
4. `A = _delaunay_dual_areas(a, b)` (extension); `w_j = A_j / (d_j·n)`
5. `length(w) == M`, all `w_j` finite and `> 0`, else `ArgumentError`
6. `ℰ_j` as above; return `PlaneWaveSpectrum(port_F, λ, d, ℰ, w)`

Trace (single ray, medium n = 1.5, λ = 1 µm, focus r₀ = 0, port A at z = −1 mm, normal z,
u = x): p = (0.5, 0, −1) mm, d = (−0.4472, 0, 0.8944), opl = 2 mm, P = 1 mW, w = 10⁻⁴ sr,
Jones (1, 0).

| step | value |
|---|---|
| e = `jones_to_global(port, d, (1, 0))` | (0.8944, 0, 0.4472), eᵀd = 0 |
| κ = 1.5/(2·376.730313668 Ω) | 1.990814·10⁻³ S |
| √(P/(κw)) / λₘ | 70.873 V / 0.6667 µm = 1.0631·10⁸ V/m per sr |
| k₀ opl − k d·(p − r₀) = k₀ (opl + n‖p‖) | k₀ · 3.6771 mm (total OPL to the focus) |
| extra phase | −π/2 (factor −i) |
| κ λₘ² w ‖ℰ‖² | 1.0000 mW |

## Workstreams

| W  | Owns | Needs | Executor | Summary |
|----|------|-------|----------|---------|
| W1 | `src/PlaneWaveSpectrum.jl`, `src/PlaneWaveSummation.jl`, `test/TestPlaneWaveSpectrum.jl`, `test/TestPlaneWaveSummation.jl` | – | opus | format + summation |
| W2 | `src/PlaneWaveDecomposition.jl`, `ext/OpticsBaseFFTWExt.jl`, `test/TestPlaneWaveDecomposition.jl` | W1 | opus | FFT decomposition |
| W3 | `src/DebyeWolf.jl`, `ext/OpticsBaseDelaunayTriangulationExt.jl`, `test/TestDebyeWolf.jl` | W1 | lead | Debye rays → plane waves, Richards–Wolf tests |
| W4 | `test/integration/TestFocus.jl`, `examples/bmo_oap_focus.jl` | W3 | lead | BMO off-axis parabola focus chain |
| –  | `src/OpticsBase.jl`, `Project.toml`, `test/runtests.jl`, `test/integration/{Project.toml, runtests.jl}`, `.github/workflows/CI.yml`, `docs/`, `CLAUDE.md` | all | lead | wiring, docs, conventions |

Before W1 starts, the lead pre-wires includes, exports, `__init__` hints, weakdeps,
extensions, test targets and stub files, and commits locally. Implementers work in the
main working copy (no worktrees; file sets are disjoint). The lead writes W3 while W1 runs.

## Acceptance

**W1**
- [ ] `TestPlaneWaveSpectrum`, `TestPlaneWaveSummation` pass.
- [ ] Each invariant (Details W1) throws `ArgumentError` when violated; `Float32` inputs give `T = Float32`; empty spectrum has `total_power == 0`.
- [ ] `total_power` of one sample equals `κ λₘ² w ‖ℰ‖²` within `rtol 10⁻¹⁴`.
- [ ] One plane wave, evaluated on a port tilted by 20° against n, on a 2D and a 3D grid, N = 1 and N = 3: every grid value equals `w ℰ exp(i k s·(r − r₀))` within `rtol 10⁻¹²`.
- [ ] Result is independent of the chunk size (M = 3·chunk + 1 against a naive double loop, `rtol 10⁻¹²`); the M × N_pts work is a matrix product, no `cis`/`exp` in that loop.
- [ ] Analytic Gaussian spectrum (w₀ = 20λ, Details W1): summation at the waist equals `E₀ exp(−ρ²/w₀²)` within `10⁻⁶·E₀`; `total_power` equals `κπw₀²E₀²/2` within `rtol 10⁻³`.
- [ ] Errors: other refractive index at the output port, `s·n_out ≤ 0` → `ArgumentError`; non-spectrum input → `MissingConverterError`.

**W2**
- [ ] `TestPlaneWaveDecomposition` passes.
- [ ] Plane wave exactly on an FFT bin (even and odd N, tilted in ξ and η, `pad_factor` 1 and 2): one sample carries `w ℰ = A` within `10⁻¹²·|A|`, all others `< 10⁻¹²·|A|`; its direction equals the analytic `s` within `10⁻¹⁴`; for N = 3 its vector is the analytic transverse amplitude.
- [ ] Round trip `SampledField` → decomposition → `PlaneWaveSummation(grid)` at the same port: scalar Gaussian (w₀ = 4λₘ, Δ = λₘ/4, 64 × 64): `max|ΔE|/max|E| < 10⁻¹⁰`, for `pad_factor` 1 and 2.
- [ ] Vectorial input with E_x only: spectrum passes the transversality invariant; the round trip returns E_x, E_y within `10⁻¹⁰`; decomposing the round-trip field again gives the same spectrum within `10⁻¹⁰` (idempotent).
- [ ] Paraxial Gaussian (w₀ = 10λₘ): `total_power(spectrum)` equals `total_power(field)` within `rtol 10⁻³`.
- [ ] Evanescent samples (k_t ≥ k) are dropped: sample count equals the analytic count for Δ = λₘ/4.
- [ ] Without FFTW the `MethodError` hint names FFTW; 3D fields → `MissingConverterError`; `pad_factor < 1` → `ArgumentError`.

**W3**
- [ ] `TestDebyeWolf` passes; the trace-table case is reproduced within `rtol 10⁻¹²`.
- [ ] `total_power(spectrum) == total_power(bundle)` within `rtol 10⁻¹²`.
- [ ] Richards–Wolf, NA 0.9, n = 1.5, x and radial polarization, rays on a Gauss–Legendre × uniform (θ, φ) grid with explicit solid angles: fields in the focal plane and at z = λₘ within 3λₘ of the axis match the 1D Bessel integrals within `10⁻⁸·max|E|`.
- [ ] Delaunay weights: rays on a regular (s_u, s_v) grid give interior weights `Δs²/s_n` within `rtol 10⁻¹²`; `Σ w_j s_n,j` equals the convex-hull area within `rtol 10⁻¹²`; a Fibonacci ray set reproduces the Richards–Wolf focus within `2·10⁻²·max|E|`.
- [ ] Errors: medium mismatch, non-converging ray, duplicate directions, fewer than 3 or collinear directions, wrong length or non-positive `solid_angle` → `ArgumentError`; without DelaunayTriangulation the `MethodError` hint names it.

**W4**
- [ ] `TestFocus` (integration): 90° off-axis parabola from BMO, 20 000 Fibonacci polarized rays, P = 1 mW → detector → `RayBundle` → `DebyeWolf` (Delaunay) → `PlaneWaveSummation` (64 × 64, ±3λ at the focus).
- [ ] Every BMO ray passes within `10⁻⁶·rfl` of the focus; optical path to the focus is equal for all rays within λ/1000.
- [ ] `total_power(spectrum)` = 1 mW within `rtol 10⁻¹²`; the focal field matches a reference from analytic mirror rays with exact solid angles (Gauss–Legendre pupil) within `2·10⁻²·max|E|`.
- [ ] `examples/bmo_oap_focus.jl` runs in the integration CI job and prints NA, peak intensity, longitudinal fraction and the deviation from the reference.

**Integration**
- [ ] Core suite (all modules incl. Aqua) and integration suite pass; docs build.
- [ ] `[deps]` unchanged; FFTW and DelaunayTriangulation are weakdeps with compat bounds and in the test targets; BeamletOptics unchanged.
- [ ] `docs/src/conventions.md` and `CLAUDE.md` describe the spectrum convention (D1, D2), the Debye phase factor −i (D6) and the new converters.

## Decisions

**D1 — Stored spectral quantity:** a) spectral density per solid angle ℰ in V/m, weights in sr; E(r) = Σ w ℰ exp(i k s·(r − r₀)) (chosen)

**D2 — `total_power` of a spectrum:** a) exact flux κ λₘ² Σ w ‖ℰ‖²; differs from the paraxial `SampledField` power by O(θ²), documented (chosen)

**D3 — Vectorial `SampledField` → spectrum:** a) the transverse components (E_u, E_v) at the port are authoritative; E_n of each plane wave follows from k·E = 0 (chosen)

**D4 — Evaluating a spectrum on a grid:** a) direct summation at any port (planar and 3D grids) as matrix products, exact, O(M · N_points); chirp-z later as a fast path; amends "Debye–Wolf via chirp-z" in `CLAUDE.md` (chosen)

**D5 — Solid angle per ray:** a) Delaunay dual cells in the direction plane via a DelaunayTriangulation.jl extension; explicit `solid_angle` vector as alternative (chosen)

**D6 — Scope of `DebyeWolf`:** a) converging bundles only (focus downstream of every ray), Debye factor −i; else `ArgumentError` (chosen)

**D7 — BeamletOptics chain:** a) 90° off-axis parabolic mirror, polarized Fibonacci rays, NA ≈ 0.4, analytic reference (chosen)

## Details

Common rules: follow `CLAUDE.md` and `docs/src/conventions.md`; model new code on
`src/RayBundle.jl` (constructor, invariants, promotion) and `src/BeamletSummation.jl`
(converter, docstring sections). Use only public OpticsBase API plus the functions named
here. Do not edit files you do not own; wiring, `Project.toml`, runtests and docs are
pre-wired. Julia ≥ 1.10 syntax only. Single-module tests:
`julia --project=. -e 'using Pkg; Pkg.test(test_args=["TestPlaneWaveSummation"])'`.

<details>
<summary>W1 — PlaneWaveSpectrum and PlaneWaveSummation</summary>

- `src/PlaneWaveSpectrum.jl`: struct as in API delta, inner constructor checks, outer
  positional constructor `PlaneWaveSpectrum(port, wavelength, direction, amplitude,
  weight)` that promotes all inputs to one float type `T` and dispatches on `N` behind a
  function barrier (as `RayBundle`). `amplitude` is a vector of numbers (N = 1) or of
  3-vectors (N = 3).
- Invariants, `tol = √eps(T)`, each violation an `ArgumentError` with the sample index:
  1. N ∈ (1, 3); all vectors equal length (M = 0 allowed)
  2. `wavelength > 0`; every weight finite and `> 0`
  3. `|‖s‖ − 1| ≤ tol`; `s·normal(port) > 0`
  4. N = 3: `|Σᵢ ℰᵢ sᵢ| ≤ tol·‖ℰ‖` (transversality without conjugation)
- Methods: `port`, `wavelength`, `is_vectorial` (N == 3), `is_coherent` (true),
  `Base.length`, `total_power = κ λₘ² Σ w ‖ℰ‖²` with `κ = power_normalization(port)`,
  `λₘ = λ/refractive_index(port)`; zero for M = 0.
- Docstring (self-sufficient, see CLAUDE.md): the model formula, units and frames of every
  field (ℰ in V/m per sr, global frame; w in sr; s unit, global; r₀ = `origin(port)`, the
  phase reference; medium = `refractive_index(port)`), exact power, the port's role
  (`n` bounds the hemisphere, `u` fixes Jones bases), `# Interface`.
- `src/PlaneWaveSummation.jl`: `struct PlaneWaveSummation{G <: AbstractGrid, P} <:
  AbstractFieldConverter` with fields `grid`, `port` (`nothing` = port of the spectrum).
  Grid D = 2 (ξ, η along u, v) or D = 3 (ζ along n). Traits: input `PlaneWaveSpectrum`,
  output `SampledField`.
- `__convert_field(conv, pws)`:
  1. output port `p` (conv.port or `port(pws)`); must be a `PlanarPort`, else `ArgumentError`
  2. refractive indices equal (rtol √eps), all `s_j·normal(p) > 0`, else `ArgumentError`
  3. with `o = origin(p)`: `a_j = s_j·(o − r₀)`, `b_j = s_j·u`, `c_j = s_j·v`, `d_j = s_j·n`
  4. per chunk of waves (fixed chunk size, e.g. 1024): `X[i,j] = cis(k ξ_i b_j)`,
     `Y[l,j] = cis(k η_l c_j)`, `A[j,c] = w_j ℰ_j[c] cis(k a_j)` (times `cis(k ζ d_j)` per ζ
     slice); accumulate `E[:, :, (ζ,) c] += X · Diagonal(A[:, c]) · Yᵀ` with `mul!`
  5. return `SampledField(E, grid, p, λ)`, `E` of size `(size(grid)..., N)`, eltype `Complex{T}`
- Docstring: no approximation beyond the discrete spectrum (exact superposition in the
  homogeneous medium); spectra from a grid give a periodic field; cost O(M · N_points).
- Tests. `TestPlaneWaveSpectrum`: every invariant, promotion, scalar and vector input,
  `total_power`. `TestPlaneWaveSummation`: Acceptance W1. Gaussian test: with
  `E(x) = ∫ Ẽ(k_t) exp(i k_t·x) d²k_t`, a waist `E₀ exp(−ρ²/w₀²)` at r₀ has
  `Ẽ = E₀ w₀²/(4π) · exp(−k_t² w₀²/4)`; sample it on a regular (k_u, k_v) grid
  (|k_t| ≤ 8/w₀, step ≤ 0.1/w₀) and set `ℰ = k k_n Ẽ`, `w = Δk_u Δk_v/(k k_n)`,
  `s = (k_u u + k_v v + k_n n)/k`.
- Do not implement FFT or Debye code here.

</details>

<details>
<summary>W2 — PlaneWaveDecomposition (FFT)</summary>

- `src/PlaneWaveDecomposition.jl`: `struct PlaneWaveDecomposition <:
  AbstractFieldConverter` with `pad_factor::Int` (≥ 1, `ArgumentError`), keyword
  constructor. Traits: input `SampledField{<:Any, 2}`, output `PlaneWaveSpectrum`.
- Backend stub `function _centered_fft end` (no methods in core). The extension
  `ext/OpticsBaseFFTWExt.jl` adds `OpticsBase._centered_fft(E::AbstractMatrix{<:Complex})
  = fftshift(fft(ifftshift(E)))` returning a new `Matrix`. Add `_register_fftw_hint()`: a
  `MethodError` hint for `_centered_fft` saying `using FFTW`, printed only while the
  extension is not loaded (see `_register_error_hints` in `src/FreeSpace.jl`); the lead
  calls it from `__init__`.
- `__convert_field(conv, field)` (formulas in Logic are for the spectrum; here the FFT):
  1. port `p` must be a `PlanarPort`; `λₘ = λ/refractive_index(p)`, `k = 2π/λₘ`
  2. N = 1: `E₁ = E[:, :, 1]`; N = 3: `E_u = Σ_c u_c E[:, :, c]`, `E_v` likewise (E_n unused, D3)
  3. zero-pad to `(M_ξ, M_η) = pad_factor .* size(grid)` so that input index `N÷2 + 1`
     lands on `M÷2 + 1` in each dimension
  4. `F = _centered_fft(E_pad) · Δξ Δη / (2π)²`; `k_ξ = (m − (M_ξ÷2 + 1))·Δk_ξ`,
     `Δk_ξ = 2π/(M_ξ Δξ)`; likewise η
  5. keep `k_t² < k²` (strict); `k_n = √(k² − k_t²)`; `s = (k_ξ u + k_η v + k_n n)/k`;
     `w = Δk_ξ Δk_η/(k k_n)`
  6. N = 1: `ℰ = k k_n F`; N = 3: `F_n = −(k_ξ F_u + k_η F_v)/k_n`,
     `ℰ = k k_n (F_u u + F_v v + F_n n)`
  7. return `PlaneWaveSpectrum(p, λ, s, ℰ, w)` in column-major (ξ, η) order
- Sign: the forward transform uses `exp(−i k_t·x)`, so that the summation with
  `exp(+i k s·r)` inverts it (exp(−iωt)). Check with the single-plane-wave test first.
- Docstring: assumes the field is periodic on the (padded) grid; evanescent components are
  dropped with their power; `total_power` of the spectrum is the exact flux and differs
  from the paraxial `SampledField` power by O(θ²); for N = 3 the normal component is
  recomputed from E_u, E_v.
- Tests: Acceptance W2. Test the hint first, guarded by
  `Base.get_extension(OpticsBase, :OpticsBaseFFTWExt) === nothing`, then `using FFTW`.
  Build all input fields analytically on the grid; the round trips use
  `PlaneWaveSummation` from W1.
- Do not touch `src/FreeSpace.jl` or `__init__`.

</details>

<details>
<summary>Deviations</summary>

- W1: `PlaneWaveSummation` rejects grids other than 2D/3D at construction; the result
  has the element type of the spectrum; the tests call the internal `_plane_wave_sum`
  with several chunk sizes. The docstring calls `total_power` the quadrature of the exact
  flux of the continuous spectrum (cross terms of distinct plane waves vanish only there).
- W2 (acceptance corrected, implementation as planned): with `pad_factor = p`, a plane
  wave filling the window becomes a sampled sinc, so its peak bin carries A/p², not A;
  the test checks that and the vanishing unpadded bins. The Gaussian round trip runs on
  192 × 192 (edges at ±6w₀) instead of 64 × 64: truncation at ±2w₀ puts content beyond
  k, which is correctly dropped as evanescent (error 1.5·10⁻³ there, 7·10⁻¹⁶ now).
  Fields at non-planar ports get an `ArgumentError`.
- W3: solid angles come from Voronoi cells clipped to the convex hull (the dual of the
  Delaunay triangulation). They are unique also for cocircular points such as regular
  grids, where barycentric Delaunay duals would depend on the arbitrary diagonals. The
  Richards–Wolf tests use sin α = 0.9 in n = 1.5 (NA 1.35). Measured: Richards–Wolf
  1.4·10⁻¹² (tolerance 10⁻⁸), E_z up to 36 % of max|E| (x polarization); Fibonacci rays
  with Voronoi weights 2.0·10⁻³ (tolerance 2·10⁻²).
- W4: measured NA 0.40 / 0.47 (u / v), rays within 0.16 nm of the focus (tolerance
  25 nm), OPL spread 1.5·10⁻⁴ λ (tolerance 10⁻³ λ), deviation from the reference
  2.0·10⁻³ (tolerance 2·10⁻²). Extra check: BMO's reflected polarization equals the
  ideal-mirror formula of the reference within 10⁻⁹.
- Integration: FFTW and DelaunayTriangulation are in the core test targets; `CLAUDE.md`
  documents this and the new conventions and converters.

</details>

<details>
<summary>Rationale</summary>

- Rejected (rejected, do not implement): D1 b) density per dk_u dk_v, D1 c) discrete
  amplitudes; D2 b) paraxial power; D3 b) 3D projection, D3 c) require consistency;
  D4 b) chirp-z now; D5 b) own Delaunay, D5 c) explicit weights only; D6 b) diverging and
  astigmatic bundles; D7 b) aspheric lens, D7 c) no BMO chain.
- Debye needs one solid angle per ray; a `RayBundle` only stores power, and BMO pure rays
  carry no tube size. Any ray → field conversion therefore needs a tessellation or known
  sampling. Delaunay dual cells cover the convex hull of the directions, so hull rays get
  about half a cell: the focus error scales as O(1/√M) (≈ 1 % at 10⁴ rays), the power is
  exact for any weights.
- The Debye factor −i is the focal phase anomaly of a converging wave (Born & Wolf 8.8,
  same time convention); diverging or astigmatic bundles need other factors (D6 b).
- Beamlet bundles are accepted by `DebyeWolf` as rays (chief ray, power); their `Q` is
  ignored and the docstring says so. For moderate NA the beamlet route
  (`GaussianBeamletSummation` at the focus) remains available.
- Richards–Wolf reference (W3), with A(θ) = −(i/λₘ)√(S₀ f² cos θ/κ) for a uniform
  entrance pupil (irradiance S₀, focal length f) and x = kρ sin θ:
  x-pol E_x = πA₀(I₀ + I₂ cos 2φ), E_y = πA₀ I₂ sin 2φ, E_z = −2πi A₀ I₁ cos φ with
  I₀ = ∫√cos θ sin θ (1 + cos θ) J₀ e^{ikz cos θ} dθ, I₁ = ∫√cos θ sin²θ J₁ …,
  I₂ = ∫√cos θ sin θ (1 − cos θ) J₂ …; radial pol E_ρ = 2πi A₀ ∫√cos θ cos θ sin θ J₁ …,
  E_z = −2π A₀ ∫√cos θ sin²θ J₀ …. Bessel functions by the trapezoid rule on
  J_n(x) = (1/π)∫₀^π cos(nτ − x sin τ) dτ (no SpecialFunctions dependency).
- W4 setup: BMO `OffAxisParabolicMirror(rfl, D; angle = 90)` reflects a beam along +y to
  the focus (−rfl sin 90°, −rfl cos 90°, 0); parent focal length rfl/2; the surface
  obeys ‖q − F‖ = rfl − q_y, so all rays have the same path to F and dΩ = dA_in/‖q − F‖².
  BMO's ideal mirror uses the Jones matrix diag(−1, 1) in (s, p); the reference uses the
  same. Detector between mirror and focus, clear of the incoming beam (D < rfl).

</details>
