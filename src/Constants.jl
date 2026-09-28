"""
    VACUUM_IMPEDANCE

Wave impedance of vacuum, Z₀ = μ₀c = 376.730313668 Ω (CODATA 2018).

The value is the same literal as `Z_vacuum` in BeamletOptics.jl, so that powers computed
in both packages agree to machine precision. CODATA 2022 gives 376.730313412 Ω; the
relative difference of 7·10⁻¹⁰ is below any practical relevance.
"""
const VACUUM_IMPEDANCE = 376.730313668

"""
    power_normalization(port::AbstractPort) -> Real

Returns the factor κ = n / (2 Z₀) in \\[1/Ω\\] that converts the squared field amplitude
into intensity at `port`, with `n = refractive_index(port)` and Z₀ =
[`VACUUM_IMPEDANCE`](@ref). The power through the port is

```math
P = \\kappa \\int |\\mathbf{E}|^2 \\,\\mathrm{d}A.
```

# Conventions

  - `E` is the physical complex field amplitude in \\[V/m\\] (peak, not RMS): the real
    field is `Re(E·exp(−iωt))`, hence the factor 1/2.
  - The formula is exact for fields that propagate along the port normal and a paraxial
    approximation otherwise; `E` carries no obliquity weighting.
  - κ is returned in the element type of the port.
"""
function power_normalization(port::AbstractPort)
    n = refractive_index(port)
    return n / (2 * oftype(n, VACUUM_IMPEDANCE))
end
