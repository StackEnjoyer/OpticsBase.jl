"""
    OpticsBase.Conformance

Shared conformance tests for packages that produce or propagate [`PlaneField`](@ref)s.
They compare a package's fields with analytic Gaussian beams and catch convention errors
that power and amplitude checks miss: a conjugated field (wrong sign of the time
convention), a wrong Gouy phase or curvature, a doubled or missing phase, swapped
polarization handedness, a misplaced plane or inconsistent `H`.

A package takes one or both roles:

- source: [`check_source`](@ref)`(make)` with `make(beam, origin, axes, size, spacing)`
  returning the package's `PlaneField` of the [`GaussianBeam`](@ref) `beam` on the given
  plane;
- propagator: [`check_propagator`](@ref)`(propagate)` with `propagate(f, L)` returning
  the field `f` propagated by `L` along its normal through a homogeneous medium of index
  `f.n`, on the plane `f.origin + L n` with the axes of `f`.

Both return a vector of [`Result`](@ref)s. In the package's test suite:

```julia
using OpticsBase: Conformance

@testset "\$(r.name)" for r in Conformance.check_source(make)
    @test r.value <= r.limit
end
```

One test set per check keeps the name of a failing check in the test output. The module
has no dependency on `Test`.
"""
module Conformance

using LinearAlgebra: cross, dot, norm, normalize
using StaticArrays: SMatrix, SVector
using ..OpticsBase: PlaneField, VACUUM_IMPEDANCE, _coordinates, backward, power,
                    reference_phase

# Default limit of the geometry checks; the orthonormality tolerance of `PlaneField`, so
# that axes computed in single precision pass
const GEOMETRY_TOL = 1e-6

"""
    GaussianBeam(; waist, direction, λ, w0, P = 1e-3, n = 1, jones = (1, 0))

Stigmatic, paraxial Gaussian beam used as the reference of the conformance tests.

# Arguments

- `waist`: global position of the waist center in \\[m\\].
- `direction`: propagation direction (normalized).
- `λ`: vacuum wavelength in \\[m\\].
- `w0`: waist radius (1/e² of the intensity) in \\[m\\].
- `P`: power in \\[W\\].
- `n`: refractive index of the medium.
- `jones`: Jones vector (normalized) in the basis `(u, v)` of the plane the beam is
  evaluated on.

# Conventions

With `z` the distance from the waist along `direction`, `ρ` the distance from the axis,
`k = 2π n/λ` and `z_R = π w0² n/λ`, the field is the `exp(−iωt)` Gaussian beam

`E = E₀ (w0/w) exp(−ρ²/w²) exp(i (k z − atan(z/z_R) + k ρ²/(2R)))`

with `w = w0 √(1 + (z/z_R)²)`, `1/R = z/(z² + z_R²)` and the real peak amplitude
`E₀ = √(4 Z₀ P/(n π w0²))`. The absolute phase is zero at the waist center.
"""
struct GaussianBeam
    waist::SVector{3, Float64}
    direction::SVector{3, Float64}
    λ::Float64
    w0::Float64
    P::Float64
    n::Float64
    jones::SVector{2, ComplexF64}
end

function GaussianBeam(; waist, direction, λ, w0, P = 1e-3, n = 1, jones = (1, 0))
    return GaussianBeam(SVector{3, Float64}(waist), normalize(SVector{3, Float64}(direction)),
        λ, w0, P, n, normalize(SVector{2, ComplexF64}(jones)))
end

rayleigh_range(b::GaussianBeam) = π * b.w0^2 * b.n / b.λ
beam_radius(b::GaussianBeam, z) = b.w0 * sqrt(1 + (z / rayleigh_range(b))^2)

"""
    field(beam::GaussianBeam, origin, axes, size, spacing)

The analytic [`PlaneField`](@ref) of `beam` on the plane with `origin`, `axes`
(columns `u`, `v`, `n`), `size = (nx, ny)` and `spacing`, built with the E-only
constructor. The plane normal `n` must be the beam direction to within `1e-6` (the
orthonormality tolerance of `PlaneField`, so axes computed in single precision are
accepted); the origin may lie off the axis.
"""
function field(b::GaussianBeam, origin, axes, size::NTuple{2, Integer}, spacing)
    A = SMatrix{3, 3, Float64}(axes)
    norm(A[:, 3] - b.direction) <= GEOMETRY_TOL ||
        throw(ArgumentError("the plane normal must be the beam direction"))
    E = samples(b, origin, A, size, spacing)
    return PlaneField(E, spacing, SVector{3, Float64}(origin), A, b.λ; n = b.n)
end

# `nx × ny × 2` array of the field of `b` along `u` and `v` at the samples of a plane. Each
# sample is evaluated at its position in space, with the direction of the beam and not
# the plane normal, so this holds on any plane. The reference of `field` and `compare`.
function samples(b::GaussianBeam, origin, axes, size, spacing)
    A = SMatrix{3, 3, Float64}(axes)
    u, v, d = A[:, 1], A[:, 2], b.direction
    o = SVector{3, Float64}(origin) - b.waist
    k, zR = 2π * b.n / b.λ, rayleigh_range(b)
    E0 = sqrt(4 * VACUUM_IMPEDANCE * b.P / (b.n * π * b.w0^2))
    ξs, ηs = _coordinates(Float64, size, spacing)
    ψ = [begin
             r = o + ξ * u + η * v
             z = dot(r, d)
             ρ2 = sum(abs2, r - z * d)
             w = beam_radius(b, z)
             E0 * b.w0 / w * exp(-ρ2 / w^2) *
             cis(k * z - atan(z, zR) + k * ρ2 * z / (2 * (z^2 + zR^2)))
         end
         for ξ in ξs, η in ηs]
    return cat(b.jones[1] .* ψ, b.jones[2] .* ψ; dims = 3)
end

"""
    Result(name, value, limit)

Outcome of one conformance check: the measured error `value` and its tolerance `limit`.
The check passes if `value <= limit`.
"""
struct Result
    name::String
    value::Float64
    limit::Float64
end

passed(r::Result) = r.value <= r.limit

function Base.show(io::IO, r::Result)
    verdict, relation = passed(r) ? ("pass", " <= ") : ("FAIL", " > ")
    print(io, verdict, ": ", r.name, " (", r.value, relation, r.limit, ")")
end

"""
    compare(f::PlaneField, beam::GaussianBeam; name = "", power_rtol = 1e-3,
            field_tol = 1e-4, phase_tol = 1e-3, backward_tol = 1e-6)

Compares `f` with the analytic field of `beam` on the plane of `f` and returns
[`Result`](@ref)s for:

- power: `|power(f)/P − 1|`;
- field: `1 − |c|`, with `c` the normalized complex correlation of the physical
  tangential fields (shape, position and polarization);
- absolute phase: `|arg c|` in rad (catches a conjugated field, a wrong Gouy phase or a
  doubled phase);
- backward light: `|power(backward(f))|/P` (`H` consistent with a forward wave).

The reference is the beam at the sample positions of `f`, with `jones` taken in the
`(u, v)` of `f`. The checks are meant for a plane normal to the beam. A field on another
plane does not throw: it fails the checks, and [`check_source`](@ref) and
[`check_propagator`](@ref) report the plane itself.
"""
function compare(f::PlaneField, b::GaussianBeam; name = "", power_rtol = 1e-3,
        field_tol = 1e-4, phase_tol = 1e-3, backward_tol = 1e-6)
    ref = samples(b, f.origin, f.axes, (size(f.E, 1), size(f.E, 2)), f.spacing)
    E = collect(f.E .* reference_phase(f))
    c = dot(ref, E) / (norm(ref) * norm(E))
    prefix = isempty(name) ? "" : name * ": "
    return [
        Result(prefix * "power", abs(power(f) / b.P - 1), power_rtol),
        Result(prefix * "field shape and polarization", 1 - abs(c), field_tol),
        Result(prefix * "absolute phase", abs(angle(c)), phase_tol),
        Result(prefix * "backward light", abs(power(backward(f))) / b.P, backward_tol)
    ]
end

# Plane axes (u, v, n) with normal `n` and `u` along the projection of `hint`
function plane_axes(n, hint)
    n = normalize(SVector{3, Float64}(n))
    u = normalize(hint - dot(hint, n) * n)
    return hcat(u, cross(n, u), n)
end

const λ_TEST = 1e-6
const W0_TEST = 20e-6

# The standard cases: (name, beam, z of the plane relative to the waist, lateral offset of
# the plane origin along u, propagation distance for propagators)
function cases(n)
    zR = π * W0_TEST^2 * n / λ_TEST
    oblique = normalize(SVector(1.0, 1.0, 1.0))
    return [
        ("axis +z, linear along u, plane before the focus",
            GaussianBeam(; waist = [0.0, 0, 0], direction = [0.0, 0, 1], λ = λ_TEST,
                w0 = W0_TEST, n), plane_axes([0.0, 0, 1], SVector(1.0, 0, 0)), -zR, 0.0,
            2zR),
        ("axis +y, right-circular, plane at the waist",
            GaussianBeam(; waist = [0.0, 0, 0], direction = [0.0, 1, 0], λ = λ_TEST,
                w0 = W0_TEST, n, jones = (1, im)), plane_axes([0.0, 1, 0], SVector(1.0, 0, 0)),
            0.0, 0.0, zR),
        ("oblique axis, linear at 45°, plane behind the focus and off the axis",
            GaussianBeam(; waist = [1e-3, 2e-3, 3e-3], direction = oblique, λ = λ_TEST,
                w0 = W0_TEST, n, jones = (1, 1)), plane_axes(oblique, SVector(0.0, 0, 1)),
            2zR, 10e-6, zR)
    ]
end

# Sampling of a plane at the distance z from the waist: 64 samples, 8 per beam radius
function plane(b::GaussianBeam, A, z, offset)
    o = b.waist + z * b.direction + offset * A[:, 1]
    Δ = beam_radius(b, z) / 8
    return o, (64, 64), (Δ, Δ)
end

"""
    check_source(make; n = 1, geometry_tol = 1e-6, kwargs...)

Runs the source conformance tests: for each standard [`GaussianBeam`](@ref) case,
`make(beam, origin, axes, size, spacing)` must return the field of `beam` on that plane as
a [`PlaneField`](@ref) with that origin, axes and sampling. The cases cover the
axes `+z`, `+y` (the BeamletOptics optical axis) and an oblique one, planes at the waist,
before and behind it and off the axis, linear and right-circular polarization
(`jones = (1, i)`, positive helicity), in a medium of index `n`. `kwargs` are tolerances,
see [`compare`](@ref). Returns the [`Result`](@ref)s, including a check that the returned
plane is the requested one.

`geometry_tol` is the limit of that check: the largest of the origin error in units of the
sample spacing, the error of the axes and the relative error of the spacing. The default
is the orthonormality tolerance of `PlaneField`, so axes computed in single precision
pass. A returned plane outside the limit fails this check; it does not throw.
"""
function check_source(make; n = 1, geometry_tol = GEOMETRY_TOL, kwargs...)
    results = Result[]
    for (name, b, A, z, offset, _) in cases(n)
        o, sz, sp = plane(b, A, z, offset)
        f = make(b, o, A, sz, sp)
        geometry = max(norm(f.origin - o) / sp[1], norm(f.axes - A),
            maximum(abs.(f.spacing ./ sp .- 1)), Float64(size(f.E)[1:2] != sz))
        push!(results, Result(name * ": plane as requested", geometry, geometry_tol))
        append!(results, compare(f, b; name, kwargs...))
    end
    return results
end

"""
    check_propagator(propagate; n = 1, geometry_tol = 1e-6, kwargs...)

Runs the propagator conformance tests: for each standard [`GaussianBeam`](@ref) case, the
analytic field on a plane normal to the beam is passed to `propagate(f, L)`, which must
return the field after the distance `L` along the plane normal in a homogeneous medium of
index `f.n`, on the plane `f.origin + L n` with the axes of `f` (sampling is up to the
propagator). The cases include propagation through a focus, circular polarization and an
oblique, off-axis plane. `kwargs` are tolerances, see [`compare`](@ref). Returns the
[`Result`](@ref)s, including a check of the output plane with the limit `geometry_tol`
(origin error in units of the input sample spacing and error of the axes, see
[`check_source`](@ref)).
"""
function check_propagator(propagate; n = 1, geometry_tol = GEOMETRY_TOL, kwargs...)
    results = Result[]
    for (name, b, A, z, offset, L) in cases(n)
        o, sz, sp = plane(b, A, z, offset)
        g = propagate(field(b, o, A, sz, sp), L)
        target = o + L * A[:, 3]
        geometry = max(norm(g.origin - target) / sp[1], norm(g.axes - A))
        push!(results, Result(name * ": output plane", geometry, geometry_tol))
        append!(results, compare(g, b; name, kwargs...))
    end
    return results
end

end # module Conformance
