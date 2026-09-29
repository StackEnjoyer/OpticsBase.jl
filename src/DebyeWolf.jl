"""
    DebyeWolf(port::AbstractPort; solid_angle = nothing) <: AbstractFieldConverter

Converter from a converging [`PolarizedRayBundle`](@ref) to a vectorial
[`PlaneWaveSpectrum`](@ref) at `port`, in the Debye approximation: every ray becomes one
plane wave along its direction, the origin of `port` is the focus (phase reference) of the
spectrum. Evaluate the focal field with [`PlaneWaveSummation`](@ref).

# Model

Ray `j` of the bundle has position `p`, unit direction `d`, optical path length `opl`,
power `P` and unit polarization `e` (see [`PolarizedRayBundle`](@ref)). It becomes the
sample

```math
\\mathbf{s}_j = \\mathbf{d}_j, \\qquad
\\boldsymbol{\\mathcal{E}}_j = -\\frac{i}{\\lambda_m}\\sqrt{\\frac{P_j}{\\kappa\\, w_j}}\\;
\\mathbf{e}_j\\, \\exp\\!\\big(i k_0\\, \\mathrm{opl}_j - i k\\, \\mathbf{d}_j\\cdot(\\mathbf{p}_j - \\mathbf{r}_0)\\big),
```

with `r₀ = origin(port)`, `λₘ = λ/n` and `k = 2π/λₘ` in the medium `n` of the port,
`k₀ = 2π/λ`, κ = [`power_normalization`](@ref)`(port)` and `w_j` the solid angle of the
ray in \\[sr\\]. The phase is the ray phase continued along the ray to the plane through
`r₀` perpendicular to `d` (absolute OPL convention); the factor `−i` is the phase anomaly
of a converging wave in the focus (Born & Wolf, 8.8, time dependence exp(−iωt)). The
polarization `e_j` stays transverse to `d_j`, so vectorial effects of high numerical
apertures (longitudinal fields) are included.

# Solid angles

A `PolarizedRayBundle` stores the power of each ray, not the size of its ray tube, so the
solid angle per ray must come from the sampling:

  - `solid_angle = nothing` (default): Voronoi cells of the ray directions in the plane of
    direction cosines `(d·u, d·v)` of `port`, clipped to their convex hull, divided by
    `d·n` (the Jacobian of `dΩ = ds_u ds_v / s_n`). Needs DelaunayTriangulation.jl
    (`using DelaunayTriangulation` loads the extension `OpticsBaseDelaunayTriangulationExt`).
  - `solid_angle`: a vector of `M` positive solid angles in \\[sr\\], one per ray in bundle
    order, for samplings whose weights are known (e.g. quadrature rules).

# Approximations and conservation

  - Debye approximation: the Fresnel number of the converging wave is large and the field
    is evaluated near the focus. The discrete spectrum reproduces the focal field within
    about `λₘ / δs` of the focus, with `δs` the spacing of neighboring ray directions.
  - Converging bundles only: every ray must propagate towards the focus,
    `(r₀ − p_j)·d_j > 0`, and into the port, `d_j·n > 0`.
  - Power is conserved exactly: `total_power(spectrum) == total_power(bundle)` for any
    solid angles. The solid angles only shape the focal field; Voronoi cells of the hull
    rays are cut at the hull, so the aperture edge is sampled to first order (the error of
    the focal field decreases like `1/√M`).

# Fields

  - `port`: the port of the spectrum; its origin is the focus, its normal bounds the
    hemisphere of directions, and its medium must equal that of the bundle's port
  - `solid_angle`: `nothing` or the per-ray solid angles in \\[sr\\]
"""
struct DebyeWolf{P <: AbstractPort, W <: Union{Nothing, AbstractVector{<:Real}}} <:
       AbstractFieldConverter
    port::P
    solid_angle::W

    function DebyeWolf(port::P, solid_angle::W) where {P <: AbstractPort,
            W <: Union{Nothing, AbstractVector{<:Real}}}
        solid_angle === nothing || all(w -> isfinite(w) && w > 0, solid_angle) ||
            throw(ArgumentError("DebyeWolf: all solid angles must be finite and > 0"))
        return new{P, W}(port, solid_angle)
    end
end

DebyeWolf(port::AbstractPort; solid_angle = nothing) = DebyeWolf(port, solid_angle)

input_representation(::DebyeWolf) = PolarizedRayBundle
output_representation(::DebyeWolf) = PlaneWaveSpectrum{3}

function __convert_field(conv::DebyeWolf, bundle::PolarizedRayBundle{T}) where {T}
    p = conv.port
    tol = sqrt(eps(T))
    n_med = T(refractive_index(p))
    isapprox(n_med, T(refractive_index(port(bundle))); rtol = tol) ||
        throw(ArgumentError("DebyeWolf: the port must be in the medium of the bundle, " *
                            "got refractive indices $(refractive_index(p)) and " *
                            "$(refractive_index(port(bundle)))"))
    r0 = SVector{3, T}(origin(p))
    A = local_axes(p)
    u = SVector{3, T}(A[:, 1])
    v = SVector{3, T}(A[:, 2])
    n = SVector{3, T}(A[:, 3])
    M = length(bundle)
    for j in 1:M
        d = bundle.direction[j]
        dot(r0 - bundle.position[j], d) > 0 ||
            throw(ArgumentError("DebyeWolf: ray $j does not converge towards the port " *
                                "origin (focus): (origin − position)·direction ≤ 0"))
        dot(d, n) > 0 ||
            throw(ArgumentError("DebyeWolf: ray $j does not propagate into the port " *
                                "(direction·normal = $(dot(d, n)) ≤ 0)"))
    end
    w = _solid_angles(conv.solid_angle, bundle.direction, u, v, n, T)

    λ = wavelength(bundle)
    λₘ = λ / n_med
    k = 2 * T(π) / λₘ
    k0 = 2 * T(π) / λ
    κ = T(power_normalization(p))
    amplitude = map(1:M) do j
        d = bundle.direction[j]
        phase = k0 * bundle.opl[j] - k * dot(d, bundle.position[j] - r0)
        a = -im / λₘ * sqrt(bundle.power[j] / (κ * w[j])) * cis(phase)
        a * bundle.polarization[j]
    end
    return PlaneWaveSpectrum(p, λ, copy(bundle.direction), amplitude, w)
end

function _solid_angles(w::AbstractVector, directions, _, _, _, ::Type{T}) where {T}
    length(w) == length(directions) ||
        throw(ArgumentError("DebyeWolf: solid_angle must have one entry per ray, got " *
                            "$(length(w)) for $(length(directions)) rays"))
    return Vector{T}(w)
end

function _solid_angles(::Nothing, directions, u, v, n, ::Type{T}) where {T}
    x = T[dot(d, u) for d in directions]
    y = T[dot(d, v) for d in directions]
    _check_triangulable(x, y)
    areas = _delaunay_dual_areas(x, y)
    return T[areas[j] / dot(directions[j], n) for j in eachindex(directions)]
end

# At least three distinct, not collinear points, so the Voronoi cells clipped to the
# convex hull have positive areas.
function _check_triangulable(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    M = length(x)
    M ≥ 3 ||
        throw(ArgumentError("DebyeWolf: computing solid angles needs at least 3 rays, " *
                            "got $M; pass `solid_angle` explicitly"))
    allunique(zip(x, y)) ||
        throw(ArgumentError("DebyeWolf: two rays have the same direction, so their solid " *
                            "angles are undefined; pass `solid_angle` explicitly"))
    # Farthest point from the first one, then the largest triangle area with both
    j = argmax(i -> (x[i] - x[1])^2 + (y[i] - y[1])^2, 1:M)
    ex, ey = x[j] - x[1], y[j] - y[1]
    area = maximum(i -> abs(ex * (y[i] - y[1]) - ey * (x[i] - x[1])), 1:M)
    area > sqrt(eps(T)) * (ex^2 + ey^2) ||
        throw(ArgumentError("DebyeWolf: the ray directions are collinear, so they span no " *
                            "solid angle; pass `solid_angle` explicitly"))
    return nothing
end

# Backend of the default solid angles, implemented by the extension
# `OpticsBaseDelaunayTriangulationExt` as
#     _delaunay_dual_areas(x::Vector{T}, y::Vector{T}) -> Vector{T}
# It returns the areas of the Voronoi cells of the points (x[j], y[j]) clipped to their
# convex hull (the cells partition the hull). The points are distinct and not all
# collinear. OpticsBase itself defines no method. (Internal, hence no docstring.)
function _delaunay_dual_areas end

function _register_delaunay_hint()
    Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
        if exc.f === _delaunay_dual_areas &&
           Base.get_extension(@__MODULE__, :OpticsBaseDelaunayTriangulationExt) === nothing
            print(io, "\nDebyeWolf computes solid angles with DelaunayTriangulation.jl: run ",
                "`using DelaunayTriangulation` to load the extension ",
                "OpticsBaseDelaunayTriangulationExt, or pass `solid_angle` explicitly.")
        end
    end
    return nothing
end
