"""
    GaussianBeamletSummation(grid::AbstractGrid{2}) <: AbstractFieldConverter

Converter from a [`RayBundle`](@ref) of Gaussian beamlets to a [`SampledField`](@ref):
the beamlet fields are summed coherently on `grid` in the plane of the bundle's port. The
result is given at the same port, with the same wavelength and number of components `N`.

# Model

Every ray `j` carries a fundamental Gaussian beamlet (see [`RayBundle`](@ref)): complex
matrix `Q` in its [`ray_basis`](@ref), power `P`, phasor `e` and optical path length
`opl`. At a point `r` of the port plane, with `s = (r − p)·d` along the ray and transverse
coordinates `x = ((r − p)·x̂, (r − p)·ŷ)`,

```math
\\mathbf{E}_j(\\mathbf{r}) = a_j\\,\\mathbf{e}_j\\, e^{i k_0\\,\\mathrm{opl}_j}\\,
\\frac{\\exp\\!\\big(i k s + \\tfrac{i k}{2}\\, x^\\mathsf{T} Q(s)\\, x\\big)}{\\sqrt{\\det(I + sQ)}},
\\qquad Q(s) = Q\\,(I + sQ)^{-1},
```

with `k = 2πn/λ` in the port medium, `k₀ = 2π/λ` and the peak amplitude
`a_j = √(P_j k √det(Im Q) / (κ π))`, κ = [`power_normalization`](@ref)`(port)`. The field
is the sum over all rays.

# Approximations

  - Paraxial beamlets: each beamlet is the fundamental Gaussian mode with complex
    curvature `Q`, propagated in the port medium along its own direction to reach points
    of the port plane off its transverse plane (`s ≠ 0`, oblique rays).
  - The polarization (`e_j`) is constant across a beamlet.
  - Each beamlet carries exactly its power `P_j` through its own transverse plane. Through
    the port plane, the sampled power of an oblique beamlet is `P_j / cos θ`, because the
    power normalization has no obliquity factor (see the conventions).
  - Overlapping beamlets interfere, so `total_power` of the result is in general not the
    sum of the beamlet powers.
  - The sampled field is exact at the grid points; `total_power` of the result is a
    rectangle rule and needs a grid spacing well below the beamlet widths and a grid that
    covers the beamlets.

Pure rays (`has_beamlets(bundle) == false`) are not accepted: `convert_field` throws a
[`MissingConverterError`](@ref).

# Fields

  - `grid`: 2D sampling grid in the port-local coordinates `(ξ, η)`, see
    [`RegularGrid`](@ref)
"""
struct GaussianBeamletSummation{G <: AbstractGrid{2}} <: AbstractFieldConverter
    grid::G
end

input_representation(::GaussianBeamletSummation) = RayBundle{<:Any, <:Any, <:Any, <:AbstractVector}
output_representation(::GaussianBeamletSummation) = SampledField

function __convert_field(conv::GaussianBeamletSummation, bundle::RayBundle{N}) where {N}
    p = port(bundle)
    grid = conv.grid
    ξs = coordinates(grid, 1)
    ηs = coordinates(grid, 2)
    T = float(promote_type(typeof(wavelength(bundle)), eltype(ξs), eltype(ηs)))
    λ = T(wavelength(bundle))
    k0 = 2 * T(π) / λ
    k = k0 * T(refractive_index(p))
    κ = T(power_normalization(p))
    A = local_axes(p)
    u = SVector{3, T}(A[:, 1])
    v = SVector{3, T}(A[:, 2])
    o = SVector{3, T}(origin(p))

    beamlets = [_beamlet_data(p, bundle, j, T, k0, k, κ, o, u, v) for j in 1:length(bundle)]

    E = zeros(Complex{T}, length(ξs), length(ηs), N)
    Threads.@threads for iη in eachindex(ηs)
        η = ηs[iη]
        for b in beamlets, (iξ, ξ) in enumerate(ξs)
            val = _beamlet_value(b, T(ξ), T(η), k)
            for c in 1:N
                E[iξ, iη, c] += val * b.amplitude[c]
            end
        end
    end
    return SampledField(E, grid, p, λ)
end

# Per-ray data for the summation. With r = o + ξu + ηv on the port plane, s and the
# transverse coordinates are affine in (ξ, η): s = s0 + ξ su + η sv, likewise x1 and x2.
function _beamlet_data(p, bundle::RayBundle{N}, j, ::Type{T}, k0, k, κ, o, u, v) where {N, T}
    pos = SVector{3, T}(bundle.position[j])
    d = SVector{3, T}(bundle.direction[j])
    Q = SMatrix{2, 2, Complex{T}, 4}(bundle.beamlet[j])
    x̂, ŷ = ray_basis(p, d)
    x̂ = SVector{3, T}(x̂)
    ŷ = SVector{3, T}(ŷ)
    Δ0 = o - pos
    a = sqrt(T(bundle.power[j]) * k * sqrt(det(imag(Q))) / (κ * T(π)))
    amplitude = a * cis(k0 * T(bundle.opl[j])) * SVector{N, Complex{T}}(bundle.phasor[j])
    # Eigenvalues of the complex symmetric Q; Im μ > 0 because Im Q is positive definite,
    # so the principal roots of 1 + sμ are continuous in s.
    h = tr(Q) / 2
    disc = sqrt(h^2 - det(Q))
    return (; Q, amplitude, μ1 = h + disc, μ2 = h - disc,
        s = SVector(dot(Δ0, d), dot(u, d), dot(v, d)),
        x1 = SVector(dot(Δ0, x̂), dot(u, x̂), dot(v, x̂)),
        x2 = SVector(dot(Δ0, ŷ), dot(u, ŷ), dot(v, ŷ)))
end

# Scalar beamlet factor at the port point (ξ, η), without amplitude and polarization.
function _beamlet_value(b, ξ, η, k)
    s = b.s[1] + ξ * b.s[2] + η * b.s[3]
    x = SVector(b.x1[1] + ξ * b.x1[2] + η * b.x1[3], b.x2[1] + ξ * b.x2[2] + η * b.x2[3])
    D = I + s * b.Q
    Qs = b.Q * inv(D)
    quad = transpose(x) * Qs * x
    return exp(im * k * (s + quad / 2)) / (sqrt(1 + s * b.μ1) * sqrt(1 + s * b.μ2))
end
