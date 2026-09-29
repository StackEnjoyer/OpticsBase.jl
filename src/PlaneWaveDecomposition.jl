"""
    PlaneWaveDecomposition(; pad_factor = 1) <: AbstractFieldConverter

Converter from a planar [`SampledField`](@ref) to a [`PlaneWaveSpectrum`](@ref): the field
in the port plane is decomposed into homogeneous plane waves by a discrete Fourier
transform. The result is given at the same port (phase reference r₀ = `origin(port)`,
the grid center), with the same wavelength and number of components `N`. The FFT is done
by [FFTW.jl](https://github.com/JuliaMath/FFTW.jl), which has to be loaded (`using FFTW`)
to activate the extension `OpticsBaseFFTWExt`; without it, `convert_field` throws a
`MethodError`.

# Model

With the port-local grid coordinates `x = (ξ, η)` along `(u, v)`, spacings `(Δξ, Δη)`, the
wavelength in the medium `λₘ = λ/refractive_index(port)` and `k = 2π/λₘ`, the transverse
field `E(x)` (zero-padded to `(M_ξ, M_η) = pad_factor .* size(grid)` samples, centered as
in [`RegularGrid`](@ref)) is transformed as

```math
F(\\mathbf{k}_t) = \\frac{\\Delta\\xi\\,\\Delta\\eta}{(2\\pi)^2}
\\sum_{\\mathbf{x}} E(\\mathbf{x})\\, e^{-i \\mathbf{k}_t\\cdot\\mathbf{x}},
\\qquad k_\\xi = m\\,\\frac{2\\pi}{M_\\xi\\,\\Delta\\xi},\\quad
k_\\eta = m'\\,\\frac{2\\pi}{M_\\eta\\,\\Delta\\eta},
```

so that `E(x) = Σ F(k_t) e^{i k_t·x} Δk_ξ Δk_η` (sign consistent with exp(−iωt)). Every
wave vector with `k_t² < k²` becomes one sample with `k_n = √(k² − k_t²)`:

  - direction `s = (k_ξ u + k_η v + k_n n)/k` (global frame, unit),
  - solid angle `w = Δk_ξ Δk_η/(k k_n)` in \\[sr\\],
  - spectral density `ℰ = k k_n F` in \\[V/m/sr\\] (`N = 1`), or for `N = 3`
    `ℰ = k k_n (F_u u + F_v v + F_n n)` with `F_n = −(k_ξ F_u + k_η F_v)/k_n` (global
    frame), such that `w ℰ = F Δk_ξ Δk_η`.

The samples are stored in column-major `(ξ, η)` order of the wave-vector grid.
[`PlaneWaveSummation`](@ref) at the same grid and port reproduces the input field.

# Approximations

  - The field is assumed periodic on the (padded) grid of extent
    `M .* (Δξ, Δη)`; it must decay towards the grid edges. Zero padding refines the
    wave-vector sampling but does not add information.
  - Evanescent components (`k_t² ≥ k²`) are dropped together with their power.
  - `N = 3`: only the transverse components `E_u = u·E`, `E_v = v·E` at the port are used;
    the normal component of each plane wave is recomputed from transversality
    (`k·E = 0`). The input `E_n` is ignored.
  - [`total_power`](@ref) of the spectrum is the exact flux through the port plane. It
    differs from the paraxial `total_power` of the [`SampledField`](@ref) by O(θ²) for
    plane waves at angle θ to `n` (no obliquity weighting in `SampledField`).

Only planar fields (`SampledField` with a 2D grid) at a [`PlanarPort`](@ref) are accepted;
3D fields are rejected by `convert_field` with a [`MissingConverterError`](@ref), other
ports with an `ArgumentError`.

# Arguments

  - `pad_factor`: padded size per dimension as a multiple of the grid size, integer
    `≥ 1` (`ArgumentError` otherwise)

# Fields

  - `pad_factor::Int`: zero-padding factor per dimension
"""
struct PlaneWaveDecomposition <: AbstractFieldConverter
    pad_factor::Int

    function PlaneWaveDecomposition(pad_factor::Integer)
        pad_factor ≥ 1 ||
            throw(ArgumentError("PlaneWaveDecomposition: pad_factor must be ≥ 1, got $pad_factor"))
        return new(Int(pad_factor))
    end
end

PlaneWaveDecomposition(; pad_factor::Integer = 1) = PlaneWaveDecomposition(pad_factor)

input_representation(::PlaneWaveDecomposition) = SampledField{<:Any, 2}
output_representation(::PlaneWaveDecomposition) = PlaneWaveSpectrum

function __convert_field(conv::PlaneWaveDecomposition, field::SampledField{N, 2, T}) where {N, T}
    p = _decomposition_port(port(field))
    λ = wavelength(field)
    k = 2 * T(π) * T(refractive_index(p)) / λ
    ax = local_axes(p)
    u = SVector{3, T}(ax[:, 1])
    v = SVector{3, T}(ax[:, 2])
    n = SVector{3, T}(ax[:, 3])
    g = grid(field)
    Δξ, Δη = T.(spacing(g))
    M = conv.pad_factor .* size(g)
    scale = Δξ * Δη / (2 * T(π))^2
    Fs = map(E -> _centered_fft(_pad_centered(E, M)) .* scale,
        _transverse_components(Val(N), field_array(field), u, v))
    Δkξ = 2 * T(π) / (M[1] * Δξ)
    Δkη = 2 * T(π) / (M[2] * Δη)
    cξ = M[1] ÷ 2 + 1
    cη = M[2] ÷ 2 + 1

    dirs = SVector{3, T}[]
    amps = SVector{N, Complex{T}}[]
    ws = T[]
    for iη in 1:M[2], iξ in 1:M[1]
        kξ = (iξ - cξ) * Δkξ
        kη = (iη - cη) * Δkη
        kt2 = kξ^2 + kη^2
        kt2 < k^2 || continue
        kn = sqrt(k^2 - kt2)
        push!(dirs, (kξ * u + kη * v + kn * n) / k)
        push!(ws, Δkξ * Δkη / (k * kn))
        push!(amps, _spectral_density(Fs, iξ, iη, k, kξ, kη, kn, u, v, n))
    end
    return PlaneWaveSpectrum(p, λ, dirs, amps, ws)
end

_decomposition_port(p::PlanarPort) = p
function _decomposition_port(p::AbstractPort)
    throw(ArgumentError("PlaneWaveDecomposition: the field must be given at a PlanarPort, got a $(typeof(p))"))
end

# Transverse field components at the port: the scalar field (N = 1) or (E_u, E_v) (N = 3).
_transverse_components(::Val{1}, E, _, _) = (E[:, :, 1],)
function _transverse_components(::Val{3}, E, u, v)
    Eu = u[1] .* view(E, :, :, 1) .+ u[2] .* view(E, :, :, 2) .+ u[3] .* view(E, :, :, 3)
    Ev = v[1] .* view(E, :, :, 1) .+ v[2] .* view(E, :, :, 2) .+ v[3] .* view(E, :, :, 3)
    return (Eu, Ev)
end

# Zero-pads `E` to size `M` so that input index n÷2 + 1 lands on M÷2 + 1 per dimension.
function _pad_centered(E::AbstractMatrix, M::NTuple{2, Int})
    size(E) == M && return Matrix(E)
    out = zeros(eltype(E), M)
    o1 = M[1] ÷ 2 - size(E, 1) ÷ 2
    o2 = M[2] ÷ 2 - size(E, 2) ÷ 2
    out[(o1 + 1):(o1 + size(E, 1)), (o2 + 1):(o2 + size(E, 2))] .= E
    return out
end

# Spectral density ℰ = k k_n F of one sample; for N = 3 with F_n from transversality.
function _spectral_density(Fs::NTuple{1}, iξ, iη, k, _, _, kn, _, _, _)
    return SVector(k * kn * Fs[1][iξ, iη])
end
function _spectral_density(Fs::NTuple{2}, iξ, iη, k, kξ, kη, kn, u, v, n)
    Fu = Fs[1][iξ, iη]
    Fv = Fs[2][iξ, iη]
    return k * (kn * Fu * u + kn * Fv * v - (kξ * Fu + kη * Fv) * n)
end

# Backend of `PlaneWaveDecomposition`, implemented by the extension `OpticsBaseFFTWExt` as
#     _centered_fft(E::AbstractMatrix{<:Complex}) -> Matrix
# the unnormalized forward DFT (kernel exp(−2πi (m − c)(l − c)/M), c = M÷2 + 1 per
# dimension) with the zero frequency and the zero coordinate both at index M÷2 + 1, i.e.
# fftshift(fft(ifftshift(E))). It returns a new array. OpticsBase itself defines no method.
# (Internal, hence no docstring: the docs build checks that every docstring is embedded.)
function _centered_fft end

# Called from `__init__`: hint on the `MethodError` of the missing FFT backend.
function _register_fftw_hint()
    Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
        if exc.f === _centered_fft &&
           Base.get_extension(@__MODULE__, :OpticsBaseFFTWExt) === nothing
            print(io, "\nPlaneWaveDecomposition needs FFTW.jl for the FFT: run ",
                "`using FFTW` to load the extension OpticsBaseFFTWExt.")
        end
    end
    return nothing
end
