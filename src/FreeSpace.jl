"""
    FreeSpace()

Optical system for [`PropagationProblem`](@ref): a homogeneous, isotropic, lossless medium
between the input port (the port of the field) and the output port `port_out`. It has no
parameters: the medium is the one of the ports (their [`refractive_index`](@ref)), and the
propagation distance and direction follow from the port geometry.

Solved by [`AngularSpectrumMethod`](@ref), which supports coaxial planar ports only.
"""
struct FreeSpace end

"""
    AngularSpectrumMethod(; padding = true, pad_factor = 2, bandlimit = true)

Propagation algorithm for a [`PropagationProblem`](@ref) with a planar
[`SampledField`](@ref) and the system [`FreeSpace`](@ref): the scalar angular-spectrum
method of plane waves, applied to each field component separately. The numerics are done
by [WaveOpticsPropagation.jl](https://github.com/JuliaPhysics/WaveOpticsPropagation.jl),
which has to be loaded (`using WaveOpticsPropagation`) to activate the extension
`OpticsBaseWaveOpticsPropagationExt`; without it, `solve` throws a `MethodError`.

# Arguments

  - `padding`: zero-pad the field before the FFT to suppress the wrap-around of the
    cyclic convolution
  - `pad_factor`: padded size per dimension as a multiple of the grid size, integer `≥ 1`;
    only used if `padding = true`
  - `bandlimit`: apply the band limit of Matsushima & Shimobaba (Opt. Express 17, 19662
    (2009)) with a smooth Hann edge, which suppresses aliasing of the undersampled
    transfer function for long distances

# Geometry

Only coaxial propagation is supported. With `port_in = port(prob.field)` and
`port_out = prob.port_out`:

  - both ports are [`PlanarPort`](@ref)s,
  - `normal(port_out) == normal(port_in)` and the first axes `u` are equal (within
    `√eps` of the port element type),
  - `origin(port_out) − origin(port_in) = z·n` with `n = normal(port_in)` and `z` of
    either sign; the lateral part must be `≤ √eps·max(‖Δ‖, λ)`,
  - both ports have the same refractive index `n_med` (within `√eps` relative).

Anything else throws an `ArgumentError` before any numerics run. The output field lives on
the same grid as the input field, at `port_out`.

# Conventions and approximations

  - Transfer function `exp(i k z √(1 − λₘ²(f_ξ² + f_η²)))` with `λₘ = λ/n_med` and
    `k = 2π/λₘ`, consistent with the time convention exp(−iωt): propagation by `z > 0`
    adds the phase `+k z` to a plane wave along `n`. Evanescent components decay for
    `z > 0`; for `z < 0` the method uses the complex conjugate kernel (back-propagation
    of the propagating part).
  - The field is periodic on the (padded) grid of extent `size(grid) .* spacing(grid)`;
    it must decay towards the grid edges. The grid follows the [`RegularGrid`](@ref)
    convention (port origin at sample `n÷2 + 1`), which is the centering WOP uses.
  - Each Cartesian component (`N = 3`: global `Ex, Ey, Ez`; `N = 1`: scalar) is propagated
    with the same scalar kernel. This is exact for homogeneous media: every plane-wave
    component keeps its field vector.
  - Power is conserved up to the band limit, the padding window and evanescent parts;
    κ = [`power_normalization`](@ref) is the same at both ports since `n_med` is equal.

# Interface

`input_representation` and `output_representation` are `SampledField{<:Any, 2}` (planar
fields, scalar or vectorial).
"""
struct AngularSpectrumMethod <: AbstractPropagationAlgorithm
    padding::Bool
    pad_factor::Int
    bandlimit::Bool

    function AngularSpectrumMethod(padding::Bool, pad_factor::Integer, bandlimit::Bool)
        pad_factor ≥ 1 ||
            throw(ArgumentError("AngularSpectrumMethod: pad_factor must be ≥ 1, got $pad_factor"))
        return new(padding, Int(pad_factor), bandlimit)
    end
end

function AngularSpectrumMethod(; padding::Bool = true, pad_factor::Integer = 2,
        bandlimit::Bool = true)
    return AngularSpectrumMethod(padding, pad_factor, bandlimit)
end

input_representation(::AngularSpectrumMethod) = SampledField{<:Any, 2}
output_representation(::AngularSpectrumMethod) = SampledField{<:Any, 2}

function __solve(prob::PropagationProblem{<:SampledField, FreeSpace},
        alg::AngularSpectrumMethod)
    field = prob.field
    λ = wavelength(field)
    z, n_med = _coaxial_geometry(port(field), prob.port_out, λ)
    λₘ = λ / n_med
    g = grid(field)
    L = size(g) .* spacing(g)
    E_out = _propagate_components(field_array(field), z, λₘ, L, alg)
    return PropagationSolution(SampledField(E_out, g, prob.port_out, λ), prob, alg)
end

# Per component: copy the slice (the backend needs a dense matrix) and propagate it.
function _propagate_components(E::AbstractArray{<:Complex, 3}, z, λₘ, L, alg)
    E_out = similar(E)
    for c in axes(E, 3)
        E_out[:, :, c] .= _angular_spectrum(E[:, :, c], z, λₘ, L, alg)
    end
    return E_out
end

# Checks the D5 geometry (coaxial planar ports, same medium) and returns the signed
# distance z along the input normal and the refractive index of the medium.
function _coaxial_geometry(port_in::PlanarPort, port_out::PlanarPort, λ::Real)
    T = float(promote_type(eltype(origin(port_in)), eltype(origin(port_out)), typeof(λ)))
    tol = sqrt(eps(T))
    n = normal(port_in)
    norm(normal(port_out) - n) ≤ tol ||
        throw(ArgumentError("AngularSpectrumMethod: the output port must be parallel to " *
                            "the input port (same normal), got normals $n and $(normal(port_out))"))
    u_in = local_axes(port_in)[:, 1]
    u_out = local_axes(port_out)[:, 1]
    norm(u_out - u_in) ≤ tol ||
        throw(ArgumentError("AngularSpectrumMethod: the output port must have the same " *
                            "axis u as the input port, got $u_in and $u_out"))
    Δ = origin(port_out) - origin(port_in)
    z = dot(Δ, n)
    lateral = norm(Δ - z * n)
    lateral ≤ tol * max(norm(Δ), λ) ||
        throw(ArgumentError("AngularSpectrumMethod: the output port origin must lie on " *
                            "the normal through the input port origin (coaxial " *
                            "propagation only), got a lateral offset of $lateral m"))
    n_in = refractive_index(port_in)
    n_out = refractive_index(port_out)
    isapprox(n_in, n_out; rtol = tol) ||
        throw(ArgumentError("AngularSpectrumMethod: both ports must be in the same " *
                            "medium, got refractive indices $n_in and $n_out"))
    return z, n_in
end

function _coaxial_geometry(port_in::AbstractPort, port_out::AbstractPort, ::Real)
    throw(ArgumentError("AngularSpectrumMethod: both ports must be PlanarPorts, got " *
                        "$(typeof(port_in)) and $(typeof(port_out))"))
end

# Backend of `AngularSpectrumMethod`, implemented by the extension
# `OpticsBaseWaveOpticsPropagationExt` as
#     _angular_spectrum(E::AbstractMatrix{<:Complex}, z, λₘ, L::NTuple{2}, alg) -> matrix
# It propagates the scalar field `E` (one component; dim 1 along u, dim 2 along v) by the
# signed distance `z` [m] in a medium with wavelength `λₘ = λ/n` [m]; `L` is the grid
# extent `size .* spacing` [m]. It returns a new array that aliases neither `E` nor a
# backend buffer. OpticsBase itself defines no method. (Internal, hence no docstring:
# the docs build checks that every docstring is embedded.)
function _angular_spectrum end

function _register_error_hints()
    Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
        if exc.f === _angular_spectrum &&
           Base.get_extension(@__MODULE__, :OpticsBaseWaveOpticsPropagationExt) === nothing
            print(io, "\nAngularSpectrumMethod needs WaveOpticsPropagation.jl for the ",
                "numerics: run `using WaveOpticsPropagation` to load the extension ",
                "OpticsBaseWaveOpticsPropagationExt.")
        end
    end
    return nothing
end
