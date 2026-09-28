# Numerics of `AngularSpectrumMethod`: the array call into WaveOpticsPropagation.jl. The
# port-geometry checks and the per-component loop live in OpticsBase (src/FreeSpace.jl).
module OpticsBaseWaveOpticsPropagationExt

using OpticsBase: OpticsBase, AngularSpectrumMethod
using WaveOpticsPropagation: WaveOpticsPropagation

# WOP's kernel is exp(i k z √(1 − (λf)²)) with k = 2π/λ, consistent with exp(−iωt), so the
# medium wavelength λₘ is passed as λ. WOP names the dimensions "y, x"; these are labels
# only: dim 1 uses L[1], dim 2 uses L[2], and both are centered at sample N÷2 + 1 like
# `RegularGrid`. `AngularSpectrum` (not `angular_spectrum`) also accepts non-square fields.
function OpticsBase._angular_spectrum(E::AbstractMatrix{<:Complex}, z::Real, λₘ::Real,
        L::NTuple{2, Real}, alg::AngularSpectrumMethod)
    as = WaveOpticsPropagation.AngularSpectrum(E, z, λₘ, L; padding = alg.padding,
        pad_factor = alg.pad_factor, bandlimit = alg.bandlimit)
    # The result is a view into the propagator's buffer: return an independent copy.
    return copy(as(E))
end

end # module OpticsBaseWaveOpticsPropagationExt
