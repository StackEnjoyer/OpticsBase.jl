# Centered FFT backend of `PlaneWaveDecomposition`. The padding, scaling and the mapping
# onto plane waves live in OpticsBase (src/PlaneWaveDecomposition.jl).
module OpticsBaseFFTWExt

using OpticsBase: OpticsBase
using FFTW: fft, fftshift, ifftshift

# Zero coordinate and zero frequency both at index M÷2 + 1 per dimension (the
# `RegularGrid` centering); forward kernel exp(−i k·x), unnormalized.
function OpticsBase._centered_fft(E::AbstractMatrix{<:Complex})
    return Matrix(fftshift(fft(ifftshift(E))))
end

end # module OpticsBaseFFTWExt
