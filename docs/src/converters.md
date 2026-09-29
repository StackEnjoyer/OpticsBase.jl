```@meta
CurrentModule = OpticsBase
```

# Converters

A converter turns optical data of one representation into another, e.g. a plane-wave spectrum into a
sampled field. Conversions are almost always approximations, so converters are explicit
objects that carry their parameters, and each converter states what it assumes and what it
conserves. `OpticsBase` never converts implicitly. Summing the Gaussian beamlets of a beamlet
tracer is not a converter here: it is the tracer's job, see [Solvers and extensions](@ref).

```julia
conv = PlaneWaveSummation(RegularGrid((256, 256), (10e-6, 10e-6)))
field = convert_field(conv, spectrum)   # PlaneWaveSpectrum → SampledField
```

`convert_field` checks that the converter accepts the data, a field or a ray bundle (otherwise a
[`MissingConverterError`](@ref) is thrown) and that the result has the promised
representation, just like `solve` does for propagation algorithms.

## Converter interface

```@docs
AbstractFieldConverter
convert_field
__convert_field
```


## Plane waves

A [`PlaneWaveSpectrum`](@ref) is the node between rays and waves. Converging rays become
plane waves with [`DebyeWolf`](@ref) (input: a [`PolarizedRayBundle`](@ref)), a sampled field is decomposed with
[`PlaneWaveDecomposition`](@ref), and [`PlaneWaveSummation`](@ref) evaluates any spectrum
on a planar or volume grid at any port in the same medium, e.g. in a focus:

```julia
using OpticsBase, DelaunayTriangulation   # solid angles per ray (extension)

spectrum = convert_field(DebyeWolf(port_focus), bundle)        # PolarizedRayBundle → spectrum
focus = convert_field(PlaneWaveSummation(grid), spectrum)      # → SampledField at port_focus

using FFTW                                                     # FFT backend (extension)
spectrum = convert_field(PlaneWaveDecomposition(), field)      # SampledField → spectrum
```

```@docs
DebyeWolf
PlaneWaveDecomposition
PlaneWaveSummation
```
