```@meta
CurrentModule = OpticsBase
```

# Converters

A converter turns a field of one representation into another, e.g. a ray bundle into a
sampled field. Conversions are almost always approximations, so converters are explicit
objects that carry their parameters, and each converter states what it assumes and what it
conserves. `OpticsBase` never converts implicitly.

```julia
conv = GaussianBeamletSummation(RegularGrid((256, 256), (10e-6, 10e-6)))
field = convert_field(conv, bundle)   # RayBundle → SampledField at the same port
```

`convert_field` checks that the converter accepts the field (otherwise a
[`MissingConverterError`](@ref) is thrown) and that the result has the promised
representation, just like `solve` does for propagation algorithms.

## Converter interface

```@docs
AbstractFieldConverter
convert_field
__convert_field
```

## Beamlet summation

```@docs
GaussianBeamletSummation
```
