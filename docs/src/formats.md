```@meta
CurrentModule = OpticsBase
```

# Exchange formats

An exchange format is a representation of an optical field that one solver hands to the
next at a port. Every format is a subtype of [`AbstractOpticalField`](@ref) and describes
one monochromatic, coherent component; the physical conventions are on the
[Conventions](@ref) page.

## Field interface

```@docs
AbstractOpticalField
port
wavelength
total_power
is_vectorial
is_coherent
```

## Ray bundles

```@docs
RayBundle
has_beamlets
```

## Sampled fields

```@docs
SampledField
field_array
grid
```

## Sampling grids

```@docs
AbstractGrid
RegularGrid
Base.size(::RegularGrid)
spacing
coordinates
```
