```@meta
CurrentModule = OpticsBase
```

# Exchange formats

An exchange format is a representation of optical data that one solver hands to the next
at a port. Formats form one hierarchy with two branches under a common root
[`AbstractOpticalData`](@ref) (interface: [`port`](@ref), [`wavelength`](@ref)):

```text
AbstractOpticalData
├── AbstractOpticalField{N}      Maxwell fields, N = 1 scalar, N = 3 vectorial
│   ├── SampledField{N, …}
│   └── PlaneWaveSpectrum{N, …}
└── AbstractRayBundle            geometrical optics
    ├── RayBundle
    └── PolarizedRayBundle
```

Every format describes one monochromatic, coherent component; the physical conventions
are on the [Conventions](@ref) page. Fields are handed to field solvers, rays to other ray
tracers or into a ray-to-field converter such as [`DebyeWolf`](@ref). `OpticsBase` knows no
solver: beamlets, detectors and other solver internals never appear in an exchange format.

## Common interface

```@docs
AbstractOpticalData
port
wavelength
```

## Fields

Fields are Maxwell fields on a screen or on the sphere of directions. The parameter `N`
is the number of field components: `N = 1` scalar, `N = 3` vectorial.

```@docs
AbstractOpticalField
AbstractVectorField
AbstractScalarField
total_power
is_vectorial
is_coherent
```

## Ray bundles

Rays are the geometrical-optics limit. [`RayBundle`](@ref) carries only geometry and is
what ray tracers exchange; [`PolarizedRayBundle`](@ref) additionally carries optical path
length, power and polarization per ray, which a ray-to-field converter needs.

```@docs
AbstractRayBundle
positions
directions
RayBundle
PolarizedRayBundle
```

## Plane-wave spectra

```@docs
PlaneWaveSpectrum
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
