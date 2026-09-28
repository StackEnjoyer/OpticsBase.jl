```@meta
CurrentModule = OpticsBase
```

# Ports

A port is the surface at which a field is handed from one solver to the next. It fixes the
local frame, the polarization basis and the medium. The binding rules are on the
[Conventions](@ref) page.

## Port types

```@docs
AbstractPort
PlanarPort
origin
local_axes
normal
refractive_index
```

## Coordinates

```@docs
to_local
to_global
```

## Polarization basis

```@docs
ray_basis
jones_to_global
circular_jones
```

## Power normalization factor

```@docs
power_normalization
VACUUM_IMPEDANCE
```
