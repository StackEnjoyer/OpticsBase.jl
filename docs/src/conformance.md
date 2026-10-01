```@meta
CurrentModule = OpticsBase.Conformance
```

# Conformance tests

```@docs
Conformance
```

## Example: BeamletOptics as a source

```julia
using BeamletOptics, OpticsBase, Test
using OpticsBase: Conformance

function make(beam, origin, axes, size, spacing)
    # a BeamletOptics beamlet with the waist, direction, power and polarization of `beam`
    u, v = axes[:, 1], axes[:, 2]
    E0 = sqrt(4 * VACUUM_IMPEDANCE * beam.P / (π * beam.w0^2)) .*
         (beam.jones[1] .* u .+ beam.jones[2] .* v)
    b = AstigmaticGaussianBeamlet(beam.waist, beam.direction, beam.λ, beam.w0, beam.w0; E0)
    detector = Detector(0.05)          # placed at `origin`, facing the beam
    # ... translate and align the detector, solve_system! ...
    return PlaneField(detector; size, spacing, origin, axes)
end

for r in Conformance.check_source(make)
    @test r.value <= r.limit
end
```

## Reference

```@docs
GaussianBeam
field
check_source
check_propagator
compare
Result
```
