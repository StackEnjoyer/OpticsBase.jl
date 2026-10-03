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
    b = AstigmaticGaussianBeamlet(collect(beam.waist), collect(beam.direction), beam.λ,
        beam.w0, beam.w0; E0 = collect(E0), support = collect(u))
    solve_system!(System(), b)         # free space: the beamlet has one segment
    # that segment on the requested plane, which may also lie before the waist
    return PlaneField(b, 1; size, spacing, origin, axes)
end

@testset "$(r.name)" for r in Conformance.check_source(make)
    @test r.value <= r.limit
end
```

The beamlet starts at its waist with a real amplitude, which is the phase reference of
[`GaussianBeam`](@ref). BeamletOptics starts beams in vacuum, so only `n = 1` applies.

## Reference

```@docs
GaussianBeam
field
check_source
check_propagator
compare
Result
```
