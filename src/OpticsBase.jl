"""
    OpticsBase

Common interface package for an ecosystem of interoperable optical solvers.

Solvers do not know each other, only `OpticsBase`. Optical propagation data is handed
from one solver to the next via exchange formats (subtypes of `AbstractOpticalData`:
fields and ray bundles)
at ports, i.e. surfaces with an explicit frame in global coordinates. `OpticsBase`
contains no solver code: only abstract types, exchange formats, traits, ports,
conventions and generic converters.

The solver interface follows [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl):
`solve`, `init`, `solve!` and `step!` are re-exported from there.

All quantities are in SI units (m, s, W). See the "Conventions" page of the
documentation for the binding time, phase, power and polarization conventions.
"""
module OpticsBase

using LinearAlgebra
using StaticArrays
using CommonSolve: CommonSolve, solve, init, solve!, step!

export solve, init, solve!, step!

# Abstract types and traits
export AbstractOpticalData, AbstractOpticalField, AbstractVectorField, AbstractScalarField
export AbstractRayBundle, AbstractPort, AbstractGrid, AbstractPropagationAlgorithm
export is_vectorial, is_coherent, total_power, input_representation, output_representation

# Ports and exchange formats
export PlanarPort, RegularGrid, SampledField, PlaneWaveSpectrum, RayBundle, PolarizedRayBundle

# Problem interface
export PropagationProblem, PropagationSolution, MissingConverterError
export is_compatible, check_compatibility

# Converters
export AbstractFieldConverter, convert_field
export PlaneWaveSummation, PlaneWaveDecomposition, DebyeWolf

# Free-space propagation (numerics via the WaveOpticsPropagation extension)
export FreeSpace, AngularSpectrumMethod

# Documented API that is not exported to avoid name clashes with solver packages
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("""
        public port, wavelength, VACUUM_IMPEDANCE, power_normalization,
               origin, local_axes, normal, refractive_index, to_local, to_global,
               ray_basis, jones_to_global, circular_jones,
               spacing, coordinates, grid, field_array, positions, directions,
               __solve, __init, __convert_field
        """))
end

include("AbstractTypes.jl")
include("Constants.jl")
include("Ports.jl")
include("Grids.jl")
include("SampledField.jl")
include("RayBundle.jl")
include("PlaneWaveSpectrum.jl")
include("Problem.jl")
include("Converters.jl")
include("FreeSpace.jl")
include("PlaneWaveSummation.jl")
include("PlaneWaveDecomposition.jl")
include("DebyeWolf.jl")

function __init__()
    _register_error_hints()
    _register_fftw_hint()
    _register_delaunay_hint()
    return nothing
end

end # module OpticsBase
