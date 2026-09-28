"""
    OpticsBase

Common interface package for an ecosystem of interoperable optical solvers.

Solvers do not know each other, only `OpticsBase`. Optical propagation data is handed
from one solver to the next via exchange formats (subtypes of `AbstractOpticalField`)
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

end # module OpticsBase
