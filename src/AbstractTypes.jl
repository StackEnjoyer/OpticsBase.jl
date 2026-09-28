"""
    AbstractOpticalField

Supertype of all exchange formats, i.e. representations of an optical field that is
handed from one solver to the next. Every field is given at a port (see
[`AbstractPort`](@ref)) and is one monochromatic component: polychromatic or mutually
incoherent light is a collection of fields.

Concrete formats: [`RayBundle`](@ref), [`SampledField`](@ref).

# Interface

A subtype must implement:

  - [`port`](@ref)`(field)`: the [`AbstractPort`](@ref) at which the field is given
  - [`wavelength`](@ref)`(field)`: vacuum wavelength in \\[m\\]
  - [`total_power`](@ref)`(field)`: power through the port in \\[W\\]
  - [`is_vectorial`](@ref)`(field)`: whether the field carries 3D field vectors
  - [`is_coherent`](@ref)`(field)`: whether all parts of the field add coherently

# Conventions

All formats follow the conventions page of the documentation: SI units, time dependence
exp(−iωt), 3D field vectors in the global frame, phase stored without the carrier
exp(−iωt), power normalization via [`power_normalization`](@ref).
"""
abstract type AbstractOpticalField end

"""
    AbstractPort

Supertype of ports. A port is the surface at which a field is handed from one solver to
the next. It fixes the local frame `(u, v, n)` in global coordinates, the polarization
basis and the refractive index of the medium the field is given in. There is no field
handover without a port, and `OpticsBase` has no global optical axis.

Concrete ports: [`PlanarPort`](@ref).

# Interface

A subtype must implement:

  - [`origin`](@ref)`(port)`: reference point in global coordinates, \\[m\\]
  - [`local_axes`](@ref)`(port)`: 3×3 matrix with the unit columns `(u, v, n)` in global
    coordinates, right-handed (`u × v = n`); `n` points downstream, i.e. fields at the
    port propagate into the half space `n·k > 0`
  - [`normal`](@ref)`(port)`: the unit normal `n`
  - [`refractive_index`](@ref)`(port)`: refractive index of the medium at the port

[`power_normalization`](@ref) is defined generically on top of `refractive_index`.
"""
abstract type AbstractPort end

"""
    AbstractGrid{D}

Supertype of sampling grids with `D` dimensions. A grid lives in the local coordinates
`(ξ, η[, ζ])` of a port (along `u`, `v` and `n`) and carries no offset of its own: the port
origin is the reference point of every grid.

Concrete grids: [`RegularGrid`](@ref).

# Interface

A subtype must implement:

  - `Base.size(grid)`: number of samples per dimension, `NTuple{D,Int}`
  - [`spacing`](@ref)`(grid)`: sample spacing per dimension in \\[m\\]
  - [`coordinates`](@ref)`(grid, d)`: local coordinates of the samples along dimension `d`
    in \\[m\\]
"""
abstract type AbstractGrid{D} end

"""
    AbstractPropagationAlgorithm

Supertype of the algorithm types that solver packages define. An algorithm is passed to
`solve(prob, alg)` together with a [`PropagationProblem`](@ref).

# Interface

A solver package defines a subtype and implements:

  - [`input_representation`](@ref)`(alg)`: the field type the algorithm accepts
  - [`output_representation`](@ref)`(alg)`: the field type the algorithm returns
  - either [`__solve`](@ref)`(prob, alg; kwargs...)`, or [`__init`](@ref)`(prob, alg; kwargs...)`
    together with `CommonSolve.solve!` on the returned integrator

`solve` and `init` themselves are owned by `OpticsBase` and must not be extended by solver
packages; they check the input (and, for `solve`, the output) and dispatch to `__solve` and
`__init`.
"""
abstract type AbstractPropagationAlgorithm end

"""
    AbstractFieldConverter

Supertype of converters: explicit objects that turn a field of one representation into
another, e.g. a [`RayBundle`](@ref) into a [`SampledField`](@ref). Conversions are
almost always approximations, so a converter carries all its parameters (grid, sampling,
number of modes, coherence assumption) and its docstring states what it assumes and what
it conserves. `OpticsBase` never converts implicitly (no `Base.convert` methods).

Concrete converters: [`GaussianBeamletSummation`](@ref).

# Interface

A subtype implements:

  - [`input_representation`](@ref)`(conv)`: the field type the converter accepts
  - [`output_representation`](@ref)`(conv)`: the field type it returns
  - [`__convert_field`](@ref)`(conv, field)`: the conversion itself

Users call [`convert_field`](@ref), which checks the input with
[`check_compatibility`](@ref) and the output against `output_representation`.
"""
abstract type AbstractFieldConverter end

"""
    port(field::AbstractOpticalField) -> AbstractPort

Returns the port at which `field` is given. Part of the [`AbstractOpticalField`](@ref)
interface.
"""
function port end

"""
    wavelength(field::AbstractOpticalField) -> Real

Returns the vacuum wavelength of `field` in \\[m\\]. The wavelength in the medium is
`wavelength(field) / refractive_index(port(field))`. Part of the
[`AbstractOpticalField`](@ref) interface.
"""
function wavelength end

"""
    total_power(field::AbstractOpticalField) -> Real

Returns the power of `field` through its port in \\[W\\]. For sampled fields the power
follows from [`power_normalization`](@ref). Part of the [`AbstractOpticalField`](@ref)
interface.
"""
function total_power end

"""
    is_vectorial(field::AbstractOpticalField) -> Bool

Trait: `true` if `field` carries 3D complex field vectors in the global frame, `false` if
it carries a scalar amplitude. Part of the [`AbstractOpticalField`](@ref) interface.
"""
function is_vectorial end

"""
    is_coherent(field::AbstractOpticalField) -> Bool

Trait: `true` if all parts of `field` (rays, samples) belong to one coherent component and
add as complex amplitudes. Part of the [`AbstractOpticalField`](@ref) interface.
"""
function is_coherent end

"""
    input_representation(stage) -> Type

Trait: the field type that a propagation algorithm or converter `stage` accepts, e.g.
`SampledField{3}` or `Union{RayBundle, SampledField}`. A field is compatible with `stage`
if `field isa input_representation(stage)`, see [`is_compatible`](@ref). There is no
default method: every [`AbstractPropagationAlgorithm`](@ref) and
[`AbstractFieldConverter`](@ref) must implement it.
"""
function input_representation end

"""
    output_representation(stage) -> Type

Trait: the field type that a propagation algorithm returns in `sol.field`, or that a
converter returns from [`convert_field`](@ref). `solve` and `convert_field` check
`field isa output_representation(stage)`. There is no default method: every
[`AbstractPropagationAlgorithm`](@ref) and [`AbstractFieldConverter`](@ref) must
implement it.
"""
function output_representation end
