"""
    AbstractOpticalData

Supertype of all exchange formats: the data handed from one solver to the next. Every
exchange format is given at a port (see [`AbstractPort`](@ref)) and describes one
monochromatic component; polychromatic or mutually incoherent light is a collection of
them. There are two branches, ordered by physical model:

  - [`AbstractOpticalField`](@ref)`{N}`: Maxwell fields (wave optics, diffraction
    included), vectorial (`N = 3`) or scalar (`N = 1`).
  - [`AbstractRayBundle`](@ref): rays (geometrical optics, the limit λ → 0), exchanged
    between ray tracers or passed to a ray-to-field converter.

[`PropagationProblem`](@ref), [`convert_field`](@ref) and the compatibility checks accept
any `AbstractOpticalData`; `input_representation` of an algorithm or converter says which
format it needs.

# Interface

A subtype must implement:

  - [`port`](@ref)`(data)`: the [`AbstractPort`](@ref) at which the data is given
  - [`wavelength`](@ref)`(data)`: vacuum wavelength in \\[m\\]

# Conventions

All formats follow the conventions page of the documentation: SI units, time dependence
exp(−iωt), 3D vectors in the global frame, phase stored without the carrier exp(−iωt),
power normalization via [`power_normalization`](@ref).
"""
abstract type AbstractOpticalData end

"""
    AbstractOpticalField{N} <: AbstractOpticalData

Supertype of the field formats: a monochromatic, coherent electromagnetic field
(Maxwell optics) given at a port. `N = 3` is a vectorial field with 3D complex field
vectors in the global frame, `N = 1` a scalar field (polarization dropped). Aliases:
[`AbstractVectorField`](@ref) `= AbstractOpticalField{3}`, [`AbstractScalarField`](@ref)
`= AbstractOpticalField{1}`. Formats differ in the encoding (samples on a screen, plane
waves, …), not in the physics they describe.

Concrete formats: [`SampledField`](@ref), [`PlaneWaveSpectrum`](@ref).

# Interface

A subtype implements the [`AbstractOpticalData`](@ref) interface ([`port`](@ref),
[`wavelength`](@ref)) and:

  - [`total_power`](@ref)`(field)`: power through the port in \\[W\\]
  - [`is_coherent`](@ref)`(field)`: whether all parts of the field add coherently

[`is_vectorial`](@ref) is defined generically from `N`.
"""
abstract type AbstractOpticalField{N} <: AbstractOpticalData end

"""
    AbstractVectorField = AbstractOpticalField{3}

Vectorial fields: 3D complex field vectors \\[V/m\\] in the global frame. Use it in
`input_representation` for algorithms that need the full Maxwell field in any encoding.
"""
const AbstractVectorField = AbstractOpticalField{3}

"""
    AbstractScalarField = AbstractOpticalField{1}

Scalar fields: one complex amplitude \\[V/m\\] per sample, polarization dropped (scalar
approximation).
"""
const AbstractScalarField = AbstractOpticalField{1}

"""
    AbstractRayBundle <: AbstractOpticalData

Supertype of the ray formats: a bundle of rays crossing a port (geometrical optics). Rays
are exchanged between ray tracers, or passed to a converter that leaves geometrical
optics (e.g. [`DebyeWolf`](@ref)); they must not be handed over at caustics (see the
conventions page, "Physical validity at the port").

Concrete formats: [`RayBundle`](@ref) (geometry only), [`PolarizedRayBundle`](@ref)
(with optical path length, power and polarization).

# Interface

A subtype implements the [`AbstractOpticalData`](@ref) interface ([`port`](@ref),
[`wavelength`](@ref)) and:

  - [`positions`](@ref)`(bundle)`: points where the rays cross the port, global frame, \\[m\\]
  - [`directions`](@ref)`(bundle)`: unit propagation directions, global frame
  - `Base.length(bundle)`: number of rays
"""
abstract type AbstractRayBundle <: AbstractOpticalData end

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

  - [`input_representation`](@ref)`(alg)`: the exchange format type the algorithm
    accepts (a field or a ray format)
  - [`output_representation`](@ref)`(alg)`: the exchange format type it returns
  - either [`__solve`](@ref)`(prob, alg; kwargs...)`, or [`__init`](@ref)`(prob, alg; kwargs...)`
    together with `CommonSolve.solve!` on the returned integrator

`solve` and `init` themselves are owned by `OpticsBase` and must not be extended by solver
packages; they check the input (and, for `solve`, the output) and dispatch to `__solve` and
`__init`.
"""
abstract type AbstractPropagationAlgorithm end

"""
    AbstractFieldConverter

Supertype of converters: explicit objects that turn one exchange format
([`AbstractOpticalData`](@ref)) into another, e.g. a [`PlaneWaveSpectrum`](@ref) into a
[`SampledField`](@ref). Conversions are almost always approximations, so a converter
carries all its parameters (grid, sampling, number of modes, coherence assumption) and its
docstring states what it assumes and what it conserves. `OpticsBase` never converts
implicitly (no `Base.convert` methods).

Concrete converters: [`PlaneWaveSummation`](@ref), [`PlaneWaveDecomposition`](@ref),
[`DebyeWolf`](@ref).

# Interface

A subtype implements:

  - [`input_representation`](@ref)`(conv)`: the format type the converter accepts
  - [`output_representation`](@ref)`(conv)`: the format type it returns
  - [`__convert_field`](@ref)`(conv, data)`: the conversion itself

Users call [`convert_field`](@ref), which checks the input with
[`check_compatibility`](@ref) and the output against `output_representation`.
"""
abstract type AbstractFieldConverter end

"""
    port(data::AbstractOpticalData) -> AbstractPort

Returns the port at which `data` (a field or a ray bundle) is given. Part of the
[`AbstractOpticalData`](@ref) interface.
"""
function port end

"""
    wavelength(data::AbstractOpticalData) -> Real

Returns the vacuum wavelength of `data` in \\[m\\]. The wavelength in the medium is
`wavelength(data) / refractive_index(port(data))`. Part of the
[`AbstractOpticalData`](@ref) interface.
"""
function wavelength end

"""
    total_power(field::AbstractOpticalField) -> Real
    total_power(bundle::PolarizedRayBundle) -> Real

Returns the power of `field` through its port in \\[W\\]. For sampled fields the power
follows from [`power_normalization`](@ref). Part of the [`AbstractOpticalField`](@ref)
interface; ray bundles that carry power ([`PolarizedRayBundle`](@ref)) implement it too.
"""
function total_power end

"""
    is_vectorial(field::AbstractOpticalField{N}) -> Bool

Trait: `true` if `field` carries 3D complex field vectors in the global frame (`N = 3`),
`false` if it carries a scalar amplitude (`N = 1`). Defined generically from `N`.
"""
function is_vectorial end

is_vectorial(::AbstractOpticalField{N}) where {N} = N == 3

"""
    is_coherent(field::AbstractOpticalField) -> Bool

Trait: `true` if all parts of `field` (samples, plane waves) belong to one coherent
component and add as complex amplitudes. Part of the [`AbstractOpticalField`](@ref)
interface.
"""
function is_coherent end

"""
    positions(bundle::AbstractRayBundle) -> AbstractVector{<:SVector{3}}

Returns the points where the rays of `bundle` cross its port, global frame, \\[m\\]. Part
of the [`AbstractRayBundle`](@ref) interface.
"""
function positions end

"""
    directions(bundle::AbstractRayBundle) -> AbstractVector{<:SVector{3}}

Returns the unit propagation directions of the rays of `bundle`, global frame. Part of the
[`AbstractRayBundle`](@ref) interface.
"""
function directions end

"""
    input_representation(stage) -> Type

Trait: the exchange format type that a propagation algorithm or converter `stage` accepts,
e.g. `SampledField{3}`, [`AbstractVectorField`](@ref) (any vectorial field) or
[`RayBundle`](@ref). Data is compatible with `stage` if `data isa
input_representation(stage)`, see [`is_compatible`](@ref). There is no default method:
every [`AbstractPropagationAlgorithm`](@ref) and [`AbstractFieldConverter`](@ref) must
implement it.
"""
function input_representation end

"""
    output_representation(stage) -> Type

Trait: the exchange format type that a propagation algorithm returns in `sol.field`, or
that a converter returns from [`convert_field`](@ref). `solve` and `convert_field` check
`data isa output_representation(stage)`. There is no default method: every
[`AbstractPropagationAlgorithm`](@ref) and [`AbstractFieldConverter`](@ref) must
implement it.
"""
function output_representation end
