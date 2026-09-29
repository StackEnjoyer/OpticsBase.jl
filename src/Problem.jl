"""
    PropagationProblem(field, system, port_out)
    PropagationProblem{F <: AbstractOpticalData, S, P <: AbstractPort}

A propagation task: propagate the input `field` through `system` and return the result at
the output port `port_out`. Solve it with `solve(prob, alg)` for an
[`AbstractPropagationAlgorithm`](@ref) `alg`, see [`__solve`](@ref).

# Fields

  - `field::F`: input, an [`AbstractOpticalData`](@ref) (field or ray bundle) given at its
    own input port
  - `system::S`: description of the optical system. Its type is defined by the solver
    package (e.g. a ray tracer's system, a propagation distance); `OpticsBase` does not
    interpret it.
  - `port_out::P`: the [`AbstractPort`](@ref) at which the output field must be given.
    `solve` checks `port(sol.field) == port_out`. The output is given in the medium of
    this port (its refractive index), in the global frame and SI units.

Units, frames and normalization of `field` follow the conventions of its format; the
field is not copied. No checks are done at construction: whether `field` fits an
algorithm depends on the algorithm and is checked by `solve`/`init` (see
[`check_compatibility`](@ref)).
"""
struct PropagationProblem{F <: AbstractOpticalData, S, P <: AbstractPort}
    field::F
    system::S
    port_out::P
end

"""
    PropagationSolution(field, prob, alg, stats = nothing)
    PropagationSolution{F <: AbstractOpticalData, Pr <: PropagationProblem,
                        A <: AbstractPropagationAlgorithm, St}

Result of `solve(prob, alg)`. Solver packages construct it at the end of
[`__solve`](@ref) or of `CommonSolve.solve!` on their integrator; `solve` throws an
`ArgumentError` unless `field` is an instance of `output_representation(alg)` and is given
at `prob.port_out`.

# Fields

  - `field::F`: output, an [`AbstractOpticalData`](@ref) given at `prob.port_out`
    and of type `output_representation(alg)`
  - `prob::Pr`: the [`PropagationProblem`](@ref) that was solved
  - `alg::A`: the [`AbstractPropagationAlgorithm`](@ref) that solved it
  - `stats::St`: solver-specific statistics or diagnostics, `nothing` if there are none
"""
struct PropagationSolution{F <: AbstractOpticalData, Pr <: PropagationProblem,
    A <: AbstractPropagationAlgorithm, St}
    field::F
    prob::Pr
    alg::A
    stats::St
end

# `stats` defaults to `nothing` (documented in the type docstring).
PropagationSolution(field::AbstractOpticalData, prob::PropagationProblem,
    alg::AbstractPropagationAlgorithm) = PropagationSolution(field, prob, alg, nothing)

"""
    MissingConverterError(field_type::Type, accepted::Type) <: Exception

Thrown by [`check_compatibility`](@ref) (and thus by `solve`/`init`) when a field of type
`field_type` is passed to an algorithm that accepts only `accepted`
(`input_representation(alg)`). The field has to be converted explicitly with a converter
before it can be propagated; `OpticsBase` never converts implicitly.

# Fields

  - `field_type::Type`: type of the given field
  - `accepted::Type`: field type the algorithm accepts
"""
struct MissingConverterError <: Exception
    field_type::Type
    accepted::Type
end

function Base.showerror(io::IO, e::MissingConverterError)
    print(io, "MissingConverterError: the algorithm accepts fields of type ", e.accepted,
        " but got a field of type ", e.field_type,
        ". A converter from ", e.field_type, " to ", e.accepted,
        " is missing; convert the field explicitly before propagating it.")
end

# Anything that takes a field as input: propagation algorithms and converters
const _FieldStage = Union{AbstractPropagationAlgorithm, AbstractFieldConverter}

"""
    is_compatible(data::AbstractOpticalData, stage) -> Bool

Returns `true` if the propagation algorithm or converter `stage` accepts `field` as input,
i.e. `field isa input_representation(stage)`. Pipelines use it to check whether stages
fit together without running them.
"""
function is_compatible(field::AbstractOpticalData, stage::_FieldStage)
    return field isa input_representation(stage)
end

"""
    check_compatibility(data::AbstractOpticalData, stage)

Throws a [`MissingConverterError`](@ref) if the propagation algorithm or converter `stage`
does not accept `field` (see [`is_compatible`](@ref)), otherwise returns `nothing`. Called
by `solve`, `init` and [`convert_field`](@ref) before any solver or converter code runs.
"""
function check_compatibility(field::AbstractOpticalData, stage::_FieldStage)
    is_compatible(field, stage) ||
        throw(MissingConverterError(typeof(field), input_representation(stage)))
    return nothing
end

"""
    __solve(prob::PropagationProblem, alg::AbstractPropagationAlgorithm; kwargs...)
        -> PropagationSolution

Solver-side entry point of `solve`. A solver package implements either this method for
its algorithm type, or [`__init`](@ref) plus `CommonSolve.solve!` on the returned
integrator; the fallback `__solve(prob, alg; kwargs...) = solve!(__init(prob, alg; kwargs...))`
connects the two.

# Contract

  - It is called by `solve(prob, alg; kwargs...)` after [`check_compatibility`](@ref), so
    `prob.field isa input_representation(alg)` holds and need not be checked again.
  - It must return a [`PropagationSolution`](@ref) whose `field` is an instance of
    `output_representation(alg)` and is given at `prob.port_out` (global frame, SI units,
    power normalization per [`power_normalization`](@ref) of that port). `solve` checks
    both and throws an `ArgumentError` otherwise.
  - Keyword arguments are passed through from `solve` unchanged.

Users call `solve`, not `__solve`. Without a method for `__solve` or `__init`, `solve`
throws a `MethodError`.
"""
function __solve(prob::PropagationProblem, alg::AbstractPropagationAlgorithm; kwargs...)
    return solve!(__init(prob, alg; kwargs...))
end

"""
    __init(prob::PropagationProblem, alg::AbstractPropagationAlgorithm; kwargs...)
        -> integrator

Solver-side entry point of `init`. A solver package that supports stepping or reuse of
precomputed data implements this method for its algorithm type and returns its own
integrator (or cache) object, together with `CommonSolve.solve!(integrator)` returning a
[`PropagationSolution`](@ref) and optionally `CommonSolve.step!(integrator)`.

# Contract

  - It is called by `init(prob, alg; kwargs...)` after [`check_compatibility`](@ref).
  - `solve!` on the integrator must return a [`PropagationSolution`](@ref) that fulfills
    the same output contract as [`__solve`](@ref).

`OpticsBase` defines no method: solvers that only implement `__solve` do not support
`init`, which then throws a `MethodError`.
"""
function __init end

"""
    solve(prob::PropagationProblem, alg::AbstractPropagationAlgorithm; kwargs...)
        -> PropagationSolution

Propagates `prob.field` through `prob.system` to `prob.port_out` with the solver
algorithm `alg`.

Steps: [`check_compatibility`](@ref)`(prob.field, alg)` (throws a
[`MissingConverterError`](@ref) before any solver code runs), then
[`__solve`](@ref)`(prob, alg; kwargs...)`, then a check that `sol.field isa
output_representation(alg)` and `port(sol.field) == prob.port_out` (`ArgumentError`
otherwise). Keyword arguments are passed to the solver unchanged. Solver packages must not
add methods to `solve`; they implement `__solve` or `__init`.
"""
function CommonSolve.solve(prob::PropagationProblem, alg::AbstractPropagationAlgorithm;
        kwargs...)
    check_compatibility(prob.field, alg)
    sol = __solve(prob, alg; kwargs...)
    _check_solution(sol, prob, alg)
    return sol
end

"""
    init(prob::PropagationProblem, alg::AbstractPropagationAlgorithm; kwargs...)
        -> integrator

Checks [`check_compatibility`](@ref)`(prob.field, alg)` (throws a
[`MissingConverterError`](@ref)) and returns [`__init`](@ref)`(prob, alg; kwargs...)`, the
solver's own integrator, unchanged. `solve!(integrator)` then returns a
[`PropagationSolution`](@ref).

Unlike `solve`, this path does not check the output: the solution returned by
`solve!(init(prob, alg))` is not validated against `output_representation(alg)` and
`prob.port_out`. Use `solve(prob, alg)` for a checked result.
"""
function CommonSolve.init(prob::PropagationProblem, alg::AbstractPropagationAlgorithm;
        kwargs...)
    check_compatibility(prob.field, alg)
    return __init(prob, alg; kwargs...)
end

# Output check of `solve`: representation and port of the returned field.
function _check_solution(sol::PropagationSolution, prob::PropagationProblem,
        alg::AbstractPropagationAlgorithm)
    out = output_representation(alg)
    sol.field isa out ||
        throw(ArgumentError("solve: the solver returned a field of type $(typeof(sol.field)), " *
                            "but output_representation(alg) is $out"))
    port(sol.field) == prob.port_out ||
        throw(ArgumentError("solve: the solver returned a field at a port that differs " *
                            "from prob.port_out"))
    return nothing
end

function _check_solution(sol, ::PropagationProblem, ::AbstractPropagationAlgorithm)
    throw(ArgumentError("solve: __solve must return a PropagationSolution, got $(typeof(sol))"))
end
