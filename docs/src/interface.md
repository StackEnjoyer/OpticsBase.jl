```@meta
CurrentModule = OpticsBase
```

# Solver interface

Solvers follow the [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl) interface. A
propagation problem bundles the input field, the optical system and the output port; a
solver package supplies the algorithm:

```julia
prob = PropagationProblem(field_in, system, port_out)
sol = solve(prob, alg)   # sol.field isa AbstractOpticalField, given at port_out
```

`OpticsBase` owns `solve` and `init` for propagation problems. They check that the input
field fits the algorithm and, for `solve`, that the solution has the promised
representation and port, then dispatch to the solver's implementation.

## Problems and solutions

```@docs
PropagationProblem
PropagationSolution
```

## Solving

```@docs
CommonSolve.solve(::PropagationProblem, ::AbstractPropagationAlgorithm)
CommonSolve.init(::PropagationProblem, ::AbstractPropagationAlgorithm)
```

## Implementing a solver

```@docs
AbstractPropagationAlgorithm
input_representation
output_representation
__solve
__init
```

## Compatibility checks

```@docs
is_compatible
check_compatibility
MissingConverterError
```
