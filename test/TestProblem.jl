# Tests for PropagationProblem, PropagationSolution and the solve/init entry points.
# Uses test-local mock fields and algorithms only (no RayBundle/SampledField).

# `solve!` and `step!` are the CommonSolve functions, re-exported by OpticsBase.

# --- Mock fields ------------------------------------------------------------------------

struct MockField{P <: AbstractPort} <: AbstractOpticalField{1}
    port::P
    value::Float64
end
OpticsBase.port(f::MockField) = f.port

struct OtherMockField{P <: AbstractPort} <: AbstractOpticalField{1}
    port::P
end
OpticsBase.port(f::OtherMockField) = f.port

const PORT_IN = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
const PORT_OUT = PlanarPort([0.0, 0.0, 0.1], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
const PORT_ELSEWHERE = PlanarPort([0.0, 0.0, 0.2], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])

# Counts __solve/__init calls to check that nothing runs before the compatibility check.
const CALLS = Ref(0)

# --- Mock algorithm implementing __solve -------------------------------------------------

struct DirectAlg <: AbstractPropagationAlgorithm end
OpticsBase.input_representation(::DirectAlg) = MockField
OpticsBase.output_representation(::DirectAlg) = MockField
function OpticsBase.__solve(prob::PropagationProblem, alg::DirectAlg; gain = 1.0)
    CALLS[] += 1
    out = MockField(prob.port_out, gain * prob.field.value)
    return PropagationSolution(out, prob, alg)
end

# --- Mock algorithm implementing __init + solve! + step! ---------------------------------

struct SteppedAlg <: AbstractPropagationAlgorithm end
OpticsBase.input_representation(::SteppedAlg) = MockField
OpticsBase.output_representation(::SteppedAlg) = MockField

mutable struct MockIntegrator{Pr, A}
    prob::Pr
    alg::A
    value::Float64
    gain::Float64
    steps::Int
end
function OpticsBase.__init(prob::PropagationProblem, alg::SteppedAlg; gain = 1.0)
    CALLS[] += 1
    return MockIntegrator(prob, alg, prob.field.value, gain, 0)
end
function OpticsBase.step!(integ::MockIntegrator)
    integ.value *= integ.gain
    integ.steps += 1
    return integ
end
function OpticsBase.solve!(integ::MockIntegrator)
    integ.steps == 0 && step!(integ)
    out = MockField(integ.prob.port_out, integ.value)
    return PropagationSolution(out, integ.prob, integ.alg, (steps = integ.steps,))
end

# --- Faulty mock algorithms ----------------------------------------------------------------

struct WrongPortAlg <: AbstractPropagationAlgorithm end
OpticsBase.input_representation(::WrongPortAlg) = MockField
OpticsBase.output_representation(::WrongPortAlg) = MockField
function OpticsBase.__solve(prob::PropagationProblem, alg::WrongPortAlg)
    return PropagationSolution(MockField(PORT_ELSEWHERE, 0.0), prob, alg)
end

struct WrongTypeAlg <: AbstractPropagationAlgorithm end
OpticsBase.input_representation(::WrongTypeAlg) = MockField
OpticsBase.output_representation(::WrongTypeAlg) = MockField
function OpticsBase.__solve(prob::PropagationProblem, alg::WrongTypeAlg)
    return PropagationSolution(OtherMockField(prob.port_out), prob, alg)
end

struct NoSolutionAlg <: AbstractPropagationAlgorithm end
OpticsBase.input_representation(::NoSolutionAlg) = MockField
OpticsBase.output_representation(::NoSolutionAlg) = MockField
OpticsBase.__solve(prob::PropagationProblem, ::NoSolutionAlg) = prob.field

# Accepts only OtherMockField; __solve counts calls so we can see it never runs.
struct OtherInputAlg <: AbstractPropagationAlgorithm end
OpticsBase.input_representation(::OtherInputAlg) = OtherMockField
OpticsBase.output_representation(::OtherInputAlg) = OtherMockField
function OpticsBase.__solve(prob::PropagationProblem, alg::OtherInputAlg)
    CALLS[] += 1
    return PropagationSolution(OtherMockField(prob.port_out), prob, alg)
end
function OpticsBase.__init(::PropagationProblem, ::OtherInputAlg)
    CALLS[] += 1
    return nothing
end

# Implements nothing but the traits.
struct EmptyAlg <: AbstractPropagationAlgorithm end
OpticsBase.input_representation(::EmptyAlg) = MockField
OpticsBase.output_representation(::EmptyAlg) = MockField

# --- Tests -------------------------------------------------------------------------------

field = MockField(PORT_IN, 2.0)
prob = PropagationProblem(field, :free_space, PORT_OUT)

@testset "types" begin
    @test prob isa PropagationProblem{MockField{typeof(PORT_IN)}, Symbol, typeof(PORT_OUT)}
    @test prob.field === field
    @test prob.system === :free_space
    @test prob.port_out == PORT_OUT
    sol = PropagationSolution(MockField(PORT_OUT, 1.0), prob, DirectAlg())
    @test sol.stats === nothing
    @test sol.prob === prob
    @test sol.alg === DirectAlg()
    sol2 = PropagationSolution(MockField(PORT_OUT, 1.0), prob, DirectAlg(), (steps = 3,))
    @test sol2.stats == (steps = 3,)
    @test_throws MethodError PropagationProblem(1.0, :sys, PORT_OUT)      # not a field
    @test_throws MethodError PropagationProblem(field, :sys, :not_a_port)
end

@testset "compatibility" begin
    @test is_compatible(field, DirectAlg())
    @test !is_compatible(field, OtherInputAlg())
    @test check_compatibility(field, DirectAlg()) === nothing
    err = try
        check_compatibility(field, OtherInputAlg())
    catch e
        e
    end
    @test err isa MissingConverterError
    @test err.field_type == typeof(field)
    @test err.accepted == OtherMockField
    msg = sprint(showerror, err)
    @test occursin("MissingConverterError", msg)
    @test occursin("converter", msg)
    @test occursin(string(typeof(field)), msg)
    @test occursin(string(OtherMockField), msg)
end

@testset "solve via __solve" begin
    CALLS[] = 0
    sol = solve(prob, DirectAlg(); gain = 3.0)
    @test CALLS[] == 1
    @test sol isa PropagationSolution
    @test sol.field isa MockField
    @test sol.field.value == 6.0
    @test OpticsBase.port(sol.field) == PORT_OUT
    @test sol.prob === prob
    @test @inferred(solve(prob, DirectAlg())) isa PropagationSolution
end

@testset "init + solve! and solve via __init fallback" begin
    CALLS[] = 0
    integ = init(prob, SteppedAlg(); gain = 5.0)
    @test integ isa MockIntegrator
    @test CALLS[] == 1
    step!(integ)
    step!(integ)
    sol = solve!(integ)
    @test sol.field.value == 50.0
    @test sol.stats == (steps = 2,)

    # solve falls back to solve!(__init(...)) and checks the output
    CALLS[] = 0
    sol = solve(prob, SteppedAlg(); gain = 5.0)
    @test CALLS[] == 1
    @test sol.field.value == 10.0
    @test sol.stats == (steps = 1,)
    @test OpticsBase.port(sol.field) == PORT_OUT

    # the fallback itself works when called directly
    @test OpticsBase.__solve(prob, SteppedAlg()).field.value == 2.0
end

@testset "incompatible input throws before the solver runs" begin
    CALLS[] = 0
    @test_throws MissingConverterError solve(prob, OtherInputAlg())
    @test_throws MissingConverterError init(prob, OtherInputAlg())
    @test CALLS[] == 0
end

@testset "output validation" begin
    @test_throws ArgumentError solve(prob, WrongPortAlg())
    @test_throws ArgumentError solve(prob, WrongTypeAlg())
    @test_throws ArgumentError solve(prob, NoSolutionAlg())
end

@testset "missing solver implementation" begin
    @test_throws MethodError solve(prob, EmptyAlg())
    @test_throws MethodError init(prob, EmptyAlg())
    @test_throws MethodError init(prob, DirectAlg())      # __solve only: no init support
end
