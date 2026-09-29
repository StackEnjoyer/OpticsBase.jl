# convert_field flow with test-local mock fields and converters.

struct MockA{P <: AbstractPort} <: AbstractOpticalField{1}
    port::P
end
OpticsBase.port(f::MockA) = f.port

struct MockB{P <: AbstractPort} <: AbstractOpticalField{1}
    port::P
end
OpticsBase.port(f::MockB) = f.port

const PORT = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
const CALLS = Ref(0)

# A → B
struct AtoB <: AbstractFieldConverter end
OpticsBase.input_representation(::AtoB) = MockA
OpticsBase.output_representation(::AtoB) = MockB
function OpticsBase.__convert_field(::AtoB, f::MockA)
    CALLS[] += 1
    return MockB(f.port)
end

# Promises B, returns A
struct Liar <: AbstractFieldConverter end
OpticsBase.input_representation(::Liar) = MockA
OpticsBase.output_representation(::Liar) = MockB
OpticsBase.__convert_field(::Liar, f::MockA) = f

# Implements only the traits
struct Hollow <: AbstractFieldConverter end
OpticsBase.input_representation(::Hollow) = MockA
OpticsBase.output_representation(::Hollow) = MockB

@testset "convert_field" begin
    a = MockA(PORT)
    CALLS[] = 0
    b = convert_field(AtoB(), a)
    @test b isa MockB
    @test CALLS[] == 1
    @test is_compatible(a, AtoB())
    @test !is_compatible(b, AtoB())
    @test check_compatibility(a, AtoB()) === nothing
end

@testset "input is checked before the converter runs" begin
    CALLS[] = 0
    err = try
        convert_field(AtoB(), MockB(PORT))
    catch e
        e
    end
    @test err isa MissingConverterError
    @test err.accepted == MockA
    @test CALLS[] == 0
end

@testset "output is checked" begin
    @test_throws ArgumentError convert_field(Liar(), MockA(PORT))
    @test_throws MethodError convert_field(Hollow(), MockA(PORT))
end
