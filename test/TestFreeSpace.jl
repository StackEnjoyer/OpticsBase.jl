# FreeSpace / AngularSpectrumMethod without a backend: construction, traits, geometry
# checks and the missing-backend error. WaveOpticsPropagation is never loaded here; the
# successful path is tested in test/integration/TestWaveOpticsPropagation.jl.

using LinearAlgebra

# A port type that is not a PlanarPort (only used to trigger the geometry check).
struct DummyPort <: OpticsBase.AbstractPort end

const λ = 1.064e-6
const N_GRID = 16
const Δ = 10e-6

port_in = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
grid2 = RegularGrid((N_GRID, N_GRID), (Δ, Δ))
field = SampledField(ones(ComplexF64, N_GRID, N_GRID), grid2, port_in, λ)
alg = AngularSpectrumMethod()

@testset "construction and traits" begin
    @test alg.padding === true
    @test alg.pad_factor === 2
    @test alg.bandlimit === true
    a = AngularSpectrumMethod(; padding = false, pad_factor = 1, bandlimit = false)
    @test (a.padding, a.pad_factor, a.bandlimit) === (false, 1, false)
    @test_throws ArgumentError AngularSpectrumMethod(; pad_factor = 0)
    @test alg isa AbstractPropagationAlgorithm
    @test input_representation(alg) == SampledField{<:Any, 2}
    @test output_representation(alg) == SampledField{<:Any, 2}
    @test is_compatible(field, alg)
    fieldv = SampledField(ones(ComplexF64, N_GRID, N_GRID, 3), grid2, port_in, λ)
    @test is_compatible(fieldv, alg)
end

@testset "volume grid is rejected" begin
    grid3 = RegularGrid((4, 4, 4), (Δ, Δ, Δ))
    f3 = SampledField(ones(ComplexF64, 4, 4, 4), grid3, port_in, λ)
    @test !is_compatible(f3, alg)
    prob = PropagationProblem(f3, FreeSpace(), port_in)
    @test_throws MissingConverterError solve(prob, alg)
end

@testset "invalid geometry throws before the backend" begin
    z = 0.1
    bad_ports = [
        # tilted normal
        PlanarPort([0.0, 0.0, z], [0.0, sind(1), cosd(1)], [1.0, 0.0, 0.0]),
        # anti-parallel normal
        PlanarPort([0.0, 0.0, z], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]),
        # same plane orientation, rotated u
        PlanarPort([0.0, 0.0, z], [0.0, 0.0, 1.0], [cosd(1), sind(1), 0.0]),
        # u flipped
        PlanarPort([0.0, 0.0, z], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]),
        # lateral shift
        PlanarPort([1e-6, 0.0, z], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        # lateral shift only (z = 0)
        PlanarPort([0.0, 1e-6, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        # different medium
        PlanarPort([0.0, 0.0, z], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]; refractive_index = 1.5),
    ]
    for p in bad_ports
        @test_throws ArgumentError solve(PropagationProblem(field, FreeSpace(), p), alg)
    end
    # Non-planar output port and non-planar input port
    @test_throws ArgumentError solve(PropagationProblem(field, FreeSpace(), DummyPort()),
        alg)
    f_dummy = SampledField(ones(ComplexF64, N_GRID, N_GRID), grid2, DummyPort(), λ)
    @test_throws ArgumentError solve(PropagationProblem(f_dummy, FreeSpace(), port_in),
        alg)
end

@testset "valid geometry reaches the backend" begin
    # Without WaveOpticsPropagation the backend has no method: a MethodError with a hint.
    n̂ = normalize([1.0, 2.0, 3.0])
    p_in = PlanarPort([0.1, -0.2, 0.3], n̂, [1.0, 0.0, 0.0]; refractive_index = 1.5)
    f_in = SampledField(ones(ComplexF64, N_GRID, N_GRID, 3), grid2, p_in, λ)
    for z in (0.05, -0.05, 0.0)
        p_out = PlanarPort(OpticsBase.origin(p_in) + z * n̂, n̂, [1.0, 0.0, 0.0];
            refractive_index = 1.5)
        prob = PropagationProblem(f_in, FreeSpace(), p_out)
        e = try
            solve(prob, alg)
            nothing
        catch err
            err
        end
        @test e isa MethodError
        @test e.f === OpticsBase._angular_spectrum
        msg = sprint(showerror, e)
        @test occursin("using WaveOpticsPropagation", msg)
    end
end
