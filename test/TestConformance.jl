using LinearAlgebra: dot, norm
const C = OpticsBase.Conformance

passes(results) = all(r -> r.value <= r.limit, results)
failing(results) = [r.name for r in results if r.value > r.limit]

# Exact angular spectrum propagation on the input grid (matrix DFT, no FFT package)
function asm_propagate(f, L)
    nx, ny = size(f.E, 1), size(f.E, 2)
    k = 2π * f.n / f.λ
    dft(N, Δ) = begin
        x = ((0:(N - 1)) .- N ÷ 2) .* Δ
        kx = 2π .* ((0:(N - 1)) .- N ÷ 2) ./ (N * Δ)
        cis.(-kx * x'), kx
    end
    Fx, kx = dft(nx, f.spacing[1])
    Fy, ky = dft(ny, f.spacing[2])
    kz = [sqrt(complex(k^2 - a^2 - b^2)) for a in kx, b in ky]
    E = similar(f.E)
    for c in 1:2
        Ê = Fx * f.E[:, :, c] * transpose(Fy)
        E[:, :, c] = Fx' * (Ê .* cis.(kz .* L)) * conj(Fy) ./ (nx * ny)
    end
    # H from the E-only constructor: the beams of the test travel along the plane normal
    return PlaneField(E, f.spacing, f.origin + L * f.axes[:, 3], f.axes, f.λ; n = f.n)
end

@testset "Analytic source and exact propagator conform" begin
    @test passes(C.check_source(C.field))
    @test passes(C.check_source(C.field; n = 1.45))
    r = C.check_propagator(asm_propagate)
    @test passes(r)
    @test length(r) == 15
end

@testset "Typical convention errors are detected" begin
    conjugated(b, o, A, sz, sp) = (f = C.field(b, o, A, sz, sp);
        PlaneField(conj.(f.E), conj.(f.H), f.spacing, f.origin, f.axes, f.λ; n = f.n))
    @test any(contains("absolute phase"), failing(C.check_source(conjugated)))

    doubled(b, o, A, sz, sp) = (f = C.field(b, o, A, sz, sp);
        PlaneField(f.E .* cis.(angle.(f.E)), f.spacing, f.origin, f.axes, f.λ; n = f.n))
    @test !passes(C.check_source(doubled))

    lefthanded(b, o, A, sz, sp) = C.field(
        C.GaussianBeam(; waist = b.waist, direction = b.direction, λ = b.λ, w0 = b.w0,
            P = b.P, n = b.n, jones = conj.(b.jones)), o, A, sz, sp)
    @test failing(C.check_source(lefthanded)) ==
          ["axis +y, right-circular, plane at the waist: field shape and polarization"]

    reversed_H(b, o, A, sz, sp) = (f = C.field(b, o, A, sz, sp);
        PlaneField(f.E, -f.H, f.spacing, f.origin, f.axes, f.λ; n = f.n))
    @test all(contains(r"power|backward"), failing(C.check_source(reversed_H)))

    shifted(b, o, A, sz, sp) = C.field(b, o + 1e-6 * A[:, 1], A, sz, sp)
    @test any(contains("plane as requested"), failing(C.check_source(shifted)))

    paraxial_power(f, L) = (g = asm_propagate(f, L);
        PlaneField(1.01 .* g.E, g.spacing, g.origin, g.axes, g.λ; n = g.n))
    @test any(contains("power"), failing(C.check_propagator(paraxial_power)))
end

@testset "Results and reference beam" begin
    r = C.Result("x", 0.5, 1.0)
    @test sprint(show, r) == "pass: x (0.5 <= 1.0)"
    @test startswith(sprint(show, C.Result("x", 2.0, 1.0)), "FAIL")
    b = C.GaussianBeam(; waist = zeros(3), direction = [0, 0, 2], λ = 1e-6, w0 = 20e-6)
    @test b.direction == [0, 0, 1]
    f = C.field(b, zeros(3), [1.0 0 0; 0 1 0; 0 0 1], (128, 128), (2.5e-6, 2.5e-6))
    @test power(f) ≈ b.P rtol = 1e-9
    @test f.E[65, 65, 1] ≈ sqrt(4 * VACUUM_IMPEDANCE * b.P / (π * b.w0^2))
    @test_throws ArgumentError C.field(b, zeros(3), [1.0 0 0; 0 0 1; 0 -1 0], (8, 8),
        (1e-6, 1e-6))
end
