using LinearAlgebra: I

const Z0 = VACUUM_IMPEDANCE
const AXES = Matrix{Float64}(I, 3, 3)
# BeamletOptics-like frame: n = +y, u = +x, v = n × u = −z
const AXES_Y = [1.0 0 0; 0 0 1; 0 -1 0]

gaussian(u, v, w, E0) = E0 * exp(-(u^2 + v^2) / w^2)

function gaussian_field(; N = 128, w = 10e-6, E0 = 3.0, n = 1.45, axes = AXES, R = Inf)
    Δ = 10w / N
    c = ((0:(N - 1)) .- N ÷ 2) .* Δ
    E = [gaussian(u, v, w, E0) for u in c, v in c]
    return PlaneField(E, (Δ, Δ), (1.0, 2.0, 3.0), axes, 1e-6; n, R)
end

@testset "Gaussian beam, E-only constructor" begin
    w, E0, n = 10e-6, 3.0, 1.45
    f = gaussian_field(; w, E0, n)
    @test f isa PlaneField{Float64, Array{ComplexF64, 3}}
    @test all(iszero, f.E[:, :, 2])                         # scalar → Eu
    @test f.H[:, :, 2] ≈ n / Z0 .* f.E[:, :, 1]
    @test power(f) ≈ n / (2Z0) * E0^2 * π * w^2 / 2 rtol = 1e-9
    fw, bw = forward(f), backward(f)
    @test fw.E ≈ f.E && fw.H ≈ f.H
    @test maximum(abs, bw.E) < 1e-12 * E0
    @test coordinates(f, 1)[65] == 0
    @test step(coordinates(f, 2)) ≈ f.spacing[2]
end

@testset "Counter-propagating waves" begin
    n, N = 1.3, 8
    Y = n / Z0
    Ef = fill(1.0 + 0im, N, N, 2);
    Ef[:, :, 2] .= 0
    Eb = fill(0.5im, N, N, 2);
    Eb[:, :, 1] .= 0
    # forward: H = Y n × E; backward: H = −Y n × E
    Hf = zeros(ComplexF64, N, N, 2);
    Hf[:, :, 2] .= Y .* Ef[:, :, 1]
    Hb = zeros(ComplexF64, N, N, 2);
    Hb[:, :, 1] .= Y .* Eb[:, :, 2]
    geo = ((1e-6, 1e-6), zeros(3), AXES, 1e-6)
    b = PlaneField(Eb, Hb, geo...; n)
    @test power(b) < 0
    @test power(b) ≈ -Y / 2 * 0.25 * N^2 * 1e-12
    f = PlaneField(Ef .+ Eb, Hf .+ Hb, geo...; n)
    @test forward(f).E ≈ Ef && forward(f).H ≈ Hf
    @test backward(f).E ≈ Eb && backward(f).H ≈ Hb
    @test power(forward(f)) + power(backward(f)) ≈ power(f)
    # standing wave with equal amplitudes carries no net power
    s = PlaneField(Ef .+ Ef, Hf .- Hf, geo...; n)
    @test power(s) == 0
    @test forward(s).E ≈ Ef && backward(s).E ≈ Ef
end

@testset "Rigid motion" begin
    f = gaussian_field()
    g = PlaneField(f.E, f.H, f.spacing, (0.0, 5.0, 0.0), AXES_Y, f.λ; n = f.n)
    @test g.E === f.E && g.H === f.H                         # no copy, no rotation
    @test power(g) == power(f)
end

@testset "Reference phase" begin
    λ, n = 1e-6, 1.2
    k = 2π / λ * n
    for R in (5e-3, -5e-3)
        f = gaussian_field(; R, n)
        P = reference_phase(f)
        @test P[65, 65] == 1
        u, v = coordinates(f, 1)[10], coordinates(f, 2)[100]
        # the reference formula cancels (√(ρ² + R²) ≈ |R|), hence the loose tolerance
        @test P[10, 100] ≈ cis(sign(R) * k * (sqrt(u^2 + v^2 + R^2) - abs(R))) rtol = 1e-9
        # the reference phase cancels in the Poynting flux of given E and H
        g = PlaneField(f.E, f.H, f.spacing, f.origin, f.axes, f.λ; n = f.n)
        @test power(f) == power(g)
    end
    @test all(==(1), reference_phase(gaussian_field()))
    # converging: phase decreases outwards
    @test angle(reference_phase(gaussian_field(; R = -1e-2))[80, 65]) < 0
end

@testset "Element and array types" begin
    E = rand(ComplexF32, 4, 6, 2)
    f = PlaneField(E, (1.0f-6, 2.0f-6), SVector{3, Float32}(0, 0, 0), Float32.(AXES), 1.0f-6)
    @test f isa PlaneField{Float32}
    @test power(f) isa Float32
    @test forward(f) isa PlaneField{Float32}
    @test length(coordinates(f, 2)) == 6
    # views are copied into the array type of `similar`, arrays of that type are not
    big = rand(ComplexF64, 5, 6, 2)
    Hv = view(copy(big), 1:4, :, :)
    fv = PlaneField(view(big, 1:4, :, :), Hv, (1e-6, 1e-6), zeros(3), AXES, 1e-6)
    @test fv.E isa Array{ComplexF64, 3}
    @test fv.E == big[1:4, :, :] && fv.H == Hv
    H = rand(ComplexF64, 5, 6, 2)
    @test PlaneField(big, H, (1e-6, 1e-6), zeros(3), AXES, 1e-6).H === H
    # real input becomes complex
    @test PlaneField(rand(4, 4, 2), (1e-6, 1e-6), zeros(3), AXES, 1e-6).E isa
          Array{ComplexF64, 3}
end

@testset "Constructor errors" begin
    E = zeros(ComplexF64, 4, 4, 2)
    geo(; axes = AXES, λ = 1e-6, spacing = (1e-6, 1e-6)) = (spacing, zeros(3), axes, λ)
    @test_throws DimensionMismatch PlaneField(E, zeros(4, 5, 2), geo()...)
    @test_throws DimensionMismatch PlaneField(zeros(4, 4, 3), geo()...)
    @test_throws DimensionMismatch PlaneField(zeros(4, 4, 3), zeros(4, 4, 3), geo()...)
    @test_throws ArgumentError PlaneField(E, geo(; axes = 2AXES)...)
    @test_throws ArgumentError PlaneField(E, geo(; axes = [1.0 0 0; 0 1 0; 0 0 -1])...)
    # orthonormal to within 1e-6: single precision axes pass, a skew of 1e-5 does not
    c, s = cos(0.3f0), sin(0.3f0)
    @test PlaneField(E, geo(; axes = Float32[c 0 s; 0 1 0; -s 0 c])...).axes isa
          SMatrix{3, 3, Float64}
    @test_throws ArgumentError PlaneField(E, geo(; axes = [1.0 1e-5 0; 0 1 0; 0 0 1])...)
    @test_throws ArgumentError PlaneField(E, geo(; λ = 0)...)
    @test_throws ArgumentError PlaneField(E, geo(; spacing = (1e-6, 0))...)
    @test_throws ArgumentError PlaneField(E, geo()...; n = 0)
    @test_throws ArgumentError PlaneField(E, geo()...; R = 0)
end

@testset "Reference directions: spherical wave, exact H" begin
    # A diverging spherical wave from the center of the reference sphere, sampled out to
    # about 40° off axis. The E-only constructor with R must reproduce its exact H, and
    # the split must find no backward part.
    λ, n, R = 1e-6, 1.3, 50e-6
    k, Y = 2π / λ * n, n / Z0
    N, Δ = 161, 0.5e-6
    c = ((0:(N - 1)) .- N ÷ 2) .* Δ
    Et = zeros(ComplexF64, N, N, 2)
    Ht = zeros(ComplexF64, N, N, 2)
    for (j, η) in enumerate(c), (i, ξ) in enumerate(c)
        p = [ξ, η, R]                             # sample relative to the sphere center
        r = norm(p)
        d = p / r
        pol = [0.3, 1.0, 0] - dot(d, [0.3, 1.0, 0]) * d    # transverse polarization
        E3 = pol / r * cis(k * (r - R))           # stored field: reference phase removed
        H3 = Y * cross(d, E3)
        Et[i, j, :] = E3[1:2]
        Ht[i, j, :] = H3[1:2]
    end
    geo = ((Δ, Δ), zeros(3), AXES, λ)
    f = PlaneField(Et, geo...; n, R)
    @test f.H ≈ Ht rtol = 1e-12
    exact = PlaneField(Et, Ht, geo...; n, R)
    @test forward(exact).E ≈ Et rtol = 1e-12
    @test maximum(abs, backward(exact).E) < 1e-12 * maximum(abs, Et)
    # without the reference sphere the plane-normal rule is visibly wrong off axis
    @test !isapprox(PlaneField(Et, geo...; n).H, Ht; rtol = 1e-2)
end

@testset "Oblique plane wave: documented split error" begin
    # s-polarized plane wave at θ to n with its exact H: forward/backward assume
    # propagation along n and leak the backward amplitude (1 − cos θ)/2.
    θ, λ, n = 0.5, 1e-6, 1.0
    k, Y = 2π / λ * n, n / Z0
    N, Δ = 32, 0.1e-6
    c = ((0:(N - 1)) .- N ÷ 2) .* Δ
    ψ = [cis(k * sin(θ) * ξ) for ξ in c, η in c]
    E = cat(zero(ψ), ψ; dims = 3)                 # E along v
    H = cat(-Y * cos(θ) .* ψ, zero(ψ); dims = 3)  # H = Y d × E with d = (sin θ, 0, cos θ)
    f = PlaneField(E, H, (Δ, Δ), zeros(3), AXES, λ; n)
    @test backward(f).E[:, :, 2] ≈ (1 - cos(θ)) / 2 .* ψ rtol = 1e-12
    @test power(f) ≈ Y / 2 * cos(θ) * N^2 * Δ^2 rtol = 1e-12
end

@testset "Array and geometry precision are independent" begin
    E = rand(ComplexF32, 8, 8, 2)
    f = PlaneField(E, (1e-6, 1e-6), (1.0, 2.0, 3.0), AXES, 1e-6; R = 1e-3)
    @test f isa PlaneField{Float64, Array{ComplexF32, 3}}
    @test f.E === E || f.E == E
    @test eltype(f.H) == ComplexF32
    @test eltype(reference_phase(f)) == ComplexF32
    @test eltype(forward(f).E) == ComplexF32
    @test f.origin == SVector(1.0, 2.0, 3.0)
end

@testset "Validation rejects non-finite input" begin
    E = zeros(ComplexF64, 4, 4, 2)
    @test_throws ArgumentError PlaneField(E, (1e-6, 1e-6), zeros(3), AXES, Inf)
    @test_throws ArgumentError PlaneField(E, (1e-6, Inf), zeros(3), AXES, 1e-6)
    @test_throws ArgumentError PlaneField(E, (1e-6, 1e-6), [0.0, NaN, 0], AXES, 1e-6)
    @test_throws ArgumentError PlaneField(E, (1e-6, 1e-6), zeros(3), AXES, 1e-6; n = Inf)
    @test_throws ArgumentError PlaneField(E, (1e-6, 1e-6), zeros(3), AXES, 1e-6; R = NaN)
    @test PlaneField(E, (1e-6, 1e-6), zeros(3), AXES, 1e-6; R = -Inf).R == -Inf
end

@testset "GPU-style arrays (JLArrays, no scalar indexing)" begin
    E = JLArray(rand(ComplexF32, 8, 6, 2))
    f = PlaneField(E, (1e-6, 2e-6), zeros(3), AXES, 1e-6; R = 1e-4)
    @test f.E isa JLArray && f.H isa JLArray
    @test reference_phase(f) isa JLArray
    @test forward(f).E isa JLArray
    @test power(f) ≈ power(PlaneField(Array(f.E), Array(f.H), f.spacing, f.origin, f.axes,
        f.λ; R = f.R)) rtol = 1e-5
    @test Array(forward(f).E) ≈ Array(f.E) rtol = 1e-4
    @test PlaneField(Array(E)[:, :, 1] |> JLArray, (1e-6, 2e-6), zeros(3), AXES, 1e-6).E isa
          JLArray
end
