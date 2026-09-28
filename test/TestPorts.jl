using LinearAlgebra
using StaticArrays

const OB = OpticsBase

@testset "Constants" begin
    @test OB.VACUUM_IMPEDANCE == 376.730313668
    p = PlanarPort([0, 0, 0], [0, 0, 1], [1, 0, 0]; refractive_index = 1.5)
    @test OB.power_normalization(p) ≈ 1.5 / (2 * 376.730313668) rtol = 1e-15
    p32 = PlanarPort(Float32[0, 0, 0], Float32[0, 0, 1], Float32[1, 0, 0])
    @test OB.power_normalization(p32) isa Float32
end

@testset "PlanarPort construction" begin
    # Oblique normal and a u that is neither normalized nor orthogonal to it
    p = PlanarPort([1.0, -2.0, 0.5], [1.0, 2.0, 2.0], [3.0, 0.0, 1.0]; refractive_index = 1.33)
    A = OB.local_axes(p)
    u, v, n = A[:, 1], A[:, 2], A[:, 3]
    @test opnorm(transpose(A) * A - I) < 1e-12
    @test cross(u, v) ≈ n atol = 1e-12
    @test n ≈ [1, 2, 2] / 3 atol = 1e-15
    @test OB.normal(p) == n
    @test OB.origin(p) == SVector(1.0, -2.0, 0.5)
    @test OB.refractive_index(p) == 1.33
    # u stays in the plane spanned by the input u and n, on the side of the input u
    @test abs(dot(cross([3.0, 0.0, 1.0], n), u)) < 1e-12
    @test dot(u, [3.0, 0.0, 1.0]) > 0

    # Element type promotion
    @test PlanarPort([0, 0, 0], [0, 0, 1], [1, 0, 0]) isa PlanarPort{Float64}
    @test PlanarPort(Float32[0, 0, 0], Float32[0, 0, 1], Float32[1, 0, 0]) isa
          PlanarPort{Float32}
    @test PlanarPort([0, 0, 0], [0, 0, 1], [1, 0, 0]; refractive_index = 1) isa
          PlanarPort{Float64}

    # Invalid inputs
    @test_throws ArgumentError PlanarPort([0, 0, 0], [0, 0, 0], [1, 0, 0])
    @test_throws ArgumentError PlanarPort([0, 0, 0], [0, 0, 1], [0, 0, 2])
    @test_throws ArgumentError PlanarPort([0, 0, 0], [0, 0, 1], [1e-10, 0, 1])
    @test_throws ArgumentError PlanarPort([0, 0, 0], [0, 0, 1], [0, 0, 0])
    @test_throws ArgumentError PlanarPort([0, 0, 0], [0, 0, 1], [1, 0, 0]; refractive_index = 0)
    @test_throws ArgumentError PlanarPort([0, 0, 0], [0, 0, 1], [1, 0, 0]; refractive_index = -1)
    @test_throws ArgumentError PlanarPort([0, 0], [0, 0, 1], [1, 0, 0])
    @test_throws ArgumentError PlanarPort([0, 0, 0], [0, 0, 1, 0], [1, 0, 0])
end

@testset "Local and global coordinates" begin
    p = PlanarPort([1.0, -2.0, 0.5], [1.0, 2.0, 2.0], [3.0, 0.0, 1.0])
    for ξ in (SVector(0.3, -1.2, 2.5), [0.3, -1.2, 2.5], SVector(1e-3, 4.0, -7.0))
        @test OB.to_local(p, OB.to_global(p, ξ)) ≈ ξ rtol = 1e-12
    end
    r = SVector(4.0, 1.0, -3.0)
    @test OB.to_global(p, OB.to_local(p, r)) ≈ r rtol = 1e-12
    # Two local coordinates give a point on the port plane
    q = OB.to_global(p, [0.7, -0.2])
    @test q ≈ OB.to_global(p, [0.7, -0.2, 0.0]) rtol = 1e-14
    @test abs(OB.to_local(p, q)[3]) < 1e-14
    @test OB.to_local(p, OB.origin(p)) == SVector(0.0, 0.0, 0.0)
    @test_throws ArgumentError OB.to_global(p, [1.0])
    @test_throws ArgumentError OB.to_global(p, [1.0, 2.0, 3.0, 4.0])
end

@testset "Ray basis" begin
    p = PlanarPort([0.0, 0.0, 0.0], [1.0, 2.0, 2.0], [3.0, 0.0, 1.0])
    A = OB.local_axes(p)
    u, v, n = A[:, 1], A[:, 2], A[:, 3]

    # Exactly the port basis along the normal, also for a non-normalized direction
    @test OB.ray_basis(p, OB.normal(p)) == (u, v)
    @test OB.ray_basis(p, 2 * OB.normal(p)) == (u, v)

    # Orthonormal and right-handed for oblique directions
    for d in (normalize(n + 0.3u - 0.5v), normalize(n + 5u + 4v), normalize(n + 1e-9u))
        x̂, ŷ = OB.ray_basis(p, d)
        @test norm(x̂) ≈ 1 atol = 1e-14
        @test norm(ŷ) ≈ 1 atol = 1e-14
        @test abs(dot(x̂, ŷ)) < 1e-14
        @test abs(dot(x̂, d)) < 1e-14
        @test cross(x̂, ŷ) ≈ d atol = 1e-14
    end

    # Richards–Wolf: ê_ρ ↦ ê_θ, ê_φ unchanged (D12)
    φ = deg2rad(45)
    θ = deg2rad(64)
    e_ρ = cos(φ) * u + sin(φ) * v
    e_φ = -sin(φ) * u + cos(φ) * v
    d = sin(θ) * e_ρ + cos(θ) * n
    e_θ = cos(θ) * e_ρ - sin(θ) * n
    x̂, ŷ = OB.ray_basis(p, d)
    @test cos(φ) * x̂ + sin(φ) * ŷ ≈ e_θ atol = 1e-12
    @test -sin(φ) * x̂ + cos(φ) * ŷ ≈ e_φ atol = 1e-12

    # Directions that do not propagate into the port
    @test_throws ArgumentError OB.ray_basis(p, -n)
    @test_throws ArgumentError OB.ray_basis(p, [0.0, 0.0, 0.0])
    # Exactly in the port plane (d·n == 0)
    p_z = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
    @test_throws ArgumentError OB.ray_basis(p_z, [1.0, 0.0, 0.0])
end

@testset "Jones vectors" begin
    p = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
    # Along the normal the Jones basis is (u, v)
    @test OB.jones_to_global(p, [0, 0, 1], [1, 0]) == SVector{3, ComplexF64}(1, 0, 0)
    @test OB.jones_to_global(p, [0, 0, 1], [0, 1]) == SVector{3, ComplexF64}(0, 1, 0)
    # Oblique ray: transverse field with the norm of the Jones vector
    d = normalize([0.4, -0.3, 1.0])
    J = SVector(0.6 + 0.2im, -0.3im)
    E = OB.jones_to_global(p, d, J)
    @test E isa SVector{3, ComplexF64}
    @test abs(sum(E .* d)) < 1e-14
    @test norm(E) ≈ norm(J) rtol = 1e-14
    @test_throws ArgumentError OB.jones_to_global(p, d, [1, 0, 0])
end

@testset "Circular polarization (D2)" begin
    R = OB.circular_jones(+1)
    L = OB.circular_jones(-1)
    @test R ≈ SVector(1, im) / sqrt(2)
    @test L == conj(R)
    @test norm(R) ≈ 1
    # Real field Re(J·exp(−iωt)) of right-circular light rotates from u to v
    field(J, ωt) = real(J * exp(-im * ωt))
    @test normalize(field(R, 0)) ≈ [1, 0] atol = 1e-15
    @test normalize(field(R, π / 2)) ≈ [0, 1] atol = 1e-15
    @test normalize(field(L, π / 2)) ≈ [0, -1] atol = 1e-15
    # Positive helicity: rotation about the propagation direction u × v = n
    p = PlanarPort([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
    n = OB.normal(p)
    E = OB.jones_to_global(p, n, R)
    E0 = real(E)
    E1 = real(E * exp(-im * π / 2))
    @test dot(cross(E0, E1), n) > 0
    @test_throws ArgumentError OB.circular_jones(0)
    @test_throws ArgumentError OB.circular_jones(2)
end
