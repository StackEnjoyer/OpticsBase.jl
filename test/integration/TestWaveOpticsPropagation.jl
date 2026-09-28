# AngularSpectrumMethod with the WaveOpticsPropagation backend against the analytic
# Gaussian beam (plan chain-bmo-fourier, acceptance W2). The input field is built directly
# as a SampledField at the waist, no RayBundle involved.

using LinearAlgebra
using OpticsBase: origin, spacing, coordinates, field_array, port

const W0 = 0.25e-3        # waist radius (1/e² intensity) [m]
const Λ = 1.064e-6        # vacuum wavelength [m]
const NS = 256            # samples per dimension
const DS = W0 / 16        # sample spacing [m]; grid extent 16 w₀

# Scalar Gaussian at its waist, amplitude 1 V/m on axis, on a square grid at `p`.
function waist_field(p)
    g = RegularGrid((NS, NS), (DS, DS))
    ξ = coordinates(g, 1)
    η = coordinates(g, 2)
    E = [complex(exp(-(x^2 + y^2) / W0^2)) for x in ξ, y in η]
    return SampledField(E, g, p, Λ)
end

# 1/e² radius from the second moment along dimension 1: w = 2√⟨ξ²⟩.
function moment_width(f)
    I = dropdims(sum(abs2, field_array(f); dims = 3); dims = 3)
    ξ = coordinates(OpticsBase.grid(f), 1)
    return 2 * sqrt(sum(ξ .^ 2 .* I) / sum(I))
end

const IC = NS ÷ 2 + 1     # index of the port origin (grid center)

@testset "scalar Gaussian, n = $n_med, z = $sgn z_R/2" for n_med in (1.0, 1.5),
    sgn in (1, -1)

    n̂ = [0.0, 0.0, 1.0]
    u = [1.0, 0.0, 0.0]
    p_in = PlanarPort([0.0, 0.0, 0.0], n̂, u; refractive_index = n_med)
    f_in = waist_field(p_in)

    k = 2π * n_med / Λ
    z_R = π * W0^2 * n_med / Λ
    z = sgn * z_R / 2
    p_out = PlanarPort(z * n̂, n̂, u; refractive_index = n_med)

    sol = solve(PropagationProblem(f_in, FreeSpace(), p_out), AngularSpectrumMethod())
    f_out = sol.field
    @test port(f_out) == p_out
    @test OpticsBase.grid(f_out) == OpticsBase.grid(f_in)
    @test !is_vectorial(f_out)

    # Beam width w(z) = w₀ √(1 + (z/z_R)²)
    w_expected = W0 * sqrt(1 + (z / z_R)^2)
    @test isapprox(moment_width(f_out), w_expected; rtol = 1e-3)

    # On-axis phase k z − atan(z/z_R) (input amplitude on axis is real)
    E_c = field_array(f_out)[IC, IC, 1]
    φ_expected = k * z - atan(z / z_R)
    @test abs(angle(E_c * cis(-φ_expected))) < 1e-3

    # Power conservation
    @test isapprox(total_power(f_out), total_power(f_in); rtol = 1e-6)
end

@testset "vectorial field matches the scalar result" begin
    # Oblique port frame, so the global components are all non-trivial.
    n̂ = normalize([1.0, 2.0, 3.0])
    p_in = PlanarPort([0.1, -0.2, 0.3], n̂, [1.0, 0.0, 0.0])
    û = OpticsBase.local_axes(p_in)[:, 1]
    v̂ = OpticsBase.local_axes(p_in)[:, 2]
    e_pol = normalize(û + 0.5im * v̂)   # elliptical, transverse to n̂

    f_scalar = waist_field(p_in)
    E_s = field_array(f_scalar)[:, :, 1]
    E_v = cat((e_pol[c] .* E_s for c in 1:3)...; dims = 3)
    f_vec = SampledField(E_v, OpticsBase.grid(f_scalar), p_in, Λ)
    @test is_vectorial(f_vec)

    z_R = π * W0^2 / Λ
    p_out = PlanarPort(origin(p_in) + (z_R / 2) * n̂, n̂, [1.0, 0.0, 0.0])
    alg = AngularSpectrumMethod()
    out_s = solve(PropagationProblem(f_scalar, FreeSpace(), p_out), alg).field
    out_v = solve(PropagationProblem(f_vec, FreeSpace(), p_out), alg).field
    @test is_vectorial(out_v)

    S = field_array(out_s)[:, :, 1]
    V = field_array(out_v)
    scale = maximum(abs, S)
    for c in 1:3
        @test maximum(abs, V[:, :, c] .- e_pol[c] .* S) ≤ 1e-12 * scale
    end
    @test isapprox(total_power(out_v), total_power(f_vec); rtol = 1e-6)

    # The input field is not modified
    @test field_array(f_vec) == E_v
end
