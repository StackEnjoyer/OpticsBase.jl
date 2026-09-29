# First real chain: BeamletOptics beamlet → Detector → SampledField (BMO's own beamlet sum)
# → AngularSpectrumMethod (WaveOpticsPropagation) → SampledField, checked against the
# analytic Gaussian beam.

import BeamletOptics as BMO
using LinearAlgebra
using StaticArrays
const OB = OpticsBase

const λ = 1.064e-6
const w0 = 0.5e-3
const P0 = 1e-3
const zR = π * w0^2 / λ
const yA = 0.1                 # handover plane (port A) behind the waist
const zAB = zR / 2             # angular-spectrum propagation distance A → B
const NS = 256
const Δ = w0 / 16

w_analytic(y) = w0 * sqrt(1 + (y / zR)^2)

function detector_at(y)
    det = BMO.Detector(0.05)
    BMO.translate3d!(det, [0, y, 0])
    return det
end

function traced_detector(y)
    agb = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0,
        support = [1.0, 0, 0])
    det = detector_at(y)
    BMO.solve_system!(BMO.System([det]), agb)
    return det
end

# Width from the second moment of |E|²: for a Gaussian ⟨ξ²⟩ = w²/4
function second_moment_width(f)
    I = dropdims(sum(abs2, OB.field_array(f); dims = 3); dims = 3)
    ξ = OB.coordinates(OB.grid(f), 1)
    η = OB.coordinates(OB.grid(f), 2)
    s = sum(I)
    wξ = 2 * sqrt(sum(I[i, j] * ξ[i]^2 for i in eachindex(ξ), j in eachindex(η)) / s)
    wη = 2 * sqrt(sum(I[i, j] * η[j]^2 for i in eachindex(ξ), j in eachindex(η)) / s)
    return wξ, wη
end

grid = RegularGrid((NS, NS), (Δ, Δ))
field_A = SampledField(traced_detector(yA), grid)

portA = OB.port(field_A)
portB = PlanarPort(OB.origin(portA) + zAB * OB.normal(portA), OB.normal(portA),
    OB.local_axes(portA)[:, 1])
sol = solve(PropagationProblem(field_A, FreeSpace(), portB), AngularSpectrumMethod())
field_B = sol.field

@testset "handover at port A" begin
    @test field_A isa SampledField{3}
    @test total_power(field_A) ≈ P0 rtol = 1e-6
    wξ, wη = second_moment_width(field_A)
    @test wξ ≈ w_analytic(yA) rtol = 1e-3
    @test wη ≈ w_analytic(yA) rtol = 1e-3
end

@testset "propagated to port B" begin
    @test OB.port(field_B) == portB
    @test OB.origin(portB) ≈ SVector(0.0, yA + zAB, 0.0) atol = 1e-15
    @test total_power(field_B) ≈ P0 rtol = 1e-6
    wξ, wη = second_moment_width(field_B)
    @test wξ ≈ w_analytic(yA + zAB) rtol = 1e-3
    @test wη ≈ w_analytic(yA + zAB) rtol = 1e-3
    # Polarization stays along x (port u): the other components are negligible
    E = OB.field_array(field_B)
    @test maximum(abs, E[:, :, 2:3]) < 1e-9 * maximum(abs, E[:, :, 1])
    # On-axis phase relative to port A: k z − Δ(Gouy)
    c = NS ÷ 2 + 1
    k = 2π / λ
    Δϕ = angle(E[c, c, 1] / OB.field_array(field_A)[c, c, 1])
    expected = k * zAB - (atan((yA + zAB) / zR) - atan(yA / zR))
    @test abs(rem2pi(Δϕ - expected, RoundNearest)) < 1e-3
end
