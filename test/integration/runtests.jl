# Integration tests with the heavy weak dependencies (BeamletOptics, WaveOpticsPropagation,
# OpticSim; the OpticSim modules load it themselves).
# They live in their own environment so the core suite stays light and runs on Julia 1.10.
#
# Run from the repository root (all modules, or a subset by name):
#   julia --project=test/integration -e 'using Pkg; Pkg.instantiate(); include("test/integration/runtests.jl")'
#   julia --project=test/integration -e '...same...' TestChain

using OpticsBase
using BeamletOptics
using WaveOpticsPropagation
using Test

const TEST_MODULES = [
    "TestWaveOpticsPropagation",
    "TestBeamletOptics",
    "TestChain",
    "TestFocus",
    "TestOpticSim",
    "TestRayExchange",
]

const SELECTED = isempty(ARGS) ? TEST_MODULES : ARGS

for name in SELECTED
    name in TEST_MODULES ||
        error("Unknown integration test module \"$name\". Available: $(join(TEST_MODULES, ", "))")
end

# Each test file is evaluated in its own module, so helper names cannot clash. BeamletOptics
# and WaveOpticsPropagation are loaded (extensions active) but not `using`-ed there; test
# files import what they need.
@testset "OpticsBase integration" begin
    for name in SELECTED
        @testset "$name" begin
            mod = Module(Symbol(name))
            Core.eval(mod, :(using OpticsBase, Test))
            Base.include(mod, joinpath(@__DIR__, "$name.jl"))
        end
    end
end
