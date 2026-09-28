using OpticsBase
using Test

# Test modules in execution order. Each entry is a file `test/<name>.jl`.
# Run a subset via `Pkg.test(test_args=["TestPorts"])`; no arguments runs all of them.
const TEST_MODULES = [
    "TestInterface",
    "TestPorts",
    "TestGrids",
    "TestSampledField",
    "TestRayBundle",
    "TestProblem",
    "TestConverters",
    "TestBeamletSummation",
    "TestFreeSpace",
    "TestAqua",
]

const SELECTED = isempty(ARGS) ? TEST_MODULES : ARGS

for name in SELECTED
    name in TEST_MODULES ||
        error("Unknown test module \"$name\". Available: $(join(TEST_MODULES, ", "))")
end

# Each test file is evaluated in its own module, so helper names cannot clash.
@testset "OpticsBase.jl" begin
    for name in SELECTED
        @testset "$name" begin
            mod = Module(Symbol(name))
            Core.eval(mod, :(using OpticsBase, Test))
            Base.include(mod, joinpath(@__DIR__, "$name.jl"))
        end
    end
end
