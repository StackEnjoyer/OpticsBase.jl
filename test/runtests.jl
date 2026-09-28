using OpticsBase
using Test

# Test modules in execution order. Each entry is a file `test/<name>.jl`.
# Run a subset via `Pkg.test(test_args=["TestAqua"])`; no arguments runs all of them.
const TEST_MODULES = [
    "TestAqua",
    "TestInterface",
]

const SELECTED = isempty(ARGS) ? TEST_MODULES : ARGS

for name in SELECTED
    name in TEST_MODULES ||
        error("Unknown test module \"$name\". Available: $(join(TEST_MODULES, ", "))")
end

@testset "OpticsBase.jl" begin
    for name in SELECTED
        @testset "$name" begin
            include(joinpath(@__DIR__, "$name.jl"))
        end
    end
end
