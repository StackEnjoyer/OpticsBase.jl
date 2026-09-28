const CommonSolve = OpticsBase.CommonSolve

@testset "CommonSolve re-exports" begin
    @test solve === CommonSolve.solve
    @test init === CommonSolve.init
    @test solve! === CommonSolve.solve!
    @test step! === CommonSolve.step!
end
