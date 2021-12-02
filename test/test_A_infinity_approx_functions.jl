using Random

@testset "testing approximate Ancestor(A∞) functions" begin
    @test ExactWrightFisher.σ(1., 0.05) ≈ 0.6493907790932706 atol=10^(-8)
    @test ExactWrightFisher.σ(1., 1) ≈ 2. / 3
    Random.seed!(0)
    @test_nowarn ExactWrightFisher.Compute_A∞_approx(sum(1:4), 0.05)
end
