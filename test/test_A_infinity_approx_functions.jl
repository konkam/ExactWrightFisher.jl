using Random

@testset "testing approximate Ancestor(A∞) functions" begin
    @test ExactWrightFisher.σ(1., 0.05) ≈ 0.6493907790932706 atol=10^(-8)
    Random.seed!(0)
    @test ExactWrightFisher.Compute_A∞_approx(sum(1:4), 0.05) == 45
end
