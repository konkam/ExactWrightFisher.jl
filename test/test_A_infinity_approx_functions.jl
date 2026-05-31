using Random

@testset "testing approximate Ancestor(A∞) functions" begin
    # σ2 returns variance; σ returns standard deviation
    @test ExactWrightFisher.σ2(1., 0.05) ≈ 0.6493907790932706 atol=10^(-8)
    @test ExactWrightFisher.σ2(1., 1) ≈ 2. / 3
    # Edge case: θ == 1 (β == 0 branch in σ2)
    @test ExactWrightFisher.σ2(0.5, 1) ≈ 2/(3*0.5) atol=10^(-8)
    @test ExactWrightFisher.σ2(2.0, 1) ≈ 2/(3*2) atol=10^(-8)
    @test ExactWrightFisher.σ(1.0, 0.5) ≈ sqrt(2/(3*0.5)) atol=10^(-8)
    @test ExactWrightFisher.σ(1.0, 2.0) ≈ sqrt(2/(3*2)) atol=10^(-8)
    Random.seed!(0)
    @test_nowarn ExactWrightFisher.Compute_A∞_approx(sum(1:4), 0.05)
end
