using Random, Distributions

@testset "testing Wright-Fisher semi exact transition and simulation functions" begin
    Random.seed!(0)
    α_vec = (1:4)/4
    res = ExactWrightFisher.Wright_Fisher_K_dim_transition_with_t005_approx(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec, sum(α_vec))
    for i in 1:length(res)
        @test res[i] ≈ [0.06865636342954226, 0.36266311932755124, 0.47374528815354655, 0.0949352290893601][i] atol=10^(-8)
    end
    Random.seed!(0)
    res = ExactWrightFisher.Wright_Fisher_K_dim_trajectory_with_t005_approx(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec)
    for i in 1:length(res[:,1])
        @test res[i, 1] ≈ [ 0.438737367587959, 0.025452755599133924, 0.16906214791905755, 0.36674772889384943][i] atol=10^(-8)
    end
    for i in 1:length(res[:,end])
        @test res[i, end] ≈ [ 0.52047555209178, 0.04761111823097684, 1.3036544865733622e-5, 0.4319002931323774][i] atol=10^(-8)
    end
end
