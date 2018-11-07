using Random, Distributions

@testset "testing Wright-Fisher transition and simulation functions" begin
    Random.seed!(0)
    @test ExactWrightFisher.Wright_Fisher_exact_transition(0.5, 0.54, 0.75/2, 0.75/2) ≈ 0.9177951585792223 atol=10^(-8)
    Random.seed!(0)
    res =  ExactWrightFisher.Wright_Fisher_exact_trajectory(0.5, range(0.1, stop = 1, length = 5), 0.75/2, 0.75/2)
    for i in 1:length(res)
        @test res[i] ≈ [0.5, 0.7391267869778926, 0.9831515518978562, 0.804459370375562, 0.779192030871938, 0.6671481291268048][i] atol=10^(-8)
    end
    Random.seed!(0)
    α_vec = (1:4)/4
    res = ExactWrightFisher.Wright_Fisher_K_dim_exact_transition(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec, sum(α_vec))
    for i in 1:length(res)
        @test res[i] ≈ [0.06865636342954226, 0.36266311932755124, 0.47374528815354655, 0.0949352290893601][i] atol=10^(-8)
    end
end
