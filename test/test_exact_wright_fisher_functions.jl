using Random, Distributions

@testset "testing Wright-Fisher transition and simulation functions" begin
    #Commented lines should be valid when using the A∞ function starting at 0
    Random.seed!(0)
    @test ExactWrightFisher.Wright_Fisher_exact_transition(0.5, 0.54, 0.75/2, 0.75/2) ≈ 0.9177951585792223 atol=10^(-8)
    Random.seed!(0)
    res =  ExactWrightFisher.Wright_Fisher_exact_trajectory(0.5, range(0.1, stop = 1, length = 5), 0.75/2, 0.75/2)
    # for i in 1:length(res)
    #     @test res[i] ≈ [0.5, 0.7391267869778926, 0.9831515518978562, 0.804459370375562, 0.779192030871938, 0.6671481291268048][i] atol=10^(-8)
    # end
    for i in 1:length(res)
        @test res[i] ≈ [ 0.5, 0.7391267869778926, 0.9850784214783245, 0.8237514810336812, 0.8039797099052189, 0.6671481291268048][i] atol=10^(-8)
    end
    @test_nowarn ExactWrightFisher.Wright_Fisher_exact_trajectory(0.5, range(0.1, stop = 1, length = 5), 0.75/2, 0.75/2; use_progress_meter = true)

    Random.seed!(0)
    α_vec = (1:4)/4
    res = ExactWrightFisher.Wright_Fisher_K_dim_exact_transition(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec, sum(α_vec))
    # for i in 1:length(res)
    #     @test res[i] ≈ [0.06865636342954226, 0.36266311932755124, 0.47374528815354655, 0.0949352290893601][i] atol=10^(-8)
    # end
    for i in 1:length(res)
        @test res[i] ≈ [ 0.02776997207332505, 0.2336806438432984, 0.07839887060526822, 0.6601505134781083][i] atol=10^(-8)
    end
    Random.seed!(0)
    res = ExactWrightFisher.Wright_Fisher_K_dim_exact_trajectory(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec)
    for i in 1:length(res[:,1])
        @test res[i, 1] ≈ [0.438737367587959, 0.025452755599133924, 0.16906214791905755, 0.36674772889384943][i] atol=10^(-8)
    end
    # for i in 1:length(res[:,end])
    #     @test res[i, end] ≈ [0.52047555209178, 0.04761111823097684, 1.3036544865733622e-5, 0.4319002931323774][i] atol=10^(-8)
    # end
    for i in 1:length(res[:,end])
        @test res[i, end] ≈ [ 0.16276396349872138, 0.5893076100430819, 0.10887682518251618, 0.1390516012756806][i] atol=10^(-8)
    end
    @test_nowarn ExactWrightFisher.Wright_Fisher_K_dim_exact_trajectory(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec; use_progress_meter = true)
end
