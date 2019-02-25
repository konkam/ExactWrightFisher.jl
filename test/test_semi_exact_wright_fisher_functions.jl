using Random, Distributions

@testset "testing Wright-Fisher semi exact transition and simulation functions" begin
     #Commented lines should be valid when using the A∞ function starting at 0
     Random.seed!(0)
     ref = ExactWrightFisher.Wright_Fisher_exact_transition_arb(0.5, 0.54, 0.75/2, 0.75/2)
     Random.seed!(0)
     @test ExactWrightFisher.Wright_Fisher_exact_transition_with_t005_approx(0.5, 0.54, 0.75/2, 0.75/2) ≈ ref
     Random.seed!(0)
     @test ExactWrightFisher.Wright_Fisher_exact_transition_with_t005_approx(0.5, 0.0054, 0.75/2, 0.75/2) ≈ 0.5670049817247872 atol=10^(-8)
     Random.seed!(0)
     res =  ExactWrightFisher.Wright_Fisher_exact_trajectory_with_t005_approx(0.5, range(0.1, stop = 0.3, step = 0.04), 0.75/2, 0.75/2)
     for i in 1:length(res)
         @test res[i] == [0.5, 0.7391267869778926, 0.8757863681090295, 0.7407268942746179, 0.7312077904902156, 0.6139664321183035, 0.708280926289272][i]
     end

    Random.seed!(0)
    α_vec = (1:4)/4
    res = ExactWrightFisher.Wright_Fisher_K_dim_transition_with_t005_approx(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec, sum(α_vec))
    # for i in 1:length(res)
    #     @test res[i] ≈ [0.06865636342954226, 0.36266311932755124, 0.47374528815354655, 0.0949352290893601][i] atol=10^(-8)
    # end
    for i in 1:length(res)
        @test res[i] ≈ [ 0.02776997207332505, 0.2336806438432984, 0.07839887060526822, 0.6601505134781083][i] atol=10^(-8)
    end
    Random.seed!(0)
    res = ExactWrightFisher.Wright_Fisher_K_dim_trajectory_with_t005_approx(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec)
    for i in 1:length(res[:,1])
        @test res[i, 1] ≈ [ 0.438737367587959, 0.025452755599133924, 0.16906214791905755, 0.36674772889384943][i] atol=10^(-8)
    end
    # for i in 1:length(res[:,end])
    #     @test res[i, end] ≈ [ 0.52047555209178, 0.04761111823097684, 1.3036544865733622e-5, 0.4319002931323774][i] atol=10^(-8)
    # end
    for i in 1:length(res[:,end])
        @test res[i, end] ≈ [ 0.16276396349872138, 0.5893076100430819, 0.10887682518251618, 0.1390516012756806][i] atol=10^(-8)
    end
end
