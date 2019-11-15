using Random, Distributions

@testset "testing Wright-Fisher semi exact transition and simulation functions" begin
     #Commented lines should be valid when using the A∞ function starting at 0
     Random.seed!(0)
     ref = ExactWrightFisher.Wright_Fisher_exact_transition_arb(0.5, 0.54, 0.75/2, 0.75/2)
     Random.seed!(0)
     @test ExactWrightFisher.Wright_Fisher_exact_transition_with_t005_approx(0.5, 0.54, 0.75/2, 0.75/2) ≈ ref

     @test_nowarn  ExactWrightFisher.Wright_Fisher_exact_transition_with_t005_approx(0.5, 0., 0.75/2, 0.75/2)

     @test_nowarn  ExactWrightFisher.Wright_Fisher_exact_transition_with_t005_approx(0.5, 0.05, 0.75/2, 0.75/2)


     # Random.seed!(0)
     # ref = ExactWrightFisher.Wright_Fisher_exact_transition_arb(0.5, 0.0054, 0.75/2, 0.75/2)
     # Random.seed!(0)
     # @test ExactWrightFisher.Wright_Fisher_exact_transition_with_t005_approx(0.5, 0.0054, 0.75/2, 0.75/2) ≈ ref  atol=10^(-2)
     # @test ExactWrightFisher.Wright_Fisher_exact_transition_with_t005_approx(0.5, 0.0054, 0.75/2, 0.75/2) ≈ 0.5670049817247872 atol=10^(-8)
     # Random.seed!(0)
     # res =  ExactWrightFisher.Wright_Fisher_exact_trajectory_with_t005_approx(0.5, range(0.1, stop = 0.3, step = 0.04), 0.75/2, 0.75/2)
     # for i in 1:length(res)
     #     @test res[i] == [0.5, 0.7391267869778926, 0.8757863681090295, 0.7407268942746179, 0.7312077904902156, 0.6139664321183035, 0.708280926289272][i]
     # end

    Random.seed!(0)
    α_vec = (1:4)/4
    res = ExactWrightFisher.Wright_Fisher_K_dim_transition_with_t005_approx(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec, sum(α_vec))
    Random.seed!(0)
    ref = ExactWrightFisher.Wright_Fisher_K_dim_exact_transition_arb(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec, sum(α_vec))
    for i in 1:length(res)
        @test res[i] ≈ ref[i]
    end

     @test_nowarn ExactWrightFisher.Wright_Fisher_K_dim_exact_transition_arb(rand(Dirichlet(α_vec |> collect)), 0.0, α_vec, sum(α_vec))
     @test_nowarn ExactWrightFisher.Wright_Fisher_K_dim_exact_transition_arb(rand(Dirichlet(α_vec |> collect)), 0.05, α_vec, sum(α_vec))

    # Random.seed!(0)
    # α_vec = (1:4)/4
    # res = ExactWrightFisher.Wright_Fisher_K_dim_transition_with_t005_approx(rand(Dirichlet(α_vec |> collect)), 0.005, α_vec, sum(α_vec))
    # Random.seed!(0)
    # ref = ExactWrightFisher.Wright_Fisher_K_dim_exact_transition_arb(rand(Dirichlet(α_vec |> collect)), 0.005, α_vec, sum(α_vec))
    # for i in 1:length(res)
    #     @test res[i] ≈ ref[i]
    # end
end
