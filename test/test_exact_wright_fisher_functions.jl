using Random, Distributions

@testset "testing Wright-Fisher transition and simulation functions" begin
    #Commented lines should be valid when using the A∞ function starting at 0
    seed = 1
    Random.seed!(seed);
    ref = ExactWrightFisher.Wright_Fisher_exact_transition_arb(0.5, 0.54, 0.75/2, 0.75/2)
    Random.seed!(seed);
    res = ExactWrightFisher.Wright_Fisher_exact_transition(0.5, 0.54, 0.75/2, 0.75/2)
    @test res ≈ ref

################ This is a case where the error plays a role

    # @test res ≈ ref In this case there is a difference

    Random.seed!(seed);
    ref = ExactWrightFisher.Wright_Fisher_exact_transition_arb(0.5, 0.05, 2, 2; verbose = true)
    Random.seed!(seed);
    res = ExactWrightFisher.Wright_Fisher_exact_transition(0.5, 0.05, 2, 2)
    @test res ≈ ref

    Random.seed!(seed)
    ref =  ExactWrightFisher.Wright_Fisher_exact_trajectory_arb(0.5, range(0.1, stop = 1, length = 5), 0.75/2, 0.75/2)
    Random.seed!(seed)
    res =  ExactWrightFisher.Wright_Fisher_exact_trajectory(0.5, range(0.1, stop = 1, length = 5), 0.75/2, 0.75/2)
    # for i in 1:length(res)
    #     @test res[i] ≈ [0.5, 0.7391267869778926, 0.9831515518978562, 0.804459370375562, 0.779192030871938, 0.6671481291268048][i] atol=10^(-8)
    # end
    for i in 1:length(res)
        @test res[i] ≈ ref[i]
    end
    # @test_nowarn ExactWrightFisher.Wright_Fisher_exact_trajectory(0.5, range(0.1, stop = 1, length = 5), 0.75/2, 0.75/2; use_progress_meter = false)
    # @test_nowarn ExactWrightFisher.Wright_Fisher_exact_trajectory_arb(0.5, range(0.1, stop = 1, length = 5), 0.75/2, 0.75/2; use_progress_meter = false)

    Random.seed!(0)
    α_vec = (1:4)/4
    ref = ExactWrightFisher.Wright_Fisher_K_dim_exact_transition_arb(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec, sum(α_vec))
    Random.seed!(0)
    res = ExactWrightFisher.Wright_Fisher_K_dim_exact_transition(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec, sum(α_vec))

    for i in 1:length(res)
        @test res[i] ≈ ref[i]
    end
    @test_nowarn ExactWrightFisher.Wright_Fisher_K_dim_exact_transition(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec)
    @test_nowarn ExactWrightFisher.Wright_Fisher_K_dim_exact_transition_arb(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec)

    Random.seed!(0)
    ref = ExactWrightFisher.Wright_Fisher_K_dim_exact_trajectory_arb(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec)

    Random.seed!(0)
    res = ExactWrightFisher.Wright_Fisher_K_dim_exact_trajectory(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec)

    for i in 1:size(ref,1)
        for j in 1:size(ref,2)
            @test res[i,j] ≈ ref[i,j]
        end
    end
    # for i in 1:length(res[:,end])
    #     @test res[i, end] ≈ [0.52047555209178, 0.04761111823097684, 1.3036544865733622e-5, 0.4319002931323774][i] atol=10^(-8)
    # end
    # for i in 1:length(res[:,end])
    #     @test res[i, end] ≈ [ 0.16276396349872138, 0.5893076100430819, 0.10887682518251618, 0.1390516012756806][i] atol=10^(-8)
    # end
    # @test_nowarn ExactWrightFisher.Wright_Fisher_K_dim_exact_trajectory_arb(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec; use_progress_meter = true)
    # @test_nowarn ExactWrightFisher.Wright_Fisher_K_dim_exact_trajectory(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec; use_progress_meter = true)



    K = 4
    α = ones(K)* 1. /K
    sα = sum(α)
    Pop_size_WF3 = 15
    Ntimes_WF3 = 3
    time_step_WF3 = 0.1
    time_grid_WF3 = [k*time_step_WF3 for k in 0:(Ntimes_WF3-1)]
    Random.seed!(4)
    wfchain_WF3 = Wright_Fisher_K_dim_exact_trajectory(rand(Dirichlet(K,0.3)), time_grid_WF3[1:(end-1)], α)
    @test isreal(sum(wfchain_WF3))
end
