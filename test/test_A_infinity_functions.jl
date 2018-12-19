using Random, Distributions

@testset "testing Ancestor(A∞) functions" begin
    @test ExactWrightFisher.log_akmθ(0.75, 6, 4) ≈ 8.026285751157305 atol=10^(-8)
    @test ExactWrightFisher.log_bk_t_θ_t(8, 0.2, 0.75, 6) ≈ 5.072288573098844 atol=10^(-8)
    @test ExactWrightFisher.C_m_t_θ(5,0.2,0.75) == 3
    ref = ExactWrightFisher.S_kvec_M_plus_logsum_arb([1,3,50, 10], 0.0002, 0.075)
    res = ExactWrightFisher.S_kvec_M_plus_logsum([1,3,50, 10], 0.0002, 0.075)
    @test res[1] == ref[1]
    @test res[2] ≈ Float64(ref[2])

    ref = ExactWrightFisher.S_kvec_M_both_logsumexp_arb([1,3,5], 0.2, 0.75)
    res = ExactWrightFisher.S_kvec_M_both_logsumexp([1,3,5], 0.2, 0.75)
    @test res[1] ≈ Float64(ref[1])
    @test res[2] ≈ Float64(ref[2])

    ref = ExactWrightFisher.S_kvec_M_both_logsumexp_arb([1,3,50], 0.2, 0.75)
    res = ExactWrightFisher.S_kvec_M_both_logsumexp([1,3,50], 0.2, 0.75)
    @test res[1] ≈ Float64(ref[1])
    @test res[2] ≈ Float64(ref[2])

    ref = ExactWrightFisher.S_kvec_M_both_logsumexp_arb([1,3,50, 10], 0.0002, 0.075)
    res = ExactWrightFisher.S_kvec_M_both_logsumexp([1,3,50, 10], 0.0002, 0.075)
    @test res[1] ≈ Float64(ref[1])
    @test res[2] ≈ Float64(ref[2])

    ref = ExactWrightFisher.S_kvec_M_both_logsumexp_arb([1,3,50, 10], 0.0002, 75.)
    res = ExactWrightFisher.S_kvec_M_both_logsumexp([1,3,50, 10], 0.0002, 75.)
    @test res[1] ≈ Float64(ref[1])
    @test res[2] ≈ Float64(ref[2])

    @test ExactWrightFisher.S_kvec_M_plus_logsum_nosign([1,3,5], 0.2, 0.75) ≈ -0.6936149059655459 atol=10^(-8)
    @test ExactWrightFisher.S_kvec_M_plus_logsum_nosign([1,3,50,100], 0.2, 0.75) ≈ -0.6945637479851614 atol=10^(-8)
    @test ExactWrightFisher.S_kvec_M_plus_logsum_nosign([1,3,50,10], 0.0002, 75.) ≈ 125.40882239903411 atol=10^(-8)
    Random.seed!(0)
    @test ExactWrightFisher.Compute_A∞(1.5, 1.) == 2
    @test ExactWrightFisher.β(1.5,0.05) == 0.0125
    @test ExactWrightFisher.β(0.,0.05) == 0.
    @test ExactWrightFisher.η_of_β(0.02) ≈ 0.9900333331111149 atol=10^(-8)
    @test ExactWrightFisher.η_of_β(1.) == 1.
    @test ExactWrightFisher.η(1.5, 0.05) ≈ 0.9937630207994219 atol=10^(-8)
    @test ExactWrightFisher.μ(1.5, 0.05) ≈ 1.864268029836541 atol=10^(-8)

    ref = ExactWrightFisher.Compute_A∞_given_U_arb(sum(1:4), 0.05, 0.6, 0, [0])
    res = ExactWrightFisher.Compute_A∞_given_U(sum(1:4), 0.05, 0.6, 0, [0])
    @test ref == res
    @test ref == 36

    Random.seed!(0);
    U = rand(Uniform())
    ref = ExactWrightFisher.Compute_A∞_given_U_arb(sum(1:4), 0.05, U, 0, [0])
    res = ExactWrightFisher.Compute_A∞_given_U(sum(1:4), 0.05, U, 0, [0])

    Random.seed!(0)
    @test ExactWrightFisher.Compute_A∞(sum((1:4)), 0.05) == 37  ###This is because there is floating point inaccuracies
    Random.seed!(0)
    @test ExactWrightFisher.Compute_A∞_arb(sum((1:4)), 0.05) == 38
end
