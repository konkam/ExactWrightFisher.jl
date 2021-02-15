using StatsFuns

@testset "helper functions" begin
    @test ExactWrightFisher.minus_1_power_i(5) == -1
    @test ExactWrightFisher.minus_1_power_i(2) == 1
    res = ExactWrightFisher.signed_logsumexp(log.(1:5), repeat([1],5))
    @test res[1] == 1
    @test logsumexp(log.(1:5)) == res[2]

    res = ExactWrightFisher.signed_logsumexp([2, 4, 2], [1, -1, 1])
    @test res[1] == -1.0

    res = ExactWrightFisher.signed_logsumexp([2, 2, 4, 4], [1, -1, 1, -1])
    @test res[1] == 1.0
    @test res[2] == -Inf

    res = ExactWrightFisher.signed_logsumexp([log(3), log(3), log(4), log(4)], [1, -1, 1, -1])
    @test res[1] == 1.0
    @test res[2] == -Inf

    ref = ExactWrightFisher.signed_logsumexp_arb([0.2, -3., 1], [1, 1, -1])[2]
    @test ExactWrightFisher.signed_logsumexp(BigFloat[0.2, -3., 1], BigFloat[1, 1, -1])[2] ≈ Float64(ref)
    @test ExactWrightFisher.signed_logsumexp([0.2, -3., 1], [1, 1, -1])[2] ≈ Float64(ref)
    @test ExactWrightFisher.signed_logsumexp([0.2, -3., 1], [-1, -1, -1])[2] ≈ logsumexp([0.2, -3., 1])

    kvec = [377, 371, 363, 354, 345, 336, 326, 316, 306, 295, 284, 273, 262, 251, 240, 229, 218, 207, 195, 183, 171, 159, 147, 135, 123, 111, 100, 89, 78, 67, 57, 46, 35, 24, 14]
    logterms, signs = ExactWrightFisher.S_kvec_M_plus_logsum_pre_logsum(kvec, 0.05, 10)

    ref = ExactWrightFisher.signed_logsumexp_arb(logterms, signs)
    res = ExactWrightFisher.signed_logsumexp(logterms, signs)
    @test res[1] ≈ Float64(ref[1])
    @test res[2] ≈ Float64(ref[2])


    @test ExactWrightFisher.lgamma_local(1) == 0.0
    @test ExactWrightFisher.lgamma_local(2) == 0.0
    @test ExactWrightFisher.lgamma_local(3) == log(2)
end
