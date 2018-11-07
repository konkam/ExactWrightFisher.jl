using StatsFuns

@testset "helper functions" begin
    @test minus_1_power_i(5) == -1
    @test minus_1_power_i(2) == 1
    res = ExactWrightFisher.signed_log_sum_exp(log.(1:5), repeat([1],5))
    @test logsumexp(log.(1:5)) == res[2]
end
