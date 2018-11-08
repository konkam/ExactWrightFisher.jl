using StatsFuns

@testset "helper functions" begin
    @test ExactWrightFisher.minus_1_power_i(5) == -1
    @test ExactWrightFisher.minus_1_power_i(2) == 1
    res = ExactWrightFisher.signed_log_sum_exp(log.(1:5), repeat([1],5))
    @test res[1] == 1
    @test logsumexp(log.(1:5)) == res[2]
    res = ExactWrightFisher.signed_log_sum_exp([2, 4, 2], [1, -1, 1])
    @test res[1] == -1.0
    res = ExactWrightFisher.signed_log_sum_exp([2, 2, 4, 4], [1, -1, 1, -1])
    @test res[1] == 1.0
    @test res[2] == 0
end
