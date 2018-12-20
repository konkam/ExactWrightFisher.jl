println("Testing...")

using Test
using ExactWrightFisher
@test 1 == 1

@time include("test_helper_functions.jl")
@time include("test_A_infinity_functions.jl")
@time include("test_exact_wright_fisher_functions.jl")
# @time include("test_semi_exact_wright_fisher_functions.jl")
