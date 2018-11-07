println("Testing...")

using Test
using ExactWrightFisher
@test 1 == 1

include("test_helper_functions.jl")
include("test_A_infinity_functions.jl")
include("test_exact_wright_fisher_functions.jl")
