using StatsFuns
"""
    minus_1_power_i(i::Int64)

Returns \$(-1)^i\$ with hopefully the least amount of work necessary.

# Examples
```julia-repl
julia> minus_1_power_i(1)
-1
```
"""
function minus_1_power_i(i::Int64)
    if iseven(i)
        return 1
    else
        return -1
    end
end

"""
    signed_logsumexp(lx, signs)

logsumexp with potentially subtraction of terms.

...
# Arguments
- `lx`: an array of logterms.
- `signs`: an array of signs for these terms.
...

# Examples
```julia-repl
julia>signed_logsumexp([0.2, -3., 1], [1, 1, -1])
2-element Array{Float64,1}:
 -1.0
  0.36955602678296184
```
"""
# function signed_logsumexp(lx, signs)
## There might potentially be a problem when the larger terms cancel each other
#   m = maximum(lx)
#   scaled_sum = sum(signs .* exp.(lx .- m))
#   if abs(scaled_sum) <= 10*eps(Float64)
#     return [1., -Inf]
#   elseif scaled_sum < 0
#     sgn = -1
#     scaled_sum = -1*scaled_sum
#   else
#     sgn = 1
#   end
#   return [sgn, m + log(scaled_sum) ]#Will give an error if the sum is negative
# end

function signed_logsumexp(lx, signs)
  # summing the positive terms together and the negative terms together should decrease the probability of cancellation of large terms
  # @assert length(lx) == length(signs)
  # pos = signs .== 1
  if all(signs .> 0)
    return [1.0, logsumexp(lx)]
  elseif all(signs .< 0)
    return [-1.0, logsumexp(lx)]
  else
    n = length(lx)
    @inbounds logsumexp_positive_terms = logsumexp(lx[i] for i in 1:n if signs[i] > 0)
    @inbounds logsumexp_negative_terms = logsumexp(lx[i] for i in 1:n if signs[i] < 0)
    if logsumexp_positive_terms > logsumexp_negative_terms
      sgn = 1
      res = log(1-exp(logsumexp_negative_terms - logsumexp_positive_terms)) + logsumexp_positive_terms
    elseif logsumexp_positive_terms < logsumexp_negative_terms
      sgn = -1
      res = log(1-exp(logsumexp_positive_terms - logsumexp_negative_terms)) + logsumexp_negative_terms
    else
      sgn = 1
      res = -Inf
    end
    return [sgn, res]
  end
end

import Base.sign

function sign(x::arb)
    if x < 0
        return -1
    else
        return 1
    end
end
function signed_logsumexp_arb(lx, signs)
    res = sum(exp.(RR.(lx)) .* signs)
    return sign(res), log(abs(res))
end
