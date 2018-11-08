
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
    signed_log_sum_exp(lx, signs)

logsumexp with potentially subtraction of terms.

...
# Arguments
- `lx`: an array of logterms.
- `signs`: an array of signs for these terms.
...

# Examples
```julia-repl
julia>signed_log_sum_exp([0.2, -3., 1], [1, 1, -1])
2-element Array{Float64,1}:
 -1.0
  0.36955602678296184
```
"""
function signed_log_sum_exp(lx, signs)
  m = maximum(lx)
  scaled_sum = sum(signs .* exp.(lx .- m))
  if scaled_sum == 0
    return [1., -Inf]
  elseif scaled_sum < 0
    sgn = -1
    scaled_sum = -1*scaled_sum
  else
    sgn = 1
  end
  return [sgn, m + log(scaled_sum) ]#Will give an error if the sum is negative
end
