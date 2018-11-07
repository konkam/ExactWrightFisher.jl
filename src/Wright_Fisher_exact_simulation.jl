function Wright_Fisher_exact_transition(x::T, t::T, theta_1::T, theta_2::T) where T<:Real
  if t == 0
    return x
  else
    A∞ = Compute_A∞(theta_1 + theta_2, t)
    L = rand(Binomial(A∞, x))
    Y = rand(Beta(theta_1 + L, theta_2 + A∞ - L))
    return Y
  end
end

function Wright_Fisher_exact_trajectory(initial_value::T, times::AbstractVector{T}, theta_1::T, theta_2::T) where T<:Real
  trajectory = Array{Float64}(undef, length(times) + 1)
  times_with_0 = [0; times]
  trajectory[1] = initial_value
  for i in 2:(length(times)+1)
    trajectory[i] = Wright_Fisher_exact_transition(trajectory[i-1], times_with_0[i] - times_with_0[i-1], theta_1, theta_2)
  end
  return trajectory
end

function Wright_Fisher_K_dim_exact_transition(xvec::AbstractVector{T}, t::Real, αvec::AbstractVector{T}, sα::Real) where T<:Real
  # sα = sum(αvec)
  if t==0
    return xvec
  else
    A∞ = Compute_A∞(sα, t)
    L = rand(Multinomial(A∞, xvec))
    Y = rand(Dirichlet(L .+ αvec))
    return Y
  end
end
