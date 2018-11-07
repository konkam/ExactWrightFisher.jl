function Wright_Fisher_K_dim_with_t005_approx(xvec::AbstractVector{T}, t::Real, αvec::AbstractVector{T}, sα::Real) where T<:Real
  # sα = sum(αvec)
  if t==0
    return xvec
  elseif t<=0.05
    A∞ = Compute_A∞_approx(sα, t)
  else
    A∞ = Compute_A∞(sα, t)
  end
  L = rand(Multinomial(A∞, xvec))
  Y = rand(Dirichlet(L .+ αvec))
  return Y
end

function Wright_Fisher_K_dim_trajectory_with_t005_approx(initial_vec::AbstractVector{T}, times::AbstractVector{T}, αvec::AbstractVector{T}) where T<:Real
  sα = sum(αvec)

  function WF_exact_transition_function(st::AbstractVector{T}, dt::Real)::Vector{Float64}  where T<:Real
    return Wright_Fisher_K_dim_with_t005_approx(st, dt, αvec, sα)
  end

  return cmp_K_dim_trajectory(initial_vec, times, WF_exact_transition_function)
end
