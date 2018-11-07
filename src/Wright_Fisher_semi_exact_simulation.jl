function Wright_Fisher_exact_transition_with_t005_approx(x::Real, t::Real, theta_1::Real, theta_2::Real)
  if t == 0
    return x
  elseif t<=0.05
    A∞ = Compute_A∞_approx(theta_1 + theta_2, t)
  else
    A∞ = Compute_A∞(theta_1 + theta_2, t)
  end
  L = rand(Binomial(A∞, x))
  Y = rand(Beta(theta_1 + L, theta_2 + A∞ - L))
  return Y
end

function Wright_Fisher_exact_trajectory_with_t005_approx(initial_value::Real, times::AbstractVector{T}, theta_1::Real, theta_2::Real; use_progress_meter = false) where T<:Real

  function  WF_1D_transition_fun(st::Real, dt::Real)::Real
    return Wright_Fisher_exact_transition_with_t005_approx(st, dt, theta_1, theta_2)
  end

  return cmp_1D_trajectory(initial_value, times, WF_1D_transition_fun; use_progress_meter = use_progress_meter)
end

function Wright_Fisher_K_dim_transition_with_t005_approx(xvec::AbstractVector{T}, t::Real, αvec::AbstractVector{T}, sα::Real) where T<:Real
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
    return Wright_Fisher_K_dim_transition_with_t005_approx(st, dt, αvec, sα)
  end

  return cmp_K_dim_trajectory(initial_vec, times, WF_exact_transition_function)
end
