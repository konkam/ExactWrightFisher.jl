using ProgressMeter

function Wright_Fisher_exact_transition(x::Real, t::Real, theta_1::Real, theta_2::Real)
  if t == 0
    return x
  else
    A∞ = Compute_A∞(theta_1 + theta_2, t)
    L = rand(Binomial(A∞, x))
    Y = rand(Beta(theta_1 + L, theta_2 + A∞ - L))
    return Y
  end
end

function cmp_1D_trajectory(initial_value::Real, times::AbstractVector{T}, transition_function::Function; use_progress_meter = false) where T<:Real
  trajectory = Array{Float64}(undef, length(times) + 1)
  times_with_0 = [0; times]
  trajectory[1] = initial_value
  imax = length(times)+1

  if use_progress_meter
    # println("using the progress meter")
    @showprogress 3 "Computing trajectory..." for i in 2:imax
    p = Progress(i, 2, "Computing initial pass...", imax)
    println("$(i-1) out of $(imax-1) iterations")
      trajectory[i] = transition_function(trajectory[i-1], times_with_0[i] - times_with_0[i-1])
    end
  else
    for i in 2:imax
      trajectory[i] = transition_function(trajectory[i-1], times_with_0[i] - times_with_0[i-1])
    end
  end
  return trajectory
end

function Wright_Fisher_exact_trajectory(initial_value::Real, times::AbstractVector{T}, theta_1::Real, theta_2::Real; use_progress_meter = false) where T<:Real

  function  WF_1D_transition_fun(st::Real, dt::Real)::Real
    return Wright_Fisher_exact_transition(st, dt, theta_1, theta_2)
  end

  return cmp_1D_trajectory(initial_value, times, WF_1D_transition_fun; use_progress_meter = use_progress_meter)
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

function cmp_K_dim_trajectory(initial_vec::AbstractVector{T}, times::AbstractVector{T}, transition_function::Function; use_progress_meter = false) where T<:Real
  trajectory = Array{Float64}(undef, length(initial_vec), length(times) + 1)
  times_with_0 = [0; times]
  trajectory[:,1] = initial_vec
  imax = length(times)+1

  if use_progress_meter
    @showprogress 3 "Computing trajectory..." for i in 2:imax
    println("$(i-1) out of $(imax-1) iterations")
      trajectory[:,i] = transition_function(trajectory[:,i-1], times_with_0[i] - times_with_0[i-1])
    end
  else
    for i in 2:imax
      trajectory[:,i] = transition_function(trajectory[:,i-1], times_with_0[i] - times_with_0[i-1])
    end
  end
  return trajectory
end

function Wright_Fisher_K_dim_exact_trajectory(initial_vec::AbstractVector{T}, times::AbstractVector{T}, αvec::AbstractVector{T}; use_progress_meter = false) where T<:Real
  sα = sum(αvec)

  function WF_exact_transition_function(st::AbstractVector{T}, dt::Real)::Vector{Float64}  where T<:Real
    return Wright_Fisher_K_dim_exact_transition(st, dt, αvec, sα)
  end

  return cmp_K_dim_trajectory(initial_vec, times, WF_exact_transition_function; use_progress_meter = use_progress_meter)
end
