using ProgressMeter

"""
      Wright_Fisher_exact_transition(x::Real, t::Real, θ_1::Real, θ_2::Real)

Sample exactly from the 1-D Wright-Fisher transition function. Function may hang for time steps <= 0.05
...
# Arguments
- `x`: the starting point.
- `t`: the time step.
- `θ_1, θ_2`: parameters of the 1-D Wright-Fisher process.
...

# Examples
```julia-repl
julia> ExactWrightFisher.Wright_Fisher_exact_transition(0.5, 0.54, 0.75/2, 0.75/2)
0.9177951585792223
```
"""
function Wright_Fisher_exact_transition(x::Real, t::Real, θ_1::Real, θ_2::Real)
  if t == 0
    return x
  else
    A∞ = Compute_A∞(θ_1 + θ_2, t)
    L = rand(Binomial(A∞, x))
    Y = rand(Beta(θ_1 + L, θ_2 + A∞ - L))
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

"""
      Wright_Fisher_exact_trajectory(initial_value::Real, times::AbstractVector{T}, θ_1::Real, θ_2::Real; use_progress_meter = false) where T<:Real

Sample exactly a 1-D Wright-Fisher trajectory. Function may hang for time steps <= 0.05
...
# Arguments
- `initial_value`: the starting point.
- `times`: the observation times.
- `θ_1, θ_2`: parameters of the 1-D Wright-Fisher process.
...

# Examples
```julia-repl
julia>  ExactWrightFisher.Wright_Fisher_exact_trajectory(0.5, range(0.1, stop = 1, length = 5), 0.75/2, 0.75/2)
[0.5, 0.7391267869778926, 0.9850784214783245, 0.8237514810336812, 0.8039797099052189, 0.6671481291268048]
```
"""
function Wright_Fisher_exact_trajectory(initial_value::Real, times::AbstractVector{T}, θ_1::Real, θ_2::Real; use_progress_meter = false) where T<:Real

  function  WF_1D_transition_fun(st::Real, dt::Real)::Real
    return Wright_Fisher_exact_transition(st, dt, θ_1, θ_2)
  end

  return cmp_1D_trajectory(initial_value, times, WF_1D_transition_fun; use_progress_meter = use_progress_meter)
end

"""
      Wright_Fisher_K_dim_exact_transition(xvec::AbstractVector{T}, t::Real, αvec::AbstractVector{T}, sα::Real) where T<:Real

Sample exactly from the K dimensional Wright-Fisher transition function. Function may hang for time steps <= 0.05
...
# Arguments
- `xvec`: the starting state.
- `t`: the time step.
- `αvec`: parameters of the K dimensional Wright-Fisher process.
...

# Examples
```julia-repl
julia> ExactWrightFisher.Wright_Fisher_K_dim_exact_transition(rand(Dirichlet(α_vec |> collect)), 0.5, α_vec)
[0.02776997207332505, 0.2336806438432984, 0.07839887060526822, 0.6601505134781083]
```
"""
Wright_Fisher_K_dim_exact_transition(xvec::AbstractVector{T}, t::Real, αvec::AbstractVector{T}) where T<:Real = Wright_Fisher_K_dim_exact_transition(xvec, t, αvec, sum(αvec))

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

"""
      Wright_Fisher_K_dim_exact_trajectory(initial_vec::AbstractVector{T}, times::AbstractVector{T}, αvec::AbstractVector{T}; use_progress_meter = false) where T<:Real

Sample exactly a K dimensional Wright-Fisher trajectory. Function may hang for time steps <= 0.05
...
# Arguments
- `initial_vec`: the starting state.
- `times`: the observation times.
- `αvec`: parameters of the K dimensional Wright-Fisher process.
...

# Examples
```julia-repl
julia> ExactWrightFisher.Wright_Fisher_K_dim_exact_trajectory(rand(Dirichlet(α_vec |> collect)), range(0, stop = 1, length = 10), α_vec)
4×11 Array{Float64,2}:
 1.06585e-5   1.06585e-5   1.08036e-5  …  0.000108104  0.00333635  0.022469
 0.000657425  0.000657425  0.00340315     0.0150088    0.0103089   7.54402e-6
 0.662003     0.662003     0.837201       0.978252     0.976614    0.943636
 0.337329     0.337329     0.159385       0.00663071   0.00974101  0.0338876
```
"""
function Wright_Fisher_K_dim_exact_trajectory(initial_vec::AbstractVector{T}, times::AbstractVector{T}, αvec::AbstractVector{T}; use_progress_meter = false) where T<:Real
  sα = sum(αvec)

  function WF_exact_transition_function(st::AbstractVector{T}, dt::Real)::Vector{Float64}  where T<:Real
    return Wright_Fisher_K_dim_exact_transition(st, dt, αvec, sα)
  end

  return cmp_K_dim_trajectory(initial_vec, times, WF_exact_transition_function; use_progress_meter = use_progress_meter)
end
