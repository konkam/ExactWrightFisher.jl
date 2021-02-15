using Distributions, SpecialFunctions


"""
      log of Eq. (5) of Jenkins, P. A., & Spano, D. (2017). Exact simulation of the Wright--Fisher diffusion. The Annals of Applied Probability, 27(3), 1478–1509.

``\$ a_{km}^θ = \frac{()θ+2k-1)(θ+m)_{(k-1)}}{m!(k-m)!}\$``
"""
function log_akmθ(θ::Real, k::Integer, m::Integer)
  if k < m
    error("in akmθ, m<=k must be satisfied")
  end
  if k == 0
    return 0
  else
    return log(θ+2*k-1) + lgamma_local(θ+m+k-1) - lgamma_local(θ+m) - logfactorial(m) - logfactorial(k-m)
  end
end

function log_bk_t_θ_t(k::Integer, t::Real, θ::Real, m::Integer)
  return log_akmθ(θ, k, m) - k*(k+θ-1)*t/2
end

function C_m_t_θ_rec(m::Integer, t::Real, θ::Real, log_bkm_current::Real, ind::Integer)
  log_b_im_p1 = log_bk_t_θ_t(m + ind + 1, t, θ, m)
  if log_bkm_current > log_b_im_p1
    return ind
  else
    return C_m_t_θ_rec(m, t, θ, log_b_im_p1, ind + 1)
  end
end

function C_m_t_θ(m::Integer, t::Real, θ::Real)
  return C_m_t_θ_rec(m, t, θ, log_bk_t_θ_t(m, t, θ, m), 0)
end

function S_kvec_M_plus_logsum_pre_logsum(kvec::Array{T, 1}, t::Real, θ::Real) where
T<:Integer
  M = length(kvec)

  two_kvec_plus_1 = sum(2*kvec .+ 1)

  U = typeof(t)

  logterms = Array{U}(undef, two_kvec_plus_1)
  signs = Array{Float64}(undef, two_kvec_plus_1)
  cnt = 1
  for m in 0:(M-1)
    for i in 0:(2*kvec[m+1])
    # for(int i = 0; i <= 2*kvec[m]; ++i) {
      logterms[cnt] = log_bk_t_θ_t(m+i, t, θ, m);
      signs[cnt] = minus_1_power_i(i);
      cnt += 1;
    end
  end
  return logterms, signs
end

function S_kvec_M_plus_logsum(kvec::Array{T, 1}, t::Real, θ::Real) where
T<:Integer
  logterms, signs = S_kvec_M_plus_logsum_pre_logsum(kvec, t, θ)
  return signed_logsumexp(logterms, signs);
end

function S_kvec_M_plus_logsum_arb(kvec::Array{T, 1}, t::Real, θ::Real) where
T<:Integer
  logterms, signs = S_kvec_M_plus_logsum_pre_logsum(kvec, t, θ)
  return signed_logsumexp_arb(logterms, signs)
end

S_kvec_M_plus_logsum_nosign(kvec::Array{T, 1}, t::Real, θ::Real) where
T<:Integer = S_kvec_M_plus_logsum(kvec, t, θ)[2]

function S_kvec_M_minus_log_newterms(kvec::Array{T, 1}, t::Real, θ::Real) where
T<:Integer

  M = length(kvec)
  logterms = Array{Float64}(undef, M)
  for m in 0:(M-1)
    logterms[m+1] = log_bk_t_θ_t(m + 2*kvec[m+1] + 1, t, θ, m)
  end
  return logterms
end


function S_kvec_M_both_logsumexp(kvec::Array{T, 1}, t::Real, θ::Real) where
T<:Integer
  logS_kvec_M_plus_res = S_kvec_M_plus_logsum(kvec, t, θ)
  sgn_logS_kvec_M_plus_res = logS_kvec_M_plus_res[1]
  sum_logS_kvec_M_plus_res = logS_kvec_M_plus_res[2]

  log_newterms = S_kvec_M_minus_log_newterms(kvec, t, θ)
  logsum_newterms =  signed_logsumexp(log_newterms, repeat([1.], length(log_newterms)))[2]

  # println(logS_kvec_M_plus_res)

  return S_kvec_M_both_logsumexp_inner(kvec, t, θ, logS_kvec_M_plus_res, log_newterms, logsum_newterms)
end

function S_kvec_M_both_logsumexp_arb(kvec::Array{T, 1}, t::Real, θ::Real) where
T<:Integer
  logS_kvec_M_plus_res = S_kvec_M_plus_logsum_arb(kvec, t, θ)
  sgn_logS_kvec_M_plus_res = logS_kvec_M_plus_res[1]
  sum_logS_kvec_M_plus_res = logS_kvec_M_plus_res[2]

  log_newterms = S_kvec_M_minus_log_newterms(kvec, t, θ)
  logsum_newterms =  signed_logsumexp_arb(log_newterms, repeat([1.], length(log_newterms)))[2]

  # println(logS_kvec_M_plus_res)

  return S_kvec_M_both_logsumexp_inner(kvec, t, θ, logS_kvec_M_plus_res, log_newterms, logsum_newterms)
end

function S_kvec_M_both_logsumexp_inner(kvec::Array{T, 1}, t::Real, θ::Real, logS_kvec_M_plus_res, log_newterms, logsum_newterms) where
T<:Integer
  sgn_logS_kvec_M_plus_res = logS_kvec_M_plus_res[1]
  sum_logS_kvec_M_plus_res = logS_kvec_M_plus_res[2]

  S_kvec_M_minus_res = 0.

  if sgn_logS_kvec_M_plus_res == -1
    if sum_logS_kvec_M_plus_res > logsum_newterms
      S_kvec_M_minus_res = -1 * exp(sum_logS_kvec_M_plus_res)*(1+exp(logsum_newterms-sum_logS_kvec_M_plus_res))
    else
      S_kvec_M_minus_res = -1 * exp(logsum_newterms)*(exp(sum_logS_kvec_M_plus_res-logsum_newterms) + 1)
    end
  else
    if sum_logS_kvec_M_plus_res > logsum_newterms
      S_kvec_M_minus_res = exp(sum_logS_kvec_M_plus_res)*(1-exp(logsum_newterms-sum_logS_kvec_M_plus_res))
    else
      S_kvec_M_minus_res = exp(logsum_newterms)*(exp(sum_logS_kvec_M_plus_res-logsum_newterms)-1)
    end
  end
  return [S_kvec_M_minus_res, exp(sum_logS_kvec_M_plus_res)]
end

function Compute_A∞_start0(θ::Real, t::Real; m::T = 0, kvec::Array{T,1} = [0]) where
T<:Integer
  U = rand(Uniform())
  return Compute_A∞_given_U(θ, t, U, m, kvec)
end

function Compute_A∞_start0_arb(θ::Real, t::Real; m::T = 0, kvec::Array{T,1} = [0]) where
T<:Integer
  U = rand(Uniform())
  return Compute_A∞_given_U_arb(θ, t, U, m, kvec)
end

function Compute_A∞_given_U(θ, t, U, m, kvec; debug = false)
  ### 0 indexing to stick with the article's notation
  while true
    kvec[m+1] = ceil(C_m_t_θ(m, t, θ)/2)
    if debug
      println("m=$m")
      # println("km=$kvec")
    end
    # kvec[m] = ceil(C_m_t_θ(m, t, θ)/2)
    # print(km)
    S_kvec_M_BOTH = S_kvec_M_both_logsumexp(kvec, t, θ)
    while (S_kvec_M_BOTH[1] < U) && (S_kvec_M_BOTH[2] > U)
      kvec = kvec .+ 1
      # print(kvec)
      S_kvec_M_BOTH = S_kvec_M_both_logsumexp(kvec, t, θ)
    end
    if S_kvec_M_BOTH[1] > U
      # print(paste('A∞ = ', m))
      return m
      # break()
    elseif (S_kvec_M_BOTH[2] < U)
      push!(kvec,0)
      m = m + 1
    end
    if debug
      println("params = $kvec, $t, $θ")
      println("S_kvec_M_BOTH=$S_kvec_M_BOTH")
    end
  end
end

function Compute_A∞_given_U_arb(θ, t, U, m, kvec; debug = false)
  ### 0 indexing to stick with the article's notation
  while true
    kvec[m+1] = ceil(C_m_t_θ(m, t, θ)/2)
    # kvec[m] = ceil(C_m_t_θ(m, t, θ)/2)
    if debug
      println("m=$m")
      # println("kvec=$kvec")
    end
    S_kvec_M_BOTH = S_kvec_M_both_logsumexp_arb(kvec, t, θ)
    while (S_kvec_M_BOTH[1] < U) && (S_kvec_M_BOTH[2] > U)
      kvec = kvec .+ 1
      # print(kvec)
      S_kvec_M_BOTH = S_kvec_M_both_logsumexp_arb(kvec, t, θ)
    end
    if S_kvec_M_BOTH[1] > U
      # print(paste('A∞ = ', m))
      return m
      # break()
    elseif (S_kvec_M_BOTH[2] < U)
      push!(kvec,0)
      m = m + 1
    end
    if debug
      println("params = $kvec, $t, $θ")
      println("S_kvec_M_BOTH=$S_kvec_M_BOTH")
    end
  end
end


function β(θ::Real, t::Real)
  if θ == 0
    return 0
  else
    1/2*(θ-1)*t
  end
end


function η_of_β(β::Real)
  if β == 1
    return 1
  else
    return β/(exp(β)-1)
  end
end

η(θ::Real, t::Real) = η_of_β(β(θ, t))

μ(t::Real, θ::Real) = 2*η(θ, t)/t

function Compute_A∞_good_m_start(θ::Real, t::Real)
  m0::Int64 = μ(t, θ) |> round
  # println(m0)
  kvec0::Array{Int64,1} = ceil.(C_m_t_θ.(0:m0, t, θ) ./ 2)
  # println(kvec0)
  # println(length(kvec0))
  return Compute_A∞_start0(θ, t, m = m0, kvec = kvec0)
end

function Compute_A∞_good_m_start_arb(θ::Real, t::Real)
  m0::Int64 = μ(t, θ) |> round
  # println(m0)
  kvec0::Array{Int64,1} = ceil.(C_m_t_θ.(0:m0, t, θ) ./ 2)
  # println(kvec0)
  # println(length(kvec0))
  return Compute_A∞_start0_arb(θ, t, m = m0, kvec = kvec0)
end

"""
      Return a sample from the ancestral process A∞(t) of Kingman’s co-
alescent with mutation.

At the moment, this function starts from a rough guess to increase the acceptance rate of the algorithm. The guess is based on a small time step approximation of the transition function.

# Examples
```julia-repl
julia> ExactWrightFisher.Compute_A∞(1.5, 1.)
2
```
"""
Compute_A∞(θ::Real, t::Real) = Compute_A∞_good_m_start(θ, t)
Compute_A∞_arb(θ::Real, t::Real) = Compute_A∞_good_m_start_arb(θ, t)
