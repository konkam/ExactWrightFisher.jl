function σ(t::Real, θ::Real)
  if θ == 1
    return  2/(3*t)
  else
    β_θ_t = β(θ, t)
    η_θ_t = η_of_β(β_θ_t)
    return 2*η_θ_t/t*(η_θ_t+β_θ_t)^2*(1+η_θ_t/(η_θ_t+β_θ_t)-2*η_θ_t)/β_θ_t^2
  end
end

"""
      Normal approximation of the transition function, valid for small time steps. In practice, used when t<=0.05. This is Theorem 1 of Jenkins, P. A., & Spano, D. (2017). Exact simulation of the Wright--Fisher diffusion. The Annals of Applied Probability, 27(3), 1478–1509.

"""
function Compute_A∞_approx(θ::Real, t::Real)
  A∞_real = rand(Normal(μ(t, θ), σ(t, θ)))
  if A∞_real<=0
    return 0
  else
    return Int64(round(A∞_real))
  end
end
