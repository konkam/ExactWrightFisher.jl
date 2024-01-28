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
Compute_A∞_approx(θ::Real, t::Real)

Normal approximation of the transition function, valid for small time steps. In practice, used when t<=0.05. This is Theorem 1 of Jenkins, P. A., & Spano, D. (2017). Exact simulation of the Wright--Fisher diffusion. The Annals of Applied Probability, 27(3), 1478–1509.

Inputs:
- θ::Real: The parameter θ of the transition function.
- t::Real: The time step.

Outputs:
- An integer approximation of A∞
"""

function Compute_A∞_approx(θ::Real, t::Real)
    A∞_real = round(rand(Normal(μ(t, θ), σ(t, θ))))
    return max(0, Int64(A∞_real))
end
