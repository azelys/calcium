"""Parameters for the Komarova et al. (2003) well-mixed bone-cell model.

The model uses the published power-law form:

    dx₁/dt = α₁ x₁^g₁₁ x₂^g₂₁ − β₁ x₁
    dx₂/dt = α₂ x₁^g₁₂ x₂^g₂₂ − β₂ x₂
    dz/dt  = −k₁ y₁ + k₂ y₂

where x₁ and x₂ are osteoclast and osteoblast populations. The y terms are
the excess populations above their non-trivial steady state, as in the
published bone-mass equation.
"""
Base.@kwdef struct KomarovaParams{T<:Real}
    α₁::T = 3.0
    α₂::T = 4.0
    β₁::T = 0.2
    β₂::T = 0.02
    g₁₁::T = 1.1
    g₂₁::T = -0.5
    g₁₂::T = 1.0
    g₂₂::T = 0.0
    k₁::T = 0.093
    k₂::T = 0.0008
    z₀::T = 100.0
end

"""Return the non-trivial steady-state cell populations."""
function komarova_equilibrium(p::KomarovaParams)
    γ = p.g₁₂ * p.g₂₁ - (1 - p.g₁₁) * (1 - p.g₂₂)
    γ ≈ 0 && throw(ArgumentError("singular Komarova exponent matrix (γ ≈ 0)"))
    x̄₁ = (p.β₁ / p.α₁)^((1 - p.g₂₂) / γ) *
          (p.β₂ / p.α₂)^(p.g₂₁ / γ)
    x̄₂ = (p.β₁ / p.α₁)^(p.g₁₂ / γ) *
          (p.β₂ / p.α₂)^((1 - p.g₁₁) / γ)
    return (x̄₁=x̄₁, x̄₂=x̄₂)
end

@inline excess(x, x̄) = max(x - x̄, zero(x))

function komarova_rhs!(du, u, p::KomarovaParams, t)
    x₁, x₂, _ = u
    x₁ <= 0 && throw(DomainError(x₁, "osteoclast population left the positive domain"))
    x₂ <= 0 && throw(DomainError(x₂, "osteoblast population left the positive domain"))

    du[1] = p.α₁ * x₁^p.g₁₁ * x₂^p.g₂₁ - p.β₁ * x₁
    du[2] = p.α₂ * x₁^p.g₁₂ * x₂^p.g₂₂ - p.β₂ * x₂

    eq = komarova_equilibrium(p)
    y₁ = excess(x₁, eq.x̄₁)
    y₂ = excess(x₂, eq.x̄₂)
    du[3] = -p.k₁ * y₁ + p.k₂ * y₂
    return nothing
end

function komarova_problem(p::KomarovaParams=KomarovaParams();
                          u₀=nothing, tspan=(0.0, 100.0))
    eq = komarova_equilibrium(p)
    initial = isnothing(u₀) ? [eq.x̄₁, eq.x̄₂, p.z₀] : collect(u₀)
    length(initial) == 3 || throw(ArgumentError("u₀ must have three states"))
    return ODEProblem(komarova_rhs!, initial, tspan, p)
end

"""Solve the Komarova model with a reproducible default ODE solver."""
function simulate_komarova(p::KomarovaParams=KomarovaParams(); kwargs...)
    sol = solve(komarova_problem(p; kwargs...), Tsit5();
                abstol=1e-10, reltol=1e-10, saveat=0.1)
    SciMLBase.successful_retcode(sol) || error("Komarova solve failed: $(sol.retcode)")
    return sol
end
