"""Lemaire et al. (2004) osteoblast–osteoclast population model.

This implementation follows the equations and parameter values in the
Physiome CellML translation of the model. The optional `ob_removal_rate`
reproduces the CellML extension that removes active osteoblasts between
days 20 and 80; its default of zero is the core model.
"""
Base.@kwdef struct LemaireParams{T<:Real}
    D_R::T = 7.0e-4
    D_B::T = 0.05 * 0.7
    D_C::T = 2.1e-3
    D_A::T = 0.7
    k_B::T = 0.189
    f₀::T = 0.05
    C_s::T = 5.0e-3
    k₁::T = 1.0e-2
    k₂::T = 10.0
    k₃::T = 5.8e-4
    k₄::T = 1.7e-2
    K::T = 10.0
    k_o::T = 0.35
    I_o::T = 0.0
    I_L::T = 0.0
    r_L::T = 1.0e3
    K_OP::T = 2.0e5
    K_LP::T = 3.0e6
    I_P::T = 0.0
    k_P::T = 86.0
    S_P::T = 250.0
    k₅::T = 0.02
    k₆::T = 3.0
    ob_removal_rate::T = 0.0
    ob_removal_start::T = 20.0
    ob_removal_stop::T = 80.0
end

@inline function lemaire_aux(u, p::LemaireParams)
    R, B, C = u
    pi_C = (C + p.f₀ * p.C_s) / (C + p.C_s)
    P = p.I_P / p.k_P
    P₀ = p.S_P / p.k_P
    P_s = p.k₆ / p.k₅
    pi_P = (P + P₀) / (P + P_s)
    pi_L = (p.k₃ / p.k₄) *
           ((p.K_LP * pi_P * B) /
            (1 + (p.k₃ * p.K / p.k₄) +
             (p.k₁ / (p.k₂ * p.k_o)) * ((p.K_OP / pi_P) * R + p.I_o))) /
           (1 + p.I_L / p.r_L)
    return (pi_C=pi_C, pi_P=pi_P, pi_L=pi_L)
end

function lemaire_rhs!(du, u, p::LemaireParams, t)
    R, B, C = u
    R <= 0 && throw(DomainError(R, "responding osteoblast population left the positive domain"))
    B <= 0 && throw(DomainError(B, "active osteoblast population left the positive domain"))
    C <= 0 && throw(DomainError(C, "active osteoclast population left the positive domain"))

    a = lemaire_aux(u, p)
    f = (p.ob_removal_rate != 0 &&
         p.ob_removal_start < t <= p.ob_removal_stop) ? -p.ob_removal_rate : 0.0

    du[1] = p.D_R * a.pi_C - (p.D_B / a.pi_C) * R
    du[2] = (p.D_B / a.pi_C) * R - p.k_B * B + f
    du[3] = p.D_C * a.pi_L - p.D_A * a.pi_C * C
    return nothing
end

function lemaire_problem(p::LemaireParams=LemaireParams();
                         u₀=[7.7e-4, 7.3e-4, 9.1e-4], tspan=(0.0, 140.0))
    length(u₀) == 3 || throw(ArgumentError("u₀ must have three states"))
    return ODEProblem(lemaire_rhs!, collect(u₀), tspan, p)
end

"""Solve the Lemaire model with a reproducible default ODE solver."""
function simulate_lemaire(p::LemaireParams=LemaireParams(); kwargs...)
    sol = solve(lemaire_problem(p; kwargs...), Tsit5();
                abstol=1e-12, reltol=1e-10, saveat=0.1)
    SciMLBase.successful_retcode(sol) || error("Lemaire solve failed: $(sol.retcode)")
    return sol
end
