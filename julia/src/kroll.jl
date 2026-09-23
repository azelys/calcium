"""Kroll (2000) temporal PTH model.

This implementation follows the public Physiome CellML translation. It is a
five-state ODE artifact with preosteoblasts (`Y`), osteoblasts (`X`), PTH
(`P`), osteoclasts (`Z`), and IL-6. The original paper motivates explicit
transit delays; the published CellML artifact encodes the reported PTH pulse
experiment as an ODE schedule, so this file preserves that artifact and makes
the schedule configurable.

Source artifact:
https://models.physiomeproject.org/workspace/kroll_2000/rawfile/5686223dbf499c1738659fde4c45c9eddbfd0dc2/kroll_2000.cellml
"""
Base.@kwdef struct KrollParams{T<:Real}
    C::T = 50.0
    k₁::T = 1.0
    k₂::T = 1.3
    k₃::T = 0.05
    k₄::T = 0.9
    k₅::T = 5.0
    k₆::T = 0.02
    kᵧ::T = 0.01
    K::T = 5.0
    K₂::T = 2.0
    C₁::T = 1.0
    C₂::T = 1.0
    C₃::T = 1.0
    input_rate::T = 10.0
    pulse_period::T = 12.0
    pulse_duration::T = 6.0
    pulse_start::T = 0.0
    mode::Symbol = :intermittent
end

@inline function kroll_pulse_on(t, p::KrollParams)
    p.mode == :continuous && return true
    p.mode == :intermittent ||
        throw(ArgumentError("Kroll mode must be :continuous or :intermittent"))
    p.pulse_period > 0 || throw(ArgumentError("pulse_period must be positive"))
    0 <= p.pulse_duration <= p.pulse_period ||
        throw(ArgumentError("pulse_duration must lie in [0, pulse_period]"))
    return t >= p.pulse_start &&
           mod(t - p.pulse_start, p.pulse_period) < p.pulse_duration
end

function kroll_aux(u, p::KrollParams, t)
    Y, X, P, Z, IL6 = u
    pth_activation = P / (p.K + P)
    return (; Y, X, P, Z, IL6, pth_activation,
            pth_input=(kroll_pulse_on(t, p) ? p.input_rate : zero(P)),
            osteoblast_to_osteoclast=(X / Z))
end

function kroll_rhs!(du, u, p::KrollParams, t)
    Y, X, P, Z, IL6 = u
    any(x -> x <= 0, (Y, X, P, Z)) &&
        throw(DomainError(u, "Kroll cell/PTH states must remain positive"))
    IL6 < 0 && throw(DomainError(IL6, "Kroll IL-6 state must remain nonnegative"))

    q = P / (p.K + P)
    input = kroll_pulse_on(t, p) ? p.input_rate : zero(P)
    du[1] = p.k₁ * p.C₁ * q * p.C -
            (p.k₂ * p.C₂ * (1 - q) + p.kᵧ) * Y
    du[2] = p.k₂ * p.C₂ * (1 - q) * Y - p.k₃ * X
    du[3] = input - p.k₄ * P
    du[4] = p.k₅ * p.C₃ * (IL6 / (p.K₂ + IL6)) - p.k₆ * Z
    du[5] = 0.1 * X - 10.0 * IL6
    return nothing
end

function kroll_problem(p::KrollParams=KrollParams();
                       u₀=[10.0, 500.0, 10.0, 200.0, 0.0],
                       tspan=(0.0, 96.0))
    length(u₀) == 5 || throw(ArgumentError("u₀ must have five states"))
    return ODEProblem(kroll_rhs!, collect(u₀), tspan, p)
end

function kroll_tstops(p::KrollParams, tspan)
    p.mode == :continuous && return Float64[]
    t₀, t₁ = tspan
    n₀ = floor(Int, (t₀ - p.pulse_start) / p.pulse_period) - 1
    n₁ = ceil(Int, (t₁ - p.pulse_start) / p.pulse_period) + 1
    candidates = Float64[]
    for n in n₀:n₁
        push!(candidates, p.pulse_start + n * p.pulse_period)
        push!(candidates, p.pulse_start + n * p.pulse_period + p.pulse_duration)
    end
    return sort(unique(filter(t -> t₀ < t < t₁, candidates)))
end

"""Solve the Kroll CellML translation with a reproducible default solver."""
function simulate_kroll(p::KrollParams=KrollParams(); kwargs...)
    tspan = haskey(kwargs, :tspan) ? kwargs[:tspan] : (0.0, 96.0)
    sol = solve(kroll_problem(p; kwargs...), Tsit5();
                abstol=1e-10, reltol=1e-10, saveat=0.05,
                tstops=kroll_tstops(p, tspan), dtmax=0.1)
    SciMLBase.successful_retcode(sol) || error("Kroll solve failed: $(sol.retcode)")
    return sol
end
