"""
Nelson et al. (2026) surgical-menopause bone-remodelling model.

This is a Julia translation of the public MATLAB repository
`ashleefv/SurgicalMenopauseBone`.  It is a seven-state ODE for pre-
osteoblasts, pre-osteoclasts, osteoclasts, osteoblasts, osteocytes,
sclerostin, and normalized bone density.  Estrogen is an imposed input;
the optional post-surgical effects modify osteocyte apoptosis and
pre-osteoclast differentiation.

The source artifact uses days as its time unit and computes the initial
cell-state vector by solving the six cell equations at steady state with
relative estrogen equal to one.  The Julia API follows that convention.
"""

const NELSON_STATE_NAMES =
    (:PB, :PC, :C, :B, :S, :Sc, :Bd)

struct NelsonParams{T<:NamedTuple}
    values::T
end

Base.getproperty(p::NelsonParams, s::Symbol) =
    s === :values ? getfield(p, :values) : getproperty(getfield(p, :values), s)

const _NELSON_BASE = (
    BMC_0 = 0.8,
    e_C = 0.990728213,
    e_PC = 0.937654795,
    e_S = 9.595406602,
    eta_B = 0.00867806,
    eta_C = 0.023815993,
    eta_S = 0.000109589,
    e_Sc = 9.595406602,
    kappa_Sc = 0.05,
    lambda_B = 1.29e-6,
    lambda_C = 3.82e-6,
    n = 1.0,
    nu_Omega = 107.5332737,
    omega_B = 0.000624353,
    omega_PB = 0.319241069,
    omega_PC = 0.930890444,
    r_C = 10.12488334,
    r_Omega = 1024.4198,
    sc_Omega = 3039.645954,
    sc_PB = 163.2824648,
    sc_PC = 8603687.044,
    t_m_years = 50.0,
    tau_E_years = 2.6,
    t_m = 9855.0,
    tau_E = 949.0,
    E_0 = 156.0,
    Eovx = (0.065 * (24 * 60) / 156.0) / (log(2) * 24 * 60 / 161),
    kappa_E = log(2) * 24 * 60 / 161,
    k_syn = 0.065 * (24 * 60) / 156.0,
    eta_ovx = 0.0,
    tau = 0.0,
    omega_ovx = 0.0,
)

function NelsonParams(; kwargs...)
    unknown = setdiff(keys(kwargs), keys(_NELSON_BASE))
    isempty(unknown) || throw(ArgumentError("unknown Nelson parameter(s): $(collect(unknown))"))
    return NelsonParams(merge(_NELSON_BASE, (; kwargs...)))
end

@inline _nelson_activation(x, sat) = x / (sat + x)
@inline _nelson_repression(x, sat) = sat / (sat + x)

"""Relative estrogen concentration used by the public MATLAB model."""
function nelson_estrogen(t, p::NelsonParams; surgical::Bool=false)
    if surgical
        return t >= p.t_m ? (1 - p.k_syn / p.kappa_E) *
            exp(-p.kappa_E * (t - p.t_m)) + p.k_syn / p.kappa_E : 1.0
    end
    return t < p.t_m ? 1.0 : 1.0 / (1.0 + (t - p.t_m) / p.tau_E)
end

"""Return the six cell-state residuals at the source-model steady state."""
function _nelson_steady_rhs!(du, x, p::NelsonParams, new_effects::Bool)
    PB, PC, C, B, S, Sc = x
    # get_initial_condition.m fixes E=1 and t=0, even when the later
    # trajectory is a surgical-menopause simulation.
    E = 1.0
    r = C
    rep_sc_pb = _nelson_repression(Sc, p.sc_PB)
    rep_e_pc = _nelson_repression(E, p.e_PC)
    act_sc_pc = _nelson_activation(Sc, p.sc_PC)
    dPB = 1.0 - PB * rep_sc_pb * p.omega_PB
    dPC = 1.0 - PC * rep_e_pc * act_sc_pc * p.omega_PC
    dB = PB * rep_sc_pb * p.omega_PB - (p.eta_B + p.omega_B) * B
    dSc = _nelson_repression(E, p.e_Sc) * S - p.kappa_Sc * Sc
    dC = PC * p.omega_PC * rep_e_pc * act_sc_pc - C * p.eta_C
    dS = p.omega_B * B - p.eta_S * S
    if new_effects
        active = 0.0 > p.t_m
        eta_surg = active ? p.eta_S + p.eta_S * p.eta_ovx *
            exp(-p.tau * abs(0.0 - p.t_m)) : p.eta_S
        omega_surg = active ? p.omega_PC + p.omega_ovx *
            exp(-p.tau * abs(0.0 - p.t_m)) : p.omega_PC
        dPC = 1.0 - PC * rep_e_pc * act_sc_pc * omega_surg
        dC = PC * omega_surg * rep_e_pc * act_sc_pc - C * p.eta_C
        dS = p.omega_B * B - eta_surg * S
    end
    du[1] = dPB
    du[2] = dPC
    du[3] = dC
    du[4] = dB
    du[5] = dS
    du[6] = dSc
    return du
end

"""Steady-state initial condition used by the MATLAB source."""
function nelson_initial_state(p::NelsonParams=NelsonParams(); new_effects::Bool=false)
    # The source uses fsolve from an all-ones guess.  DynamicSS is used here
    # rather than hard-coding a state so parameter overrides remain useful.
    x0 = ones(6)
    f!(du, x, pp) = _nelson_steady_rhs!(du, x, pp[1], pp[2])
    prob = NonlinearProblem(f!, x0, (p, new_effects))
    sol = solve(prob, DifferentialEquations.NewtonRaphson(); abstol=1e-12, reltol=1e-12)
    SciMLBase.successful_retcode(sol) || error("Nelson steady-state solve failed: $(sol.retcode)")
    return vcat(Float64.(sol.u), 1.0)
end

@inline function nelson_rhs!(du, u, p::NelsonParams, t;
                             surgical::Bool=false, new_effects::Bool=false)
    PB, PC, C, B, S, Sc, Bd = u
    E = nelson_estrogen(t, p; surgical=surgical)
    r = C
    rep_sc_pb = _nelson_repression(Sc, p.sc_PB)
    rep_e_pc = _nelson_repression(E, p.e_PC)
    act_sc_pc = _nelson_activation(Sc, p.sc_PC)
    dPB = 1.0 - PB * rep_sc_pb * p.omega_PB
    dPC = 1.0 - PC * rep_e_pc * act_sc_pc * p.omega_PC
    dB = PB * rep_sc_pb * p.omega_PB - (p.eta_B + p.omega_B) * B
    dSc = _nelson_repression(E, p.e_Sc) * S - p.kappa_Sc * Sc
    dBd = B * p.lambda_B * _nelson_repression(Sc, p.sc_Omega) *
          (1.0 + p.nu_Omega * _nelson_activation(r, p.r_Omega)) - C * p.lambda_C
    dC = PC * p.omega_PC * rep_e_pc * act_sc_pc - C * p.eta_C
    dS = p.omega_B * B - p.eta_S * S
    if new_effects
        active = t > p.t_m
        eta_surg = active ? p.eta_S + p.eta_S * p.eta_ovx *
            exp(-p.tau * abs(t - p.t_m)) : p.eta_S
        omega_surg = active ? p.omega_PC + p.omega_ovx *
            exp(-p.tau * abs(t - p.t_m)) : p.omega_PC
        dPC = 1.0 - PC * rep_e_pc * act_sc_pc * omega_surg
        dC = PC * omega_surg * rep_e_pc * act_sc_pc - C * p.eta_C
        dS = p.omega_B * B - eta_surg * S
    end
    du[1] = dPB
    du[2] = dPC
    du[3] = dC
    du[4] = dB
    du[5] = dS
    du[6] = dSc
    du[7] = dBd
    return du
end

function nelson_problem(p::NelsonParams=NelsonParams();
                        u₀=nelson_initial_state(p), tspan=(0.0, 365.0),
                        surgical::Bool=false, new_effects::Bool=false,
                        kwargs...)
    f! = (du, u, pp, t) -> nelson_rhs!(du, u, pp, t;
                                        surgical=surgical, new_effects=new_effects)
    return ODEProblem(f!, u₀, tspan, p; kwargs...)
end

function simulate_nelson(p::NelsonParams=NelsonParams();
                         u₀=nothing, tspan=(0.0, 365.0),
                         surgical::Bool=false, new_effects::Bool=false,
                         kwargs...)
    initial = isnothing(u₀) ? nelson_initial_state(p; new_effects) : collect(u₀)
    prob = nelson_problem(p; u₀=initial, tspan, surgical, new_effects)
    sol = solve(prob, Rodas5P(autodiff=false);
                abstol=1e-10, reltol=1e-10, saveat=1.0, kwargs...)
    SciMLBase.successful_retcode(sol) || error("Nelson solve failed: $(sol.retcode)")
    return sol
end
