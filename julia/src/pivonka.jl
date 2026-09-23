"""Pivonka et al. (2008) RANKL/OPG/TGF-beta bone-cell model.

This is the four-state CellML translation used to reproduce the published
model-structure simulations: responding osteoblasts (`OB_p`), active
osteoblasts (`OB_a`), active osteoclasts (`OC_a`), and bone volume (`BV`).
PTH, OPG, and RANKL are algebraic quantities in this particular artifact.

Source artifact:
https://models.physiomeproject.org/workspace/pivonka_zimak_smith_gardiner_dunstan_sims_martin_mundy_2008/rawfile/9848d85152c7ff7bae01d8ed69a9464992d8acf5/pivonka_2008_matlab.cellml
"""
Base.@kwdef struct PivonkaParams{T<:Real}
    D_OB_u::T = 7.0e-4
    D_OB_p::T = 2.674e-1
    A_OB_a::T = 1.890e-1
    A_OC_a::T = 7.000e-1
    BV₀::T = 100.0
    OC_a_initial::T = 1.0
    OB_a_initial::T = 1.0
    k_form::T = 1.571
    k_res::T = 200.0
    KD_TGF_repress::T = 1.416e-3
    KD_TGF_activate::T = 4.545e-3
    α_TGF::T = 1.0
    β_PTH::T = 2.500e2
    Deg_PTH::T = 8.600e1
    P_PTH_d::T = 0.0
    KD_PTH_repress::T = 1.5e2
    KD_PTH_activate::T = 2.226e-1
    β_OPG::T = 1.625e8
    β1_OB_p_OPG::T = 0.0
    β2_OB_a_OPG::T = 1.0
    OPG_max::T = 2.000e8
    P_OPG_d::T = 0.0
    Deg_OPG::T = 3.500e-1
    K_A2_RANKL::T = 3.412e-2
    K_A1_RANKL::T = 1.000e-3
    RANK::T = 1.000e1
    R_RANKL::T = 2.703e6
    P_RANKL_d::T = 0.0
    R1_OB_p_RANKL::T = 1.0
    R2_OB_a_RANKL::T = 0.0
    β_RANKL::T = 1.684e4
    Deg_RANKL::T = 1.013e1
    KD_RANKL_activate::T = 4.457
    D_OC_p_early::T = 2.100e-2
    D_OC_p_late::T = 2.100e-3
    early_end::T = 100.0
end

function pivonka_aux(u, p::PivonkaParams, t)
    OB_p, OB_a, OC_a, _ = u

    pi_TGF_beta_act = (p.α_TGF * OC_a) /
                      (p.KD_TGF_activate + p.α_TGF * OC_a)
    pi_TGF_beta_rep = 1 /
                      (1 + (p.α_TGF * OC_a) / p.KD_TGF_repress)

    PTH_tot = (p.β_PTH + p.P_PTH_d) / p.Deg_PTH
    pi_PTH_act = PTH_tot / (PTH_tot + p.KD_PTH_activate)
    pi_PTH_rep = 1 / (1 + PTH_tot / p.KD_PTH_repress)

    OPG_eff = (p.β1_OB_p_OPG * p.β_OPG * OB_p +
               p.β2_OB_a_OPG * p.β_OPG * OB_a) * pi_PTH_rep
    OPG = (p.P_OPG_d + OPG_eff) /
          (OPG_eff / p.OPG_max + p.Deg_OPG)

    K_RANKL_OB_p = p.R1_OB_p_RANKL * p.R_RANKL
    K_RANKL_OB_a = p.R2_OB_a_RANKL * p.R_RANKL
    RANKL_eff = (K_RANKL_OB_p * OB_p + K_RANKL_OB_a * OB_a) * pi_PTH_act
    RANKL_tot = 1 + p.K_A2_RANKL * p.RANK + p.K_A1_RANKL * OPG
    RANKL = (p.β_RANKL + p.P_RANKL_d) /
             (RANKL_tot * (p.β_RANKL / RANKL_eff + p.Deg_RANKL))
    RANKL_RANK = p.K_A2_RANKL * RANKL * p.RANK
    pi_RANKL_act = RANKL_RANK /
                   (p.KD_RANKL_activate + RANKL_RANK)

    D_OC_p = (0 < t <= p.early_end) ? p.D_OC_p_early : p.D_OC_p_late
    return (; pi_TGF_beta_act, pi_TGF_beta_rep, PTH_tot, pi_PTH_act,
            pi_PTH_rep, OPG_eff, OPG, K_RANKL_OB_p, K_RANKL_OB_a,
            RANKL_eff, RANKL_tot, RANKL, RANKL_RANK, pi_RANKL_act,
            D_OC_p)
end

function pivonka_rhs!(du, u, p::PivonkaParams, t)
    OB_p, OB_a, OC_a, _ = u
    any(x -> x <= 0, (OB_p, OB_a, OC_a)) &&
        throw(DomainError(u, "Pivonka cell populations must remain positive"))

    a = pivonka_aux(u, p, t)
    du[1] = p.D_OB_u * a.pi_TGF_beta_act -
            p.D_OB_p * a.pi_TGF_beta_rep * OB_p
    du[2] = p.D_OB_p * a.pi_TGF_beta_rep * OB_p - p.A_OB_a * OB_a
    du[3] = a.D_OC_p * a.pi_RANKL_act -
            p.A_OC_a * a.pi_TGF_beta_act * OC_a
    du[4] = p.k_form * (OB_a - p.OB_a_initial) -
            p.k_res * (OC_a - p.OC_a_initial)
    return nothing
end

function pivonka_problem(p::PivonkaParams=PivonkaParams();
                         u₀=[6.194e-4, 5.584e-4, 8.070e-4, p.BV₀],
                         tspan=(0.0, 140.0))
    length(u₀) == 4 || throw(ArgumentError("u₀ must have four states"))
    return ODEProblem(pivonka_rhs!, collect(u₀), tspan, p)
end

"""Solve the Pivonka CellML translation with a reproducible default solver."""
function simulate_pivonka(p::PivonkaParams=PivonkaParams(); kwargs...)
    sol = solve(pivonka_problem(p; kwargs...), Tsit5();
                abstol=1e-12, reltol=1e-10, saveat=0.1)
    SciMLBase.successful_retcode(sol) || error("Pivonka solve failed: $(sol.retcode)")
    return sol
end
