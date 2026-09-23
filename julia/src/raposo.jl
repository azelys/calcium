"""Raposo, Sobrinho & Ferreira (2002) minimal calcium-homeostasis model.

This is a direct Julia translation of the public Physiome CellML artifact. It
has 11 states and represents serum calcium/phosphate, intracellular
phosphate, PTH, exchangeable bone calcium/phosphate, renal 1-alpha-hydroxylase,
calcitriol, intestinal calcium transport, parathyroid secretory mass, and
phosphate kidney transport.

Source CellML artifact:
https://models.cellml.org/workspace/raposo_sobrinho_ferreira_2002/@@rawfile/3dcc0383f02438af3b98d3e74b025c870ff95bf1/raposo_2002.cellml
"""
Base.@kwdef struct RaposoParams{T<:Real}
    J_bp_Ca::T = 3.3
    J_pu_Ca::T = 0.21
    J_pu_P::T = 1.27
    C_PTH_p::T = 3.85
    k_PTH::T = 100.0
    S_E::T = 1.0
    k_E::T = 0.05
    k_D::T = 0.1
    k2_P::T = 1.0
    J_P_ing::T = 1.0
    Stoic_Ca_P::T = 0.464
    k_Ca_i::T = 1.0
    C_Ca_i::T = 1.0
    C_Ca_p::T = 2.4
    C_P_p::T = 1.2
    k_Ca_b::T = 1.0
    k3_P::T = 1.0
    k4_P::T = 51.8
    C_D_p::T = 90.0
    Y_Ca_i_2plus_Max::T = 1.0
    Y_PTH_p_2plus_Max::T = 1.0
    Y_PTH_p_2minus_Max::T = 1.0
    Y_Ca_p_1plus_minus_Max::T = 0.02
    Y_i_D_p_1plus_minus_Max::T = 0.02
    Y_PT_D_p_1plus_minus_Max::T = 0.01
    Y_k_Jext_1plus_minus_Max::T = 0.01
    X_R_PTH_Ca::T = 1.0
    X_R_i_Ca::T = 1.0
    X_R_PTH_D::T = 90.0
    X_R_E_PTH::T = 3.85
    X_R_i_D::T = 90.0
    X_R_P_k_PTH::T = 90.0
    b_PTH_Ca::T = 0.05
    b_PTH_D::T = 0.03
    b_E_PTH::T = 0.55
    b_i_D::T = 0.03
    a::T = 0.85
    d::T = 0.15
end

"""Initial state in the state order used by the CellML artifact."""
function raposo_initial_state(::RaposoParams{T}=RaposoParams{Float64}()) where {T}
    return T[1.0, 1.0, 3226.0, 1.0, 100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
end

@inline function raposo_aux(u, p::RaposoParams)
    Q_Ca_p, Q_P_p, Q_P_c, Q_PTH_p, Q_Ca_b, Q_P_b,
    Q_E_k, Q_D_p, Q_TCa_i, Q_C_PT, Q_TP_k = u

    J_bp_P = p.Stoic_Ca_P * p.J_bp_Ca
    Y_Ca_i_2plus = p.Y_Ca_i_2plus_Max *
                   (p.C_Ca_i / (p.C_Ca_i + p.X_R_i_Ca))
    Y_Ca_p_1minus = p.Y_Ca_p_1plus_minus_Max *
                    (p.a * (1 - tanh(p.b_PTH_Ca *
                                     (p.C_Ca_p - p.X_R_PTH_Ca))) + p.d)
    Y_i_D_p_1plus = p.Y_i_D_p_1plus_minus_Max *
                    (p.a * (1 + tanh(p.b_i_D *
                                     (p.C_D_p - p.X_R_i_D))) + p.d)
    Y_PT_D_p_1plus = p.Y_PT_D_p_1plus_minus_Max *
                     (p.a * (1 + tanh(p.b_PTH_D *
                                      (p.C_D_p - p.X_R_PTH_D))) + p.d)
    Y_k_Jext_1plus = p.Y_k_Jext_1plus_minus_Max *
                     (p.a * (1 + tanh(p.b_E_PTH *
                                      (p.C_P_p - p.X_R_E_PTH))) + p.d)
    Y_PTH_p_2plus = p.Y_PTH_p_2plus_Max *
                    (p.C_PTH_p / (p.C_PTH_p + p.X_R_P_k_PTH))
    Y_PTH_p_2minus = p.Y_PTH_p_2minus_Max *
                     (p.X_R_P_k_PTH / (p.C_PTH_p + p.X_R_P_k_PTH))
    Y_i_D_p_1minus = p.Y_i_D_p_1plus_minus_Max *
                     (p.a * (1 - tanh(p.b_i_D *
                                      (p.C_D_p - p.X_R_i_D))) + p.d)
    Y_PT_D_p_1minus = p.Y_PT_D_p_1plus_minus_Max *
                      (p.a * (1 - tanh(p.b_PTH_D *
                                       (p.C_D_p - p.X_R_PTH_D))) + p.d)

    J_pb_Ca = p.k_Ca_b * Q_Ca_b
    J_pb_P = p.Stoic_Ca_P * J_pb_Ca
    J_ip_Ca = (Y_Ca_i_2plus * (1 - Q_TCa_i) * Y_i_D_p_1plus +
               p.k_Ca_i * (p.C_Ca_i - p.C_Ca_p)) - Q_TCa_i * Y_i_D_p_1minus
    S_PTH = Y_Ca_p_1minus *
            (Y_PT_D_p_1plus * (1 - Q_C_PT) - Y_PT_D_p_1minus * Q_C_PT)
    J_pc_P = p.k4_P * p.C_P_p
    J_cp_P = p.k3_P * Q_P_c
    P_thr = Y_PTH_p_2minus * Q_TP_k
    P_T = P_thr / p.C_P_p
    Ca_thr = 1.95 + Y_PTH_p_2plus
    Ca_T = Ca_thr / p.C_Ca_p

    return (; Q_Ca_p, Q_P_p, Q_P_c, Q_PTH_p, Q_Ca_b, Q_P_b, Q_E_k,
            Q_D_p, Q_TCa_i, Q_C_PT, Q_TP_k,
            J_bp_Ca=p.J_bp_Ca, J_pb_Ca, J_pu_Ca=p.J_pu_Ca,
            J_ip_Ca, J_bp_P, J_pb_P, J_pu_P=p.J_pu_P,
            J_ip_P=p.J_P_ing, J_pc_P, J_cp_P, S_PTH,
            Y_Ca_i_2plus, Y_Ca_p_1minus, Y_i_D_p_1plus,
            Y_PT_D_p_1plus, Y_k_Jext_1plus, Y_PTH_p_2plus,
            Y_PTH_p_2minus, Y_i_D_p_1minus, Y_PT_D_p_1minus,
            P_thr, P_T, Ca_thr, Ca_T)
end

function raposo_rhs!(du, u, p::RaposoParams, t)
    a = raposo_aux(u, p)
    du[1] = (a.J_bp_Ca + a.J_ip_Ca) - (a.J_pb_Ca + a.J_pu_Ca)
    du[2] = (a.J_bp_P + a.J_ip_P + a.J_cp_P) -
            (a.J_pb_P + a.J_pu_P + a.J_pc_P)
    du[3] = a.J_pc_P - a.J_cp_P
    du[4] = a.S_PTH - p.k_PTH * p.C_PTH_p
    du[5] = a.J_pb_Ca - a.J_bp_Ca
    du[6] = a.J_pb_P - a.J_bp_P
    du[7] = p.S_E - p.k_E * u[7]
    du[8] = u[7] - p.k_D * u[8]
    du[9] = (1 - u[9]) * a.Y_i_D_p_1plus - u[9] * a.Y_i_D_p_1minus
    du[10] = (1 - u[10]) * a.Y_PT_D_p_1plus - u[10] * a.Y_PT_D_p_1minus
    du[11] = (1 - u[11]) * a.Y_k_Jext_1plus - u[11] * p.k2_P
    return nothing
end

function raposo_problem(p::RaposoParams=RaposoParams();
                        u₀=raposo_initial_state(p), tspan=(0.0, 168.0))
    length(u₀) == 11 || throw(ArgumentError("u₀ must have eleven states"))
    return ODEProblem(raposo_rhs!, collect(u₀), tspan, p)
end

"""Solve the Raposo CellML translation with a reproducible default solver."""
function simulate_raposo(p::RaposoParams=RaposoParams(); kwargs...)
    sol = solve(raposo_problem(p; kwargs...), Tsit5();
                abstol=1e-11, reltol=1e-10, saveat=0.1)
    SciMLBase.successful_retcode(sol) || error("Raposo solve failed: $(sol.retcode)")
    return sol
end
