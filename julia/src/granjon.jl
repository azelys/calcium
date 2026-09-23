"""
Granjon et al. (2017) calcium/phosphate/FGF23 model.

This is a direct Julia translation of the public CaPO4Sim compiled core
(`compiled_core.c`) and its accompanying R core.  The model is a rat
whole-body mineral-homeostasis model, with 22 dynamic states and explicit
chemical-complex pools.  Time is measured in minutes.

The public CaPO4Sim package is the source artifact for this port.  The
translation intentionally retains its conventions, including quantities for
plasma PTH (converted to concentration inside the RHS), the disease-state
initial conditions, and the intervention switches.
"""

const GRANJON_STATE_NAMES = (
    :PTH_g, :PTH_p, :D3_p, :FGF_p, :Ca_p, :Ca_f, :Ca_b,
    :PO4_p, :PO4_f, :PO4_b, :PO4_c, :CaHPO4_p, :CaH2PO4_p,
    :CPP_p, :CaHPO4_f, :CaH2PO4_f, :CaProt_p, :NaPO4_p,
    :Ca_tot, :PO4_tot, :EGTA_p, :CaEGTA_p,
)

"""Parameters for the Granjon/CaPO4Sim model.

The values are kept in a NamedTuple so that the 133 parameters from the
public C artifact remain visible and can be overridden without changing a
large positional parameter vector.  For example,
`GranjonParams(k_prod_PTHg=0.0)` reproduces the hypoparathyroidism switch.
"""
struct GranjonParams{T<:NamedTuple}
    values::T
end

Base.getproperty(p::GranjonParams, s::Symbol) =
    s === :values ? getfield(p, :values) : getproperty(getfield(p, :values), s)

const _GRANJON_FIXED = (
    # Parameters omitted from cap_fixed_parameters.R because they are disease
    # inputs in the app, not immutable constants in the core.
    k_prod_PTHg = 4.192,
    D3_inact = 2.5e-5,
    Vp = 1e-2,
    Vc = 1e-7,
    k_deg_PTHg = 3.5e-2,
    k_deg_PTHp = 1.3,
    K_Ca = 1.16,
    beta_exo_PTHg = 5.9e-2,
    gamma_exo_PTHg = 5.8e-2,
    gamma_prod_D3 = 8.33e5,
    n1_exo = 1e2,
    n2_exo = 3e1,
    rho_exo = 1e1,
    R = 1.1,
    K_prod_PTH_P = 2.4,
    n_prod_Pho = 3.0,
    k_conv_min = 8.8e-6,
    delta_conv_max = 14.04e-5,
    gamma_ca_conv = 3e-1,
    gamma_D3_conv = 3e6,
    gamma_P_conv = 2e-1,
    gamma_FGF_conv = 2e7,
    K_conv = 1.575e-8,
    n_conv = 6.0,
    k_deg_D3 = 1e-3,
    gamma_deg_PTH = 1e9,
    gamma_deg_FGF = 1.265e8,
    k_prod_FGF = 6.902e-11,
    delta_max_prod_D3 = 10.0,
    K_prod_D3 = 5.64e-7,
    n_prod_FGF = 5.0,
    K_prod_P = 1.6,
    k_deg_FGF = 1.4e-2,
    K_abs_D3 = 6.4e-7,
    I_Ca = 2.2e-3,
    n_abs = 2.0,
    I_P = 1.55e-3,
    Lambda_ac_Ca = 5.5e-4,
    Lambda_ac_P = 2.75e-4,
    Lambda_res_min = 1e-4,
    delta_res_max = 6e-4,
    K_res_D3 = 6e-7,
    K_res_PTH = 6.12e-9,
    n_res = 2.0,
    k_p_Ca = 4.4e-1,
    k_f_Ca = 2.34e-3,
    k_p_P = 13.5,
    k_f_P = 2.5165e-1,
    k_pc = 1.875e-1,
    k_cp = 1e-3,
    k_f_CaHPO4 = 1.1373e3,
    k_d_CaHPO4 = 1.6667e3,
    k_f_CaH2PO4 = 5.294e1,
    k_d_CaH2PO4 = 1.6667e3,
    f1 = 7.6e-1,
    f2 = 3.3e-1,
    f3 = 8e-2,
    k_fet = 3e-1,
    k_c_CPP = 3.0,
    k_f_CaProt = 1.901976e2,
    k_d_CaProt = 1.6667e3,
    N_Prot = 20.0,
    Prot_tot_p = 6e-1,
    k_f_NaHPO4 = 7.8432,
    k_f_NaH2PO4 = 4.9020,
    Na = 142.0,
    k_d_NaHPO4 = 1.6667e3,
    k_d_NaH2PO4 = 1.6667e3,
    c = 6.2e-1,
    d = 3.8e-1,
    GFR = 2e-3,
    lambda_reabs_PT_0 = 6.4e-1,
    delta_PT_max = 1e-2,
    PTH_ref = 1.5e-8,
    n_PT = 5.0,
    lambda_TAL_0 = 2.25e-1,
    delta_CaSR_max = 1.75e-2,
    Ca_ref = 1.33,
    n_TAL = 4.0,
    K_TAL_PTH = 4.2e-9,
    delta_PTH_max = 7.5e-3,
    lambda_DCT_0 = 9e-2,
    delta_DCT_max = 1e-2,
    K_DCT_PTH = 6.3e-9,
    K_DCT_D3 = 4.8e-7,
    n_reabs_P = 5.0,
    lambda_PT_0 = 5.5e-1,
    delta_PTH_max_P = 1e-1,
    K_PT_PTH = 2e-8,
    delta_P_max = 5e-2,
    K_PT_P = 1.6,
    delta_FGF_max = 2e-1,
    K_PT_FGF = 2e-8,
    lambda_DCT_P = 1e-1,
    pH = 7.4,
    pKa = 6.8,
    D3_norm = 1e-9,
    PTH_g_norm = 1e-9,
    PTH_p_norm = 1e-9,
    FGF_p_norm = 1e-9,
    Ca_p_norm = 1.0,
    Ca_f_norm = 1.0,
    Ca_b_norm = 1.0,
    Pho_p_norm = 1.0,
    Pho_c_norm = 1.0,
    Pho_f_norm = 1.0,
    Pho_b_norm = 1.0,
    HPO4_norm = 1.0,
    H2PO4_norm = 1.0,
    CaHPO4_norm = 1.0,
    CaH2PO4_norm = 1.0,
    CPP_norm = 1.0,
    EGTA_norm = 1.0,
    Ca_EGTA_norm = 1.0,
    k_on_egta = 0.0,
    k_off_egta = 0.0,
    k_inject_egta = 0.0,
    K_sp_DCPD = 1.87e-7,
    # Calculated parameters in the public compiled core.
    r = 10.0^(7.4 - 6.8),
    a = 10.0^(7.4 - 6.8) / (1.0 + 10.0^(7.4 - 6.8)),
    b = (10.0^(7.4 - 6.8) / (1.0 + 10.0^(7.4 - 6.8))) / 10.0^(7.4 - 6.8),
    # Simulated parameters / interventions.
    PTX_coeff = 1.0,
    t_start = 0.0,
    t_stop = 0.0,
    Ca_inject = 0.0,
    Ca_food = 0.0,
    D3_inject = 0.0,
    P_inject = 0.0,
    P_food = 0.0,
    D3_intake_reduction = 0.0,
    Bispho = 1.0,
    Furo = 1.0,
    Cinacal = 0.0,
)

"""Construct model parameters, optionally overriding any source parameter."""
function GranjonParams(; kwargs...)
    unknown = setdiff(keys(kwargs), keys(_GRANJON_FIXED))
    isempty(unknown) || throw(ArgumentError("unknown Granjon parameter(s): $(collect(unknown))"))
    return GranjonParams(merge(_GRANJON_FIXED, (; kwargs...)))
end

"""Initial states distributed with the public CaPO4Sim model-engine artifact."""
function granjon_initial_state(case::Symbol = :hypopara)
    states = Dict(
        :hyperD3 => [515.915169364118, 0.00397042635508138, 1461.92617285272,
                     25.3337407792577, 1.71020100044328, 2.60376622892459,
                     295.278505778431, 1.14633769338327, 0.614292303690962,
                     90.8807303779351, 2.14938317804532, 0.116413743350791,
                     0.004124765017, 0.0119548564764044, 0.094994782741542,
                     0.0033658496446347, 1.95951821864912, 0.708344738497554,
                     3.51511258393703, 2.31007579672542, 0.0, 0.0],
        :hypoD3 => [536.407163768036, 0.243446328145703, 0.708899499542378,
                    4.9300000000001, 0.695793617379871, 1.05933933849388,
                    289.342588368956, 1.23669041585042, 0.662709983082683,
                    101.690376661555, 2.318817031696, 0.0511444769808592,
                    0.00180357293015262, 0.00525151925412956, 0.0417343912771445,
                    0.00147173307475654, 0.882726807166423, 0.764410997131242,
                    1.34961999371141, 2.38220098214677, 0.0, 0.0],
        :hypopara => [3.95252516672997e-323, 0.0, 135.276975798419,
                      4.95093399450975, 0.58930437901727, 0.897227282346784,
                      262.109976789696, 1.84336513502324, 0.987810354796314,
                      96.7227951807534, 3.45561176579289, 0.0645666633064686,
                      0.0022768966074445, 0.00662971023082098, 0.0526879760189357,
                      0.00185800330553375, 0.756140786128736, 1.13940284709278,
                      1.13181843529073, 3.37914125226075, 0.0, 0.0],
        :php1 => [52212.6087730569, 0.401636317436607, 3430.96913479457,
                  23.2427383681041, 2.05154828974457, 3.12343953331389,
                  279.309958624461, 0.945744142062095, 0.506799423833081,
                  88.2446482030736, 1.77322669634131, 0.115322144533028,
                  0.00406675188413527, 0.0118412874367777, 0.094103263668052,
                  0.00331848342202493, 2.27643333718403, 0.584574128970582,
                  4.17211181078241, 1.98444845488664, 0.0, 0.0],
    )
    haskey(states, case) || throw(ArgumentError("unknown Granjon initial-state case: $case"))
    return copy(states[case])
end

@inline _pow(x, n) = x^n

"""Return Granjon fluxes and the state derivative for one time point."""
function _granjon_core(t, u, p::GranjonParams)
    # The public C core treats plasma PTH as a quantity in the state vector
    # and converts it to concentration locally.
    PTHp = u[2] / p.Vp

    k_inject_P = 0.0
    k_inject_Ca = 0.0
    k_inject_D3 = 0.0
    k_inject_FGF = 0.0
    k_inject_PTH = 0.0
    I_Ca = p.I_Ca
    I_P = p.I_P
    if p.t_stop != 0.0 && t >= p.t_start && t <= p.t_stop
        if p.Ca_inject != 0.0
            k_inject_Ca = p.Ca_inject
        elseif p.Ca_food != 0.0
            I_Ca = p.Ca_food
        elseif p.D3_inject != 0.0
            k_inject_D3 = p.D3_inject * p.Vp
        elseif p.P_inject != 0.0
            k_inject_P = p.P_inject
        elseif p.P_food != 0.0
            I_P = p.P_food
        end
    end

    PTHg_basal_synthesis = p.k_prod_PTHg * p.Vc / p.PTH_g_norm
    PTHg_synthesis_D3 = 1.0 / (1.0 + p.gamma_prod_D3 * p.D3_norm * u[3])
    PTHg_synthesis_PO4 = _pow(u[8], p.n_prod_Pho) /
        (_pow(p.K_prod_PTH_P / p.Pho_p_norm, p.n_prod_Pho) + _pow(u[8], p.n_prod_Pho))
    PTHg_synthesis = p.PTX_coeff * PTHg_basal_synthesis * PTHg_synthesis_D3 * PTHg_synthesis_PO4
    PTHg_degradation = p.k_deg_PTHg * u[1]
    n_Ca = p.n1_exo / (1.0 + exp(-p.rho_exo * p.Ca_p_norm * (p.R / p.Ca_p_norm - u[5]))) + p.n2_exo
    F_Ca = if p.Cinacal == 0.0
        p.beta_exo_PTHg - p.gamma_exo_PTHg * _pow(u[5], n_Ca) /
            (_pow(u[5], n_Ca) + _pow(p.K_Ca / p.Ca_p_norm, n_Ca))
    else
        p.beta_exo_PTHg - p.gamma_exo_PTHg
    end
    PTHg_exocytosis = F_Ca * u[1]
    PTHp_influx = PTHg_exocytosis * p.PTH_g_norm / p.PTH_p_norm
    PTHp_degradation = p.k_deg_PTHp * PTHp

    D3_basal_synthesis = (1.0 - p.D3_intake_reduction / 100.0) * p.k_conv_min * p.D3_inact / p.D3_norm
    D3_conv_PTH = (p.delta_conv_max * (p.D3_inact / p.D3_norm) * _pow(PTHp, p.n_conv)) /
        (_pow(PTHp, p.n_conv) + _pow(p.K_conv / p.PTH_p_norm, p.n_conv))
    D3_conv_Ca = 1.0 / (1.0 + p.gamma_ca_conv * p.Ca_p_norm * u[5])
    D3_conv_D3 = 1.0 / (1.0 + p.gamma_D3_conv * p.D3_norm * u[3])
    D3_conv_P = 1.0 / (1.0 + p.gamma_P_conv * p.Pho_p_norm * u[8])
    D3_conv_FGF = 1.0 / (1.0 + p.gamma_FGF_conv * p.FGF_p_norm * u[4])
    D3_synthesis = D3_basal_synthesis + D3_conv_PTH * D3_conv_Ca * D3_conv_D3 * D3_conv_P * D3_conv_FGF
    D3_degradation = (p.k_deg_D3 * (1.0 + p.gamma_deg_FGF * p.FGF_p_norm * u[4]) * u[3]) /
        (1.0 + p.gamma_deg_PTH * p.PTH_p_norm * PTHp)

    FGF_basal_synthesis = p.k_prod_FGF / p.FGF_p_norm
    FGF_D3_activ = p.delta_max_prod_D3 * _pow(u[3], p.n_prod_FGF) /
        (_pow(u[3], p.n_prod_FGF) + _pow(p.K_prod_D3 / p.D3_norm, p.n_prod_FGF))
    FGF_P_activ = u[8] / (u[8] + p.K_prod_P / p.Pho_p_norm)
    FGF_synthesis = FGF_basal_synthesis * (1.0 + FGF_D3_activ * FGF_P_activ)
    FGF_degradation = p.k_deg_FGF * u[4]

    Abs_intest_basal = 0.25 * I_Ca
    Abs_intest_D3 = (0.45 * I_Ca * _pow(u[3], p.n_abs)) /
        (_pow(u[3], p.n_abs) + _pow(p.K_abs_D3 / p.D3_norm, p.n_abs))
    Abs_intest = (Abs_intest_basal + Abs_intest_D3) / p.Ca_p_norm
    Abs_intest_basal_P = 0.4 * I_P
    Abs_intest_D3_P = (0.3 * I_P * _pow(u[3], p.n_abs)) /
        (_pow(u[3], p.n_abs) + _pow(p.K_abs_D3 / p.D3_norm, p.n_abs))
    Abs_intest_P = (Abs_intest_basal_P + Abs_intest_D3_P) / p.Pho_p_norm

    Rapid_storage_Ca = p.k_p_Ca * u[5] * p.Vp
    Rapid_release_Ca = p.k_f_Ca * u[6]
    Accretion_Ca = p.Lambda_ac_Ca * u[6]
    Rapid_storage_P = p.k_p_P * u[8] * p.Vp
    Rapid_release_P = p.k_f_P * u[9]
    Accretion_P = p.Lambda_ac_P * u[9]

    Resorption_basal = p.Lambda_res_min
    Resorption_PTH = (p.delta_res_max * 0.2 * _pow(PTHp, p.n_res)) /
        (_pow(PTHp, p.n_res) + _pow(p.K_res_PTH / p.PTH_p_norm, p.n_res))
    Resorption_D3 = (p.delta_res_max * 0.8 * _pow(u[3], p.n_res)) /
        (_pow(u[3], p.n_res) + _pow(p.K_res_D3 / p.D3_norm, p.n_res))
    Resorption = p.Bispho * (Resorption_basal + Resorption_PTH + Resorption_D3)
    Resorption_P = 0.3 * Resorption

    Reabs_PT_basal = p.lambda_reabs_PT_0
    Reabs_PT_PTH = p.delta_PT_max / (1.0 + _pow(PTHp * p.PTH_p_norm / p.PTH_ref, p.n_PT))
    Reabs_PT = Reabs_PT_basal + Reabs_PT_PTH
    Reabs_TAL_basal = p.lambda_TAL_0
    Reabs_TAL_CaSR = if p.Cinacal == 0.0
        p.delta_CaSR_max / (1.0 + _pow(u[5] * p.Ca_p_norm / p.Ca_ref, p.n_TAL))
    else
        p.delta_CaSR_max
    end
    Reabs_TAL_PTH = p.delta_PTH_max * PTHp / (PTHp + p.K_TAL_PTH / p.PTH_p_norm)
    Reabs_DCT_basal = p.lambda_DCT_0
    Reabs_DCT_PTH = (p.delta_DCT_max * 0.8 * PTHp) / (PTHp + p.K_DCT_PTH / p.PTH_p_norm)
    Reabs_DCT_D3 = (p.delta_DCT_max * 0.2 * u[3]) / (u[3] + p.K_DCT_D3 / p.D3_norm)
    Excretion = p.Furo * (1.0 - (Reabs_PT + Reabs_TAL_basal + Reabs_TAL_CaSR + Reabs_TAL_PTH +
                                  Reabs_DCT_basal + Reabs_DCT_PTH + Reabs_DCT_D3)) *
                p.GFR * (u[5] + u[12] + u[13])
    Reabs = (Reabs_PT + Reabs_TAL_basal + Reabs_TAL_CaSR + Reabs_TAL_PTH +
             Reabs_DCT_basal + Reabs_DCT_PTH + Reabs_DCT_D3) * p.GFR * (u[5] + u[12] + u[13])

    Reabs_PT_basal_P = p.lambda_PT_0
    Reabs_PT_PTH_P = (p.delta_PTH_max_P * _pow(p.K_PT_PTH / p.PTH_p_norm, p.n_reabs_P)) /
        (_pow(PTHp, p.n_reabs_P) + _pow(p.K_PT_PTH / p.PTH_p_norm, p.n_reabs_P))
    Reabs_PT_FGF_P = (p.delta_FGF_max * _pow(p.K_PT_FGF / p.FGF_p_norm, p.n_reabs_P)) /
        (_pow(u[4], p.n_reabs_P) + _pow(p.K_PT_FGF / p.FGF_p_norm, p.n_reabs_P))
    Reabs_PT_P = (p.delta_P_max * _pow(p.K_PT_P / p.Pho_p_norm, p.n_reabs_P)) /
        (_pow(u[8], p.n_reabs_P) + _pow(p.K_PT_P / p.Pho_p_norm, p.n_reabs_P))
    Reabs_DCT_basal_P = p.lambda_DCT_P
    Excretion_P = (1.0 - (Reabs_PT_basal_P + Reabs_PT_PTH_P + Reabs_PT_FGF_P + Reabs_PT_P + Reabs_DCT_basal_P)) *
                  p.GFR * (u[8] + u[12] + u[13] + u[18])
    Reabs_P = (Reabs_PT_basal_P + Reabs_PT_PTH_P + Reabs_PT_FGF_P + Reabs_PT_P + Reabs_DCT_basal_P) *
              p.GFR * (u[8] + u[12] + u[13] + u[18])

    Plasma_intra_Flux = p.k_pc * u[8] * p.Vp
    Intra_plasma_Flux = p.k_cp * u[11]

    k_form_CaHPO4 = p.k_f_CaHPO4 * u[5] * p.a * u[8] * p.f2^2
    k_diss_CaHPO4 = p.k_d_CaHPO4 * u[12]
    k_form_CaH2PO4 = p.k_f_CaH2PO4 * u[5] * p.b * u[8] * p.f2 * p.f1
    k_diss_CaH2PO4 = p.k_d_CaH2PO4 * u[13] * p.f1
    k_form_CaHPO4f = p.k_f_CaHPO4 * u[6] * p.a * u[9] * p.f2^2
    k_diss_CaHPO4f = p.k_d_CaHPO4 * u[15]
    k_form_CaH2PO4f = p.k_f_CaH2PO4 * u[6] * p.b * u[9] * p.f2 * p.f1
    k_diss_CaH2PO4f = p.k_d_CaH2PO4 * u[16] * p.f1
    k_fet_CaHPO4 = p.k_fet * u[12]
    k_fet_CaH2PO4 = p.k_fet * u[13] * p.f1
    CPP_degradation = p.k_c_CPP * u[14]
    k_form_CaProt = p.k_f_CaProt * u[5] * (p.N_Prot * p.Prot_tot_p - u[17])
    k_diss_CaProt = p.k_d_CaProt * u[17]
    k_form_NaPO4 = (p.a * p.k_f_NaHPO4 + p.b * p.k_f_NaH2PO4) * p.Na * u[8]
    k_diss_NaPO4 = (p.c * p.k_d_NaHPO4 + p.d * p.k_d_NaH2PO4) * u[18]
    EGTA_form = p.k_on_egta * u[5] * u[21]
    EGTA_diss = p.k_off_egta * u[22]

    du = zeros(eltype(u), 22)
    du[1] = PTHg_synthesis - PTHg_degradation - PTHg_exocytosis
    du[2] = k_inject_PTH + PTHp_influx - PTHp_degradation
    du[3] = k_inject_D3 + D3_synthesis - D3_degradation
    du[4] = k_inject_FGF + FGF_synthesis - FGF_degradation
    du[5] = (k_inject_Ca + Abs_intest + Resorption / p.Ca_p_norm - Rapid_storage_Ca + Rapid_release_Ca - Excretion) / p.Vp -
            k_form_CaProt + k_diss_CaProt - k_form_CaHPO4 * p.HPO4_norm +
            k_diss_CaHPO4 * p.CaHPO4_norm / p.Ca_p_norm - k_form_CaH2PO4 / p.Ca_p_norm +
            k_diss_CaH2PO4 * p.CaH2PO4_norm / p.Ca_p_norm - EGTA_form + EGTA_diss / p.Ca_p_norm
    du[6] = Rapid_storage_Ca - Rapid_release_Ca - Accretion_Ca -
            k_form_CaHPO4f * p.Ca_p_norm * p.HPO4_norm / p.CaHPO4_norm + k_diss_CaHPO4f -
            k_form_CaH2PO4f * p.Ca_p_norm * p.H2PO4_norm / p.CaH2PO4_norm + k_diss_CaH2PO4f
    du[7] = (Accretion_Ca - Resorption) / p.Ca_b_norm
    du[8] = (k_inject_P + Abs_intest_P + Resorption_P / p.Pho_p_norm - Rapid_storage_P + Rapid_release_P -
             Excretion_P - Plasma_intra_Flux + Intra_plasma_Flux / p.Pho_p_norm) / p.Vp -
            k_form_CaHPO4 * p.Ca_p_norm + k_diss_CaHPO4 * p.CaHPO4_norm / p.Pho_p_norm -
            k_form_CaH2PO4 * p.Ca_p_norm + k_diss_CaH2PO4 * p.CaH2PO4_norm / p.Pho_p_norm -
            k_form_NaPO4 + k_diss_NaPO4
    du[9] = Rapid_storage_P - Rapid_release_P - Accretion_P -
            k_form_CaHPO4f * p.Ca_p_norm * p.HPO4_norm / p.CaHPO4_norm + k_diss_CaHPO4f -
            k_form_CaH2PO4f * p.Ca_p_norm * p.H2PO4_norm / p.CaH2PO4_norm + k_diss_CaH2PO4f
    du[10] = (Accretion_P - Resorption_P) / p.Pho_b_norm
    du[11] = Plasma_intra_Flux / p.Pho_c_norm - Intra_plasma_Flux
    du[12] = k_form_CaHPO4 * p.Ca_p_norm * p.HPO4_norm / p.CaHPO4_norm - k_diss_CaHPO4 - k_fet_CaHPO4
    du[13] = k_form_CaH2PO4 * p.Ca_p_norm * p.H2PO4_norm / p.CaH2PO4_norm - k_diss_CaH2PO4 - k_fet_CaH2PO4
    du[14] = k_fet_CaHPO4 + k_fet_CaH2PO4 - CPP_degradation
    du[15] = k_form_CaHPO4f * p.Ca_p_norm * p.HPO4_norm / p.CaHPO4_norm - k_diss_CaHPO4f
    du[16] = k_form_CaH2PO4f * p.Ca_p_norm * p.H2PO4_norm / p.CaH2PO4_norm - k_diss_CaH2PO4f
    du[17] = k_form_CaProt - k_diss_CaProt
    du[18] = k_form_NaPO4 - k_diss_NaPO4
    du[19] = du[5] + du[12] + du[13] + du[17] + du[14]
    du[20] = du[8] + du[12] + du[13] + du[18] + du[14]
    du[21] = p.k_inject_egta / (p.Vp * p.EGTA_norm) - EGTA_form + EGTA_diss / p.EGTA_norm
    du[22] = EGTA_form / p.Ca_EGTA_norm - EGTA_diss

    aux = (
        U_Ca = Excretion, U_PO4 = Excretion_P,
        Abs_int_Ca = Abs_intest, Abs_int_PO4 = Abs_intest_P,
        Res_Ca = Resorption, Res_PO4 = Resorption_P,
        Ac_Ca = Accretion_Ca, Ac_PO4 = Accretion_P,
        Reabs_Ca = Reabs, Reabs_PO4 = Reabs_P,
        Ca_pf = Rapid_storage_Ca, Ca_fp = Rapid_release_Ca,
        PO4_pf = Rapid_storage_P, PO4_fp = Rapid_release_P,
        PO4_pc = Plasma_intra_Flux, PO4_cp = Intra_plasma_Flux,
        PTHg_synth = PTHg_synthesis, PTHg_synth_D3 = PTHg_synthesis_D3,
        PTHg_synth_PO4 = PTHg_synthesis_PO4, PTHg_exo_CaSR = F_Ca,
        PTHg_deg = PTHg_degradation, PTHg_exo = PTHg_exocytosis,
        PTHp_deg = PTHp_degradation, Reabs_PT_PTH = Reabs_PT_PTH,
        Reabs_TAL_CaSR = Reabs_TAL_CaSR, Reabs_TAL_PTH = Reabs_TAL_PTH,
        Reabs_DCT_PTH = Reabs_DCT_PTH, Reabs_DCT_D3 = Reabs_DCT_D3,
        Abs_int_D3 = Abs_intest_D3, Res_PTH = Resorption_PTH,
        Res_D3 = Resorption_D3, Reabs_PT_PO4_PTH = Reabs_PT_PTH_P,
        Reabs_PT_PO4_FGF = Reabs_PT_FGF_P,
    )
    return du, aux
end

function granjon_rhs!(du, u, p::GranjonParams, t)
    du .= first(_granjon_core(t, u, p))
    return nothing
end

function granjon_aux(u, p::GranjonParams, t = 0.0)
    return last(_granjon_core(t, u, p))
end

function granjon_problem(; case::Symbol = :hypopara, u0 = granjon_initial_state(case),
                         p::GranjonParams = GranjonParams(), tspan = (0.0, 1440.0), kwargs...)
    return ODEProblem(granjon_rhs!, u0, tspan, p; kwargs...)
end

function simulate_granjon(; case::Symbol = :hypopara, u0 = granjon_initial_state(case),
                          p::GranjonParams = GranjonParams(), tspan = (0.0, 1440.0),
                          saveat = nothing, kwargs...)
    prob = granjon_problem(; case, u0, p, tspan)
    alg = get(kwargs, :alg, Rodas5P(autodiff = false))
    solver_kwargs = (; (k => v for (k, v) in pairs(kwargs) if k != :alg)...)
    solve_kwargs = (; (k => v for (k, v) in pairs(solver_kwargs) if k != :dtmax)...,
                    dtmax = get(solver_kwargs, :dtmax, 1.0))
    return saveat === nothing ?
        solve(prob, alg; solve_kwargs...) :
        solve(prob, alg; saveat, solve_kwargs...)
end
