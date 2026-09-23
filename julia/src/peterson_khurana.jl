"""A 32-state Peterson/Riggs–Khurana-style integrated model.

This source-faithful Julia variant preserves the currently documented
equations, including the calcitriol balance term `+J14`, the hard-coded 99%
parathyroid-loss target, and the PTH/dose naming. A repaired model belongs in
a separate variant after numerical regression tests are in place.
"""

const _PETERSON_KHURANA_BASE = (
    k1=0.00000624,
    k2=0.112013,
    k3=0.00000624,
    k4=0.112013,
    V1=14.0,
    CaDay=88.0,
    FracJ14=0.107763,
    J14OCmax=0.543488,
    J14OCgam=1.6971,
    FracJ15=0.114376,
    kinRNKgam=0.151825,
    koutRNK=0.00323667,
    MOCratioGam=0.603754,
    Da=0.7 / 24,
    OBtgfGAM=0.0111319,
    koutTGF0=0.0000298449,
    koutTGFGam=0.919131,
    OCtgfGAM=0.593891,
    EmaxPicROB=3.9745,
    PicROBgam=1.80968,
    FracPicROB=0.883824,
    PicOBgam=0.122313,
    FracPicOB=0.000244818,
    EmaxPicOB=0.251636,
    E0Meff=0.388267,
    EmaxMeffOC=3.15667,
    kinOCgam=8.53065,
    EmaxPicOC=1.9746,
    FracPicOC=0.878215,
    PicOCgam=1.0168,
    E0RANKL=3.80338,
    EmaxL=0.469779,
    GFR=100 / 16.667,
    T16=1.06147,
    T64=0.05,
    T65=6.3,
    T67=1.54865,
    AlphOHgam=0.111241,
    k14a=0.0000244437,
    HApMRT=3.60609,
    koutL=0.00293273,
    OsteoEffectGam=0.173833,
    TESTPOWER=1.0,
    opgPTH50=3.85,
    IO=0.0,
    RX2Kout0=0.693,
    E0rx2Kout=0.125,
    EmaxPTHRX2x=5.0,
    E0crebKin=0.5,
    EmaxPTHcreb=3.39745,
    crebKout=0.00279513,
    bcl2Kout=0.693,
    ScaEffGam=0.9,
    PhosEff0=1.52493,
    PhosEff50=1.3021,
    PhosEffGam=8.25229,
    PO4inhPTHgam=0.0,
    T69=0.10,
    Reabs50=1.57322,
    T7=2.0,
    T9=90.0,
    T70=0.01,
    T71=0.03,
    T33=0.003,
    T34=0.037,
    T35=90.0,
    CaPOgam=1.0,
    T46=1.142,
    T52=0.365,
    OralPhos=10.5 / 24,
    F12=0.7,
    T49=51.8,
    T55=0.019268,
    PicOBgamkb=2.92375,
    MultPicOBkb=3.11842,
    FracPic0kb=0.764028,
    E0RUNX2kbEffFACT=1.01,
    RUNkbGAM=3.67798,
    RUNkbMaxFact=0.638114,
    RUNX20=10.0,
    Frackb=0.313186,
    T81=0.75,
    T87=0.0495,
    T0=1.58471,
    T28=0.9,
    T77=0.909359,
    T80=4.0,
    CtriolPTgam=12.5033,
    CtriolMax=4.1029,
    CtriolMin=0.9,
    PTout=0.0001604,
    FPTH=0.55,
    KAVITD=0.081,
    FVITD=0.9,
    OralCa=50 / 24,
    T57=75.0,
    T58=624.909,
    T59=11.7387,
    T61=96.25,
    IPTHint=0.0,
    IPTHinf=0.0,
    Pic0=0.228142,
    LsurvOCCgam=3.0923,
    EmaxLpth=1.30721,
    kO=15.8885,
    kb=0.000605516,
    LsurvOCgam=3.09023,
    FracOBfast=0.797629,
    TERIVC=94.0,
    TERICL=40.0,
    TERIKA=16.4,
)

"""Container for the 32-state variant parameters.

The fields are stored as a named tuple so source parameter names remain
visible (`p.k1`, `p.FracJ14`, and so on), while callers can override any
parameter with `PetersonKhuranaParams(k1=..., ...)`.
"""
struct PetersonKhuranaParams{N<:NamedTuple}
    values::N
end

function Base.getproperty(p::PetersonKhuranaParams, name::Symbol)
    name === :values && return getfield(p, :values)
    return getproperty(getfield(p, :values), name)
end

PetersonKhuranaParams(; kwargs...) =
    PetersonKhuranaParams(merge(_PETERSON_KHURANA_BASE, (; kwargs...)))

function peterson_khurana_derived_params(p::PetersonKhuranaParams, u₀)
    return PetersonKhuranaParams(merge(p.values, (
        Q0=u₀[25],
        OC0=u₀[16],
        RNK0=u₀[21],
        RANKL0=u₀[20],
        RNKL0=u₀[20],
        OB0=u₀[12] + u₀[13],
        ROB0=u₀[17],
        QboneInit=u₀[26],
        OPG0=u₀[24],
        RX20=u₀[27],
        CREB0=u₀[28],
        M0=u₀[22],
        TGFBact0=u₀[19],
        TGFB0=u₀[18],
    )))
end

function peterson_khurana_initial_state(p::PetersonKhuranaParams=PetersonKhuranaParams())
    OB = 0.00501324
    OBfast = OB * p.FracOBfast
    OBslow = OB * (1 - p.FracOBfast)
    M = p.k3 * 10.0 * 0.4 / p.k4
    N = p.k1 * 4.0 * 0.4 / p.k2
    return Float64[
        53.96, 0.5, 1.0, 1260.0, 0.0, 126.0, 32.90, 16.8,
        1.2375, 0.50, 1.0, OBfast, OBslow, 0.839, 3226.0,
        0.00115398, 0.00104122, p.Pic0 * 1000.0, p.Pic0, 0.4, 10.0,
        M, N, 4.0, 100.0, 24900.0, 10.0, 10.0, 100.0, 0.0, 0.0, 0.0,
    ]
end

@inline function peterson_khurana_rhs!(du, y, p::PetersonKhuranaParams, t)
    PTH, S, PTmax, B, SC, A, P, ECCPhos, T, R, HAp, OBfast, OBslow,
    PhosGut, IntraPO, OC, ROB1, TGFB, TGFBact, L, RNK, M, N, O, Q,
    Qbone, RX2, CREB, BCL2, TERISC, UCAL, VDORAL = y

    T13 = (p.CaDay / 24) / p.Q0
    T15 = p.CaDay / (2.35 * 14 * 24)
    T17 = 3.85 * p.T16 - 3.85
    Osteoclast = OC
    J14OC50 = exp(log((p.J14OCmax * p.OC0^p.J14OCgam / T13) -
                      p.OC0^p.J14OCgam) / p.J14OCgam)
    OCeqn = (p.J14OCmax * Osteoclast^p.J14OCgam) /
            (Osteoclast^p.J14OCgam + J14OC50^p.J14OCgam)
    kinRNK = (p.koutRNK * p.RNK0 + p.k3 * p.RNK0 * p.RANKL0 -
              p.k4 * p.M0) / p.TGFBact0^p.kinRNKgam
    MOCratio = M / Osteoclast
    MOCratio0 = p.M0 / p.OC0
    MOCratioEff = (MOCratio / MOCratio0)^p.MOCratioGam
    J14OCdepend = OCeqn * p.Q0 * p.FracJ14 * MOCratioEff
    J14 = T13 * p.Q0 * (1 - p.FracJ14) + J14OCdepend
    J41 = 0.464 * J14
    PicOCkin = p.Pic0
    bigDb = p.kb * p.OB0 * p.Pic0 / p.ROB0
    kinTGF = p.koutTGF0 * p.TGFB0
    koutTGFact = p.koutTGF0 * 1000
    koutTGFeqn = p.koutTGF0 * TGFB * (Osteoclast / p.OC0)^p.OCtgfGAM
    E0PicROB = p.FracPicROB * p.Pic0
    EC50PicROBparen = (p.EmaxPicROB * p.TGFBact0^p.PicROBgam /
                       (p.Pic0 - E0PicROB)) - p.TGFBact0^p.PicROBgam
    EC50PicROB = exp(log(EC50PicROBparen) / p.PicROBgam)
    PicROB = E0PicROB + p.EmaxPicROB * TGFBact^p.PicROBgam /
             (TGFBact^p.PicROBgam + EC50PicROB^p.PicROBgam)
    ROBin = (p.kb * p.OB0 / p.Pic0) * PicROB
    E0PicOB = p.FracPicOB * p.Pic0
    EC50PicOBparen = (p.EmaxPicOB * p.TGFBact0^p.PicOBgam /
                      (p.Pic0 - E0PicOB)) - p.TGFBact0^p.PicOBgam
    EC50PicOB = exp(log(EC50PicOBparen) / p.PicOBgam)
    PicOB = E0PicOB + p.EmaxPicOB * TGFBact^p.PicOBgam /
            (TGFBact^p.PicOBgam + EC50PicOB^p.PicOBgam)
    KPT = bigDb / PicOB
    D = ROB1
    EC50MeffOC = exp(log(p.M0^p.kinOCgam * p.EmaxMeffOC /
                         (1 - p.E0Meff) - p.M0^p.kinOCgam) / p.kinOCgam)
    MeffOC = p.E0Meff + p.EmaxMeffOC * M^p.kinOCgam /
             (M^p.kinOCgam + EC50MeffOC^p.kinOCgam)
    kinOC2 = p.Da * PicOCkin * MeffOC * p.OC0
    E0PicOC = p.FracPicOC * p.Pic0
    EC50PicOCparen = (p.EmaxPicOC * p.TGFBact0^p.PicOCgam /
                      (p.Pic0 - E0PicOC)) - p.TGFBact0^p.PicOCgam
    EC50PicOC = exp(log(EC50PicOCparen) / p.PicOCgam)
    PicOC = E0PicOC + p.EmaxPicOC * TGFBact^p.PicOCgam /
            (TGFBact^p.PicOCgam + EC50PicOC^p.PicOCgam)
    PiL0 = (p.k3 / p.k4) * p.RANKL0
    PiL = M / 10
    EC50survInPar = (p.E0RANKL - p.EmaxL) *
                    (PiL0^p.LsurvOCgam / (p.E0RANKL - 1)) -
                    PiL0^p.LsurvOCgam
    EC50surv = exp(log(EC50survInPar) / p.LsurvOCgam)
    LsurvOC = p.E0RANKL - (p.E0RANKL - p.EmaxL) *
              (PiL^p.LsurvOCgam / (PiL^p.LsurvOCgam + EC50surv^p.LsurvOCgam))
    KLSoc = p.Da * PicOC * LsurvOC
    C4 = PTH / p.V1
    T66 = (p.T67^p.AlphOHgam + 3.85^p.AlphOHgam) / 3.85^p.AlphOHgam
    k15a = p.k14a * p.QboneInit / p.Q0
    J14a = p.k14a * Qbone
    J15a = k15a * Q
    kLShap = 1 / p.HApMRT
    kHApIn = kLShap / p.OB0
    J15 = T15 * P * (1 - p.FracJ15) + T15 * P * p.FracJ15 * HAp
    J42 = 0.464 * J15
    OBfast0 = p.OB0 * p.FracOBfast
    Osteoblast = OBfast + OBslow
    kinLbase = p.koutL * p.RANKL0
    OsteoEffect = (Osteoblast / p.OB0)^p.OsteoEffectGam
    PTH50 = p.EmaxLpth * 3.85 - 3.85
    PTHconc = C4
    LpthEff = p.EmaxLpth * PTHconc /
             (PTH50 * OsteoEffect^p.TESTPOWER + PTHconc)
    kinL = kinLbase * OsteoEffect * LpthEff
    pObase = p.kO * p.OPG0
    pO = pObase * (D / p.ROB0) *
         ((PTHconc + p.opgPTH50 * (D / p.ROB0)) / (2 * PTHconc)) + p.IO
    RX2Kin = p.RX2Kout0 * p.RX20
    EC50PTHRX2x = (p.EmaxPTHRX2x * 3.85 /
                   (p.RX2Kout0 - p.E0rx2Kout)) - 3.85
    RX2Kout = p.E0rx2Kout + p.EmaxPTHRX2x * PTHconc /
              (PTHconc + EC50PTHRX2x)
    EC50PTHcreb = (p.EmaxPTHcreb * 3.85 /
                   (1 - p.E0crebKin)) - 3.85
    crebKin0 = p.crebKout * p.CREB0
    crebKin = crebKin0 * (p.E0crebKin + p.EmaxPTHcreb * PTHconc /
                          (PTHconc + EC50PTHcreb))
    bcl2Kin = RX2 * CREB * 0.693
    CaConc = P / 14
    C2 = ECCPhos / p.V1
    PO4inhPTH = (C2 / 1.2)^p.PO4inhPTHgam
    PhosEffTop = (p.PhosEff0 - 1) *
                 (1.2^p.PhosEffGam + p.PhosEff50^p.PhosEffGam)
    PhosEffBot = p.PhosEff0 * 1.2^p.PhosEffGam
    PhosEffMax = PhosEffTop / PhosEffBot
    PhosEff = p.PhosEff0 - PhosEffMax * p.PhosEff0 * C2^p.PhosEffGam /
              (C2^p.PhosEffGam + p.PhosEff50^p.PhosEffGam)
    PhosEffect = C2 > 1.2 ? PhosEff : 1.0
    T68 = T66 * C4^p.AlphOHgam /
          (p.T67^p.AlphOHgam * PO4inhPTH + C4^p.AlphOHgam)
    SE = p.T65 * T68 * PhosEffect
    C8 = B / p.V1
    C1 = P / p.V1
    T36 = p.T33 + (p.T34 - p.T33) *
          (C8^p.CaPOgam / (p.T35^p.CaPOgam + C8^p.CaPOgam))
    T37 = p.T34 - (p.T34 - p.T33) *
          (C8^p.CaPOgam / (p.T35^p.CaPOgam + C8^p.CaPOgam))
    CaFilt = 0.6 * 0.5 * p.GFR * C1
    ReabsMax = (0.3 * p.GFR * 2.35 - 0.149997) *
               (p.Reabs50 + 2.35) / 2.35
    ReabsPTHeff = (p.T16 * C4) / (C4 + T17)
    CaReabsActive = (ReabsMax * C1 / (p.Reabs50 + C1)) * ReabsPTHeff
    T20 = CaFilt - CaReabsActive
    T10 = p.T7 * C8 / (C8 + p.T9)
    J27a = (2 - T10) * T20
    J27 = J27a < 0 ? 0.0 : J27a
    ScaEff = (2.35 / CaConc)^p.ScaEffGam
    T72 = 90 * ScaEff
    T73 = p.T71 * (C8 - T72)
    T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73))
    T75 = p.T70 * (0.85 * (1 + T74) + 0.15)
    T76 = p.T70 * (0.85 * (1 - T74) + 0.15)
    T47 = p.T46 * 0.88 * p.GFR
    J48a = 0.88 * p.GFR * C2 - T47
    J48 = J48a < 0 ? 0.0 : J48a
    J53 = p.T52 * PhosGut
    J54 = p.T49 * C2
    J56 = p.T55 * IntraPO
    E0PicOBkb = p.MultPicOBkb * p.Pic0
    EmaxPicOBkb = p.FracPic0kb * p.Pic0
    EC50PicOBparenKb = ((E0PicOBkb - EmaxPicOBkb) * p.TGFBact0^p.PicOBgamkb /
                        (E0PicOBkb - p.Pic0)) - p.TGFBact0^p.PicOBgamkb
    EC50PicOBkb = exp(log(EC50PicOBparenKb) / p.PicOBgamkb)
    PicOBkb = E0PicOBkb - (E0PicOBkb - EmaxPicOBkb) *
             TGFBact^p.PicOBgamkb /
             (TGFBact^p.PicOBgamkb + EC50PicOBkb^p.PicOBgamkb)
    PicOBkbEff = PicOBkb / p.Pic0
    E0RUNX2kbEff = p.E0RUNX2kbEffFACT * p.kb
    RUNX2 = BCL2 > 105 ? BCL2 - 90 : 10.0
    RUNkbMax = E0RUNX2kbEff * p.RUNkbMaxFact
    INparen = RUNkbMax * p.RUNX20^p.RUNkbGAM /
              (E0RUNX2kbEff - p.kb) - p.RUNX20^p.RUNkbGAM
    RUNkb50 = exp(log(INparen) / p.RUNkbGAM)
    RUNX2kbPrimeEff = RUNkbMax * RUNX2^p.RUNkbGAM /
                      (RUNX2^p.RUNkbGAM + RUNkb50^p.RUNkbGAM)
    kbprime = E0RUNX2kbEff * PicOBkbEff - RUNX2kbPrimeEff
    kbslow = kbprime * p.Frackb
    kbfast = (p.kb * p.OB0 + kbslow * OBfast0 - kbslow * p.OB0) / OBfast0
    Frackb2 = kbfast / kbprime
    T29 = (p.T28 * p.T0 - 0.17533 * p.T0) / 0.17533
    T31 = p.T28 * T / (T + T29)
    T83 = R / 0.5
    J40 = T31 * T * T83 / (T + p.T81) + p.T87 * T
    T85Rpart = R^p.T80 / (R^p.T80 + p.T81^p.T80)
    T85 = p.T77 * T85Rpart
    F11 = T85
    INparenCtriol = ((p.CtriolMax - p.CtriolMin) * C8^p.CtriolPTgam /
                     (p.CtriolMax - 1)) - C8^p.CtriolPTgam
    Ctriol50 = exp(log(INparenCtriol) / p.CtriolPTgam)
    CtriolPTeff = p.CtriolMax - (p.CtriolMax - p.CtriolMin) *
                  C8^p.CtriolPTgam /
                  (C8^p.CtriolPTgam + Ctriol50^p.CtriolPTgam)
    PTin = p.PTout * CtriolPTeff
    INparenCa = (p.T58 - p.T61) * 2.35^p.T59 /
                (p.T58 - 385) - 2.35^p.T59
    T60 = exp(log(INparenCa) / p.T59)
    FCTD = (S / 0.5) * PTmax
    T63 = p.T58 - (p.T58 - p.T61) * CaConc^p.T59 /
          (CaConc^p.T59 + T60^p.T59)
    EPTH = T63 * FCTD
    IPTH = 0.693 * SC + p.IPTHinf
    SPTH = EPTH + IPTH
    kout = p.T57 / 14
    TERIPK = TERISC * p.TERICL / p.TERIVC

    du[1] = SPTH - kout * PTH + TERIPK
    du[2] = (0.01 - S) * T76 - S * T75
    du[3] = PTin - p.PTout * PTmax
    du[4] = A + J14 + p.KAVITD * (p.FVITD * VDORAL) - p.T69 * B
    du[5] = p.IPTHint - 0.693 * SC
    du[6] = SE - p.T64 * A
    du[7] = J14 - J15 - J27 + J40
    du[8] = J41 - J42 - J48 + p.T52 * PhosGut - J54 + J56
    du[9] = p.OralCa * F11 - J40
    du[10] = T36 * (1 - R) - T37 * R
    du[11] = kHApIn * Osteoblast - kLShap * HAp
    du[12] = (bigDb / PicOB) * D * p.FracOBfast * Frackb2 - kbfast * OBfast
    du[13] = (bigDb / PicOB) * D * (1 - p.FracOBfast) * p.Frackb - kbslow * OBslow
    du[14] = p.OralPhos * p.F12 - J53
    du[15] = J54 - J56
    du[16] = kinOC2 - KLSoc * OC
    du[17] = ROBin - KPT * ROB1
    du[18] = kinTGF * (Osteoblast / p.OB0)^p.OBtgfGAM - koutTGFeqn
    du[19] = koutTGFeqn - koutTGFact * TGFBact
    du[20] = kinL - p.koutL * L - p.k1 * O * L + p.k2 * N - p.k3 * RNK * L + p.k4 * M
    du[21] = kinRNK * TGFBact^p.kinRNKgam - p.koutRNK * RNK - p.k3 * RNK * L + p.k4 * M
    du[22] = p.k3 * RNK * L - p.k4 * M
    du[23] = p.k1 * O * L - p.k2 * N
    du[24] = pO - p.k1 * O * L + p.k2 * N - p.kO * O
    du[25] = J15 - J14 + J14a - J15a
    du[26] = J15a - J14a
    du[27] = RX2Kin - RX2Kout * RX2
    du[28] = crebKin - p.crebKout * CREB
    du[29] = bcl2Kin - p.bcl2Kout * BCL2
    du[30] = -TERISC * p.TERICL / p.TERIVC
    du[31] = J27
    du[32] = -(VDORAL * p.FVITD) * p.KAVITD
    return nothing
end

function peterson_khurana_problem(p::PetersonKhuranaParams=PetersonKhuranaParams();
                                  u₀=peterson_khurana_initial_state(p),
                                  tspan=(0.0, 210.0 * 24))
    length(u₀) == 32 || throw(ArgumentError("u₀ must have 32 states"))
    p_eff = peterson_khurana_derived_params(p, u₀)
    return ODEProblem(peterson_khurana_rhs!, collect(u₀), tspan, p_eff)
end

function _peterson_event_times(tspan, dose_start, dose_stop)
    t₀, t₁ = tspan
    reset_first = ceil(Int, t₀ / 24) * 24
    reset_last = floor(Int, t₁ / 24) * 24
    resets = reset_first <= reset_last ? collect(Float64.(reset_first:24:reset_last)) : Float64[]
    dose_first = max(dose_start, ceil(Int, max(t₀, dose_start) / 24) * 24)
    dose_last = min(dose_stop, t₁)
    doses = dose_first <= dose_last ? collect(Float64.(dose_first:24:dose_last)) : Float64[]
    return sort(unique(vcat(resets, doses)))
end

const _PETERSON_POSITIVE_STATES = (1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13,
                                   14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                   24, 27, 28, 29)

@inline function _peterson_out_of_domain(u, p, t)
    return any(i -> !isfinite(u[i]) || u[i] <= 0, _PETERSON_POSITIVE_STATES)
end

"""Solve the local 32-state model with the R script's default dosing schedule."""
function simulate_peterson_khurana(p::PetersonKhuranaParams=PetersonKhuranaParams();
                                   teriparatide_dose_mcg=100.0,
                                   vitamin_d_dose_mcg=0.5000000002,
                                   dose_start=720.0,
                                   dose_stop=nothing,
                                   with_events=true,
                                   kwargs...)
    tspan = haskey(kwargs, :tspan) ? kwargs[:tspan] : (0.0, 210.0 * 24)
    dose_stop_eff = dose_stop === nothing ? tspan[2] : dose_stop
    teri_dose = teriparatide_dose_mcg * 1e6 / 9424.8
    vitamin_d_dose = vitamin_d_dose_mcg * 1e6 / 417.0
    callback = nothing
    if with_events
        event_times = _peterson_event_times(tspan, dose_start, dose_stop_eff)
        function affect!(integrator)
            t = integrator.t
            integrator.u[31] = 0.0
            if t >= dose_start - 1e-8 && t <= dose_stop_eff + 1e-8
                integrator.u[30] += teri_dose
                integrator.u[32] += vitamin_d_dose
            end
            return nothing
        end
        # deSolve reports the state immediately before an event at an event
        # time; retain that convention so daily UCAL is comparable to the R
        # output while still applying the event to the continued trajectory.
        callback = PresetTimeCallback(event_times, affect!; save_positions=(true, false))
    end
    prob = peterson_khurana_problem(p; kwargs...)
    solve_kwargs = (; abstol=1e-10, reltol=1e-10, saveat=0.5, dtmax=1.0,
                    isoutofdomain=_peterson_out_of_domain)
    # Disable automatic differentiation: the source equations contain
    # fractional powers, and trial stages can temporarily cross zero even
    # when the accepted trajectory remains in the real positive domain.
    solver = Rodas5P(autodiff=AutoFiniteDiff())
    sol = callback === nothing ? solve(prob, solver; solve_kwargs...) :
          solve(prob, solver; callback=callback, solve_kwargs...)
    SciMLBase.successful_retcode(sol) ||
        error("Peterson/Khurana solve failed: $(sol.retcode)")
    return sol
end

"""Parameters for the Peterson–Riggs nonlinear lumbar-spine BMD extension.

The extension is the single indirect BMD equation from Peterson & Riggs
(2012), driven by the existing model's BSAP proxy (`OBfast + OBslow`) and
CTx proxy (`OC`).
"""
struct PetersonKhuranaBMDParams{N<:NamedTuple}
    values::N
end

Base.getproperty(p::PetersonKhuranaBMDParams, name::Symbol) =
    name === :values ? getfield(p, :values) : getproperty(getfield(p, :values), name)

function PetersonKhuranaBMDParams(p::PetersonKhuranaParams, u₀;
                                  BMD0=100.0, gamma_OB=0.0739,
                                  gamma_OC=0.0679, kout_BMD=0.000145)
    length(u₀) == 32 || throw(ArgumentError("the BMD extension requires the 32-state base model"))
    return PetersonKhuranaBMDParams(merge(p.values, (
        BMD0=BMD0,
        BSAP0=u₀[12] + u₀[13],
        CTx0=u₀[16],
        gamma_OB=gamma_OB,
        gamma_OC=gamma_OC,
        kout_BMD=kout_BMD,
        kin_BMD=kout_BMD * BMD0,
    )))
end

function peterson_khurana_bmd_initial_state(p::PetersonKhuranaParams=PetersonKhuranaParams();
                                            BMD0=100.0)
    return vcat(peterson_khurana_initial_state(p), BMD0)
end

function peterson_khurana_bmd_rhs!(du, u, p::PetersonKhuranaBMDParams, t)
    base_p = PetersonKhuranaParams(p.values)
    peterson_khurana_rhs!(@view(du[1:32]), @view(u[1:32]), base_p, t)
    BSAP = u[12] + u[13]
    CTx = u[16]
    du[33] = p.kin_BMD * (BSAP / p.BSAP0)^p.gamma_OB -
             p.kout_BMD * (CTx / p.CTx0)^p.gamma_OC * u[33]
    return nothing
end

"""Create the 33-state Peterson/Khurana model with the 2012 BMD output."""
function peterson_khurana_bmd_problem(p::PetersonKhuranaParams=PetersonKhuranaParams();
                                      u₀=peterson_khurana_bmd_initial_state(p),
                                      tspan=(0.0, 210.0 * 24),
                                      BMD0=100.0, gamma_OB=0.0739,
                                      gamma_OC=0.0679, kout_BMD=0.000145)
    length(u₀) == 33 || throw(ArgumentError("u₀ must have 33 states"))
    base_u₀ = collect(@view u₀[1:32])
    p_eff = peterson_khurana_derived_params(p, base_u₀)
    p_bmd = PetersonKhuranaBMDParams(p_eff, base_u₀;
                                     BMD0, gamma_OB, gamma_OC, kout_BMD)
    return ODEProblem(peterson_khurana_bmd_rhs!, collect(u₀), tspan, p_bmd)
end

"""Solve the local integrated model with the Peterson–Riggs BMD extension."""
function simulate_peterson_khurana_bmd(p::PetersonKhuranaParams=PetersonKhuranaParams();
                                       teriparatide_dose_mcg=100.0,
                                       vitamin_d_dose_mcg=0.5000000002,
                                       dose_start=720.0,
                                       dose_stop=nothing,
                                       with_events=true,
                                       BMD0=100.0,
                                       gamma_OB=0.0739,
                                       gamma_OC=0.0679,
                                       kout_BMD=0.000145,
                                       kwargs...)
    tspan = haskey(kwargs, :tspan) ? kwargs[:tspan] : (0.0, 210.0 * 24)
    dose_stop_eff = dose_stop === nothing ? tspan[2] : dose_stop
    teri_dose = teriparatide_dose_mcg * 1e6 / 9424.8
    vitamin_d_dose = vitamin_d_dose_mcg * 1e6 / 417.0
    u₀ = haskey(kwargs, :u₀) ? kwargs[:u₀] : peterson_khurana_bmd_initial_state(p; BMD0)
    callback = nothing
    if with_events
        event_times = _peterson_event_times(tspan, dose_start, dose_stop_eff)
        function affect!(integrator)
            t = integrator.t
            integrator.u[31] = 0.0
            if t >= dose_start - 1e-8 && t <= dose_stop_eff + 1e-8
                integrator.u[30] += teri_dose
                integrator.u[32] += vitamin_d_dose
            end
            return nothing
        end
        callback = PresetTimeCallback(event_times, affect!; save_positions=(true, false))
    end
    prob = peterson_khurana_bmd_problem(p; u₀, tspan, BMD0, gamma_OB, gamma_OC, kout_BMD)
    solve_kwargs = (; abstol=1e-10, reltol=1e-10, saveat=0.5, dtmax=1.0,
                    isoutofdomain=(u, p, t) -> any(i -> !isfinite(u[i]) || u[i] <= 0,
                                                   (1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13,
                                                    14, 15, 16, 17, 18, 19, 20, 21, 22,
                                                    23, 24, 27, 28, 29, 33)))
    solver = Rodas5P(autodiff=AutoFiniteDiff())
    sol = callback === nothing ? solve(prob, solver; solve_kwargs...) :
          solve(prob, solver; callback=callback, solve_kwargs...)
    SciMLBase.successful_retcode(sol) ||
        error("Peterson/Khurana BMD solve failed: $(sol.retcode)")
    return sol
end
