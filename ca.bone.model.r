
ca.bone.model <- function(t,y,p) {
  
  PTH <- y[1]
  S <- y[2]
  PTmax <- y[3]
  B <- y[4]
  SC <- y[5]
  A  <- y[6]
  P <- y[7]
  ECCPhos <- y[8]
  T <- y[9]
  R <- y[10]
  HAp <- y[11]
  OBfast <- y[12]
  OBslow <- y[13]
  PhosGut <- y[14]
  IntraPO <- y[15]
  OC <- y[16]
  ROB1 <- y[17]
  TGFB <- y[18]
  TGFBact <- y[19]
  L <- y[20]
  RNK <- y[21]
  M <- y[22]   
  N <- y[23]
  O <- y[24]
  Q <- y[25]
  Qbone <- y[26]
  RX2 <- y[27]
  CREB <- y[28]
  BCL2 <- y[29]
  TERISC <- y[30]
  UCAL <- y[31]
  VDORAL <- y[32]

  yn <- c()
  with(p, {
    
    T13 = (CaDay/24)/Q0      
    
    T15 = CaDay/(2.35*14*24)
    
    T17 = 3.85*T16 - 3.85
       
    Osteoclast =  OC
    
    J14OC50= exp(log((J14OCmax*OC0^J14OCgam/T13) - OC0^J14OCgam)/J14OCgam)
    
    OCeqn = (J14OCmax*Osteoclast^J14OCgam)/(Osteoclast^J14OCgam + J14OC50^J14OCgam)
    
    kinRNK = (koutRNK*RNK0 + k3*RNK0*RANKL0 - k4*M0) / TGFBact0^kinRNKgam 
    
    MOCratio = M/Osteoclast
    
    MOCratio0 = M0/OC0
    
    MOCratioEff = (MOCratio/MOCratio0)^MOCratioGam       
    
    J14OCdepend = OCeqn*Q0*FracJ14*MOCratioEff                 
    
    J14 = T13*Q0*(1-FracJ14) + J14OCdepend
    
    J41 = 0.464*J14  
    
    PicOCkin = Pic0
    
    bigDb = kb*OB0*Pic0/ROB0
       
    kinTGF = koutTGF0*TGFB0
    
    koutTGF = koutTGF0
    
    koutTGFact = koutTGF0*1000
    
    koutTGFeqn = koutTGF*TGFB*((Osteoclast/OC0)^OCtgfGAM)  	  
    
    E0PicROB = FracPicROB*Pic0
    
    EC50PicROBparen= (EmaxPicROB*TGFBact0^PicROBgam / (Pic0 - E0PicROB)) - TGFBact0^PicROBgam
    
    EC50PicROB = exp(log(EC50PicROBparen)/PicROBgam) 
    
    Dr = kb*OB0/Pic0      	  
    
    PicROB = E0PicROB + EmaxPicROB*TGFBact^PicROBgam/(TGFBact^PicROBgam + EC50PicROB^PicROBgam)
    
    ROBin2 = Dr*PicROB
    
    ROBin = ROBin2
    
    E0PicOB = FracPicOB*Pic0
    
    
    EC50PicOBparen = (EmaxPicOB*TGFBact0^PicOBgam/(Pic0 - E0PicOB)) - TGFBact0^PicOBgam
    
    EC50PicOB = exp(log(EC50PicOBparen)/PicOBgam)
    
    
    PicOB = E0PicOB + EmaxPicOB*TGFBact^PicOBgam / (TGFBact^PicOBgam + EC50PicOB^PicOBgam)
    
    KPT =1*(bigDb/PicOB)  
    
    
    D = ROB1
    
    EC50MeffOC = exp(log(M0^kinOCgam*EmaxMeffOC/(1-E0Meff) - M0^kinOCgam)/kinOCgam)
    
    
    MeffOC = E0Meff + (EmaxMeffOC * M^kinOCgam/(M^kinOCgam + EC50MeffOC^kinOCgam))
    
    kinOC2 = Da*PicOCkin*MeffOC*OC0
    
    
    E0PicOC = FracPicOC*Pic0
    
    EC50PicOCparen = (EmaxPicOC*TGFBact0^PicOCgam/(Pic0 - E0PicOC)) - TGFBact0^PicOCgam
    
    
    EC50PicOC = exp(log(EC50PicOCparen)/PicOCgam)
       
    PicOC = E0PicOC + ((EmaxPicOC*TGFBact^PicOCgam)/(TGFBact^PicOCgam + EC50PicOC^PicOCgam))
    
    PiL0 = (k3/k4)*RANKL0
    
    PiL = M/10
    
    EC50survInPar = (E0RANKL - EmaxL)*(PiL0^LsurvOCgam/(E0RANKL - 1)) - PiL0^LsurvOCgam
    
    
    EC50surv = exp(log(EC50survInPar)/LsurvOCgam)
    
    LsurvOC = E0RANKL - (E0RANKL - EmaxL)*(PiL^LsurvOCgam/(PiL^LsurvOCgam + EC50surv^LsurvOCgam))
       
       
    KLSoc = Da*PicOC*LsurvOC#*3
    
    C4 = PTH/V1
    
    T66 = (T67^AlphOHgam + 3.85^ AlphOHgam )/3.85^ AlphOHgam 
    
    k15a = k14a*QboneInit/Q0 
    
    J14a = k14a*Qbone
    
    J15a = k15a*Q 
    
    
    kLShap = 1/HApMRT
    
    kHApIn = kLShap/OB0

    J15 = T15*P*(1-FracJ15) + T15*P*FracJ15*HAp
    
    J42 = 0.464*J15 
    
    OBfast0 = OB0*FracOBfast
    
    
    Osteoblast = OBfast + OBslow
    
    kinLbase = koutL*RANKL0
       
       
    OsteoEffect = (Osteoblast/OB0)^OsteoEffectGam 
    
    PTH50 = EmaxLpth*3.85 - 3.85 
    
    
    PTHconc = C4
    
    LpthEff = EmaxLpth*(PTHconc) / ((PTH50*OsteoEffect^TESTPOWER) + (PTHconc))  
    
    kinL = kinLbase*(OsteoEffect)*LpthEff
    
    pObase = kO*OPG0
    
    
    pO = pObase*(D/ROB0)*((PTHconc+(opgPTH50*(D/ROB0)))/(2*PTHconc))+ IO
    
    RX2Kin = RX2Kout0*RX20
    
    
    EC50PTHRX2x = ((EmaxPTHRX2x*3.85)/(RX2Kout0 - E0rx2Kout)) - 3.85
    
    RX2Kout = E0rx2Kout + EmaxPTHRX2x*PTHconc/(PTHconc+EC50PTHRX2x)
    
    
    EC50PTHcreb = ((EmaxPTHcreb*3.85)/(1-E0crebKin)) -  3.85
    
    crebKin0= crebKout*CREB0
    
    
    crebKin =crebKin0* (E0crebKin + EmaxPTHcreb*PTHconc/(PTHconc+EC50PTHcreb))
    
    bcl2Kin = RX2*CREB*0.693
    
    
    CaConc = P/14
    
    C2 = ECCPhos/V1 
    
    
    PO4inhPTH = (C2/1.2)^PO4inhPTHgam      
    
    PhosEffTop = (PhosEff0 - 1)*( 1.2^PhosEffGam + PhosEff50^PhosEffGam ) 
    
    
    PhosEffBot = PhosEff0 * 1.2^PhosEffGam 
    
    PhosEffMax =  PhosEffTop / PhosEffBot
    
    
    PhosEff = PhosEff0 - (PhosEffMax*PhosEff0 * C2^PhosEffGam /(C2^ PhosEffGam  + PhosEff50^PhosEffGam))
    
    
    if(C2 > 1.2) { PhosEffect = PhosEff } else { PhosEffect = 1 }
    
    T68 = T66*C4^AlphOHgam/(T67^AlphOHgam*PO4inhPTH+C4^AlphOHgam) 
    
    
    SE = T65*T68*PhosEffect
    
    C8 = B/V1

    C1 = P/V1 
    
    T36 = T33 + (T34-T33)*(C8^CaPOgam/(T35^CaPOgam+ C8^CaPOgam))
    
    
    T37 = T34 - (T34-T33)*(C8^CaPOgam/(T35^CaPOgam+ C8^CaPOgam))
    
    CaFilt = 0.6*0.5*GFR*C1 
        
    ReabsMax = (0.3*GFR*2.35 - 0.149997)*(Reabs50 + 2.35) / 2.35
    
    ReabsPTHeff = (T16*C4)/(C4 + T17) 
    
    CaReabsActive = (ReabsMax*C1/(Reabs50 + C1))*ReabsPTHeff 
    
    T20 = CaFilt - CaReabsActive
    
    T10 = T7*C8/(C8+T9)
    
    J27a = (2-T10)*T20
    
    if(J27a<0) {J27 = 0} else{J27 = J27a}
    
    ScaEff =  (2.35/CaConc)^ ScaEffGam 
    
    T72 = 90 * ScaEff  
    
    T73 = T71 * (C8 - T72) 
    
    
    T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73))
    
    T75 = T70 * (0.85 * (1 + T74) + 0.15)  
    
    
    T76 = T70 * (0.85 * (1 - T74) + 0.15)
    
    T47 = T46*0.88*GFR
    
    
    J48a = 0.88*GFR*C2 - T47 
    
    if(J48a < 0){J48 = 0} else{J48 = J48a}
    
    
    J53 = T52*PhosGut
    
    J54 = T49*C2
    
    
    J56 = T55*IntraPO
    
    E0PicOBkb = MultPicOBkb*Pic0
    
    
    EmaxPicOBkb = FracPic0kb*Pic0
    
    EC50PicOBparenKb = ((E0PicOBkb - EmaxPicOBkb)*TGFBact0^PicOBgamkb) / (E0PicOBkb - Pic0)  - TGFBact0^PicOBgamkb
    
    
    EC50PicOBkb = exp(log(EC50PicOBparenKb)/PicOBgamkb)
    
    PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb)*TGFBact^PicOBgamkb / (TGFBact^PicOBgamkb + EC50PicOBkb^PicOBgamkb)
    
    PicOBkbEff = PicOBkb/Pic0
    
    E0RUNX2kbEff= E0RUNX2kbEffFACT*kb  
    
    
    if(BCL2 > 105) {RUNX2 = BCL2 - 90}else {RUNX2 = 10}
    
    RUNkbMax = E0RUNX2kbEff*RUNkbMaxFact
    
    
    INparen = (RUNkbMax * RUNX20^RUNkbGAM) / (E0RUNX2kbEff - kb) - RUNX20^RUNkbGAM
    
    RUNkb50 = exp(log(INparen)/RUNkbGAM)
    
    
    RUNX2kbPrimeEff = RUNkbMax*RUNX2^RUNkbGAM / (RUNX2^RUNkbGAM + RUNkb50^RUNkbGAM)
    
    kbprime = E0RUNX2kbEff*PicOBkbEff - RUNX2kbPrimeEff
    
    
    kbslow = kbprime*Frackb
    
    kbfast = (kb*OB0 + kbslow*OBfast0 - kbslow*OB0) / OBfast0 
    
    
    Frackb2 = kbfast/kbprime
    
    T29 = (T28*T0 - 0.17533*T0)/0.17533
    
    T31 = T28*T/(T+T29)	
    
    T83 = R/0.5
    
    J40 = T31*T*T83/(T + T81) + T87*T
    
    T85Rpart = R^T80/(R^T80 + T81^T80) 
    
    T85 = T77*T85Rpart
    
    F11 = T85
    
    INparenCtriol =((CtriolMax - CtriolMin) * C8^CtriolPTgam) / (CtriolMax - 1) - C8^CtriolPTgam 
    
    Ctriol50 = exp(log(INparenCtriol) / CtriolPTgam) 
    
    CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * C8^CtriolPTgam / (C8^CtriolPTgam + Ctriol50^CtriolPTgam) 
    
    PTin= PTout * CtriolPTeff
    
    INparenCa = (T58 - T61) * 2.35^T59 / (T58 - 385) - 2.35^T59  
    
    T60 = exp(log(INparenCa) / T59)  
    
    FCTD = (S / 0.5) * PTmax
    
    T63 =  T58 - (T58 - T61) * (CaConc)^T59 / ((CaConc)^T59 + T60^T59)
    
    EPTH = T63 * FCTD
    
    IPTH = 0.693*SC +  IPTHinf
    
    SPTH = EPTH + IPTH 
    
    kout = T57/14


    TERIPK = TERISC*TERICL/TERIVC
    
    ############################################
    ## Differential Equations
    ############################################
    
   
    ## Parathyroid (PTH)
    yn[1] = SPTH - kout*PTH + TERIPK
    
    ##  S - PT gland pool 
    #For 99% loss in PTH Gland Function##
    #yn[2] = (0.01 - S) * T76 - (S* T75) 
    #For 50% loss in PTH Gland Function##
    #yn[2] = (0.5 - S) * T76 - (S* T75) 
    
    yn[2] = (0.01 - S) * T76 - (S* T75)  
    
    ## PTmax - PT maximum capacity
    yn[3] = PTin - PTout * PTmax       
    
    ## B  -  Plasma calcitriol unclear what authors mean by D14 changed to J14
    yn[4] = A + J14 + KAVITD*(FVITD*VDORAL) - T69 * B
    
    ## SC - SubQ
    yn[5] = IPTHint - 0.693*SC
    
    ## A - 1-alpha hydroxylase
    yn[6] = SE - T64*A  
    
    ## P - plasma calcium
    yn[7] = J14 - J15 - J27 + J40   
    
    ## ECCphos - phosphate
    yn[8] = J41  - J42 - J48 + J53 - J54 + J56
    
    ## T - oral calcium gut
    yn[9] = OralCa*F11 - J40
    
    ## R - intestinal calcium
    yn[10] = T36*(1- R) - T37*R
    
    ## HAp - Hydroxyapatite
    yn[11] = kHApIn*Osteoblast  - kLShap*HAp
    
    ## OBfast - Fast osteoblast
    yn[12] = (bigDb/PicOB)*D*FracOBfast*Frackb2  - kbfast*OBfast
    
    ## OBslow - Slow osteoblast
    yn[13] = (bigDb/PicOB)*D*(1-FracOBfast)*Frackb - kbslow*OBslow
    
    ## PhosGut - oral phosphate
    yn[14] = OralPhos *F12 - J53
    
    ## IntraPO - intraceulular phosphate
    yn[15] = J54 - J56
    
    ## OC = Osteoclast
    yn[16] = kinOC2 - KLSoc*OC
    
    ## ROB1 - responding osteoblast
    yn[17] = ROBin - KPT*ROB1
    
    ## TGFB - Latent TGFbeta
    yn[18] = kinTGF*(((Osteoblast/OB0)^OBtgfGAM)) - koutTGFeqn
    
    ## TGFBact - active TGFbeta
    yn[19] = koutTGFeqn - koutTGFact*TGFBact
    
    ## L - RANKL
    yn[20] = kinL- koutL*L - k1*O*L + k2*N - k3*RNK*L + k4*M
    
    ## RNK - RANK
    yn[21] = kinRNK*TGFBact^kinRNKgam - koutRNK*RNK - k3*RNK*L  + k4*M 
    
    ## M - RANK-RANKL
    yn[22] = k3*RNK*L - k4*M 
    
    ## N - RANKL-OPG
    yn[23] = k1*O*L - k2*N 
    
    ## O - OPG
    yn[24] = pO - k1*O*L + k2*N - kO*O 
    
    ## Q - Bone Ca (IC)
    yn[25] = J15 - J14 + J14a - J15a 
    
    ## Qbone - Bone Ca (non-IC)
    yn[26] = J15a - J14a 
    
    ## RX2 - RunX2
    yn[27] = RX2Kin - RX2Kout*RX2 
    
    ## CREB
    yn[28] = crebKin - crebKout*CREB
    
    ## BCL2 - Bcl-2
    yn[29] = bcl2Kin - bcl2Kout*BCL2

    ## Teriparatide SQ dosing compartment
    yn[30] = - TERISC*TERICL/TERIVC
    
    ## UCAL
    yn[31] = J27
    
    ## VDORAL (Linked to # B - plasma clacitriol)
    yn[32] = -(VDORAL*FVITD)*KAVITD

    list(c(yn))
  })
  
}




