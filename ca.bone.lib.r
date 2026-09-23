

ca.bone.load.model <- function() {
  ## Sources R files to load functions for
  ## generating lists of compartments, initials,
  ## parameters, and the model function
  ##
  ## Returns: list with the following items:
  ##   $init named list of initial values
  ##   $param named list of model parameters
  ##   $cmt ordered vector of compartment names
  ## Also, the following functions are made available:
  ##   ca.bone.init()
  ##   ca.bone.cmt()
  ##   ca.bone.param()
  
  source("ca.bone.init.r")
  source("ca.bone.param.r")
  source("ca.bone.cmt.r",)
  source("ca.bone.model.r",)

  ca.bone.derive.init(list(init=ca.bone.init(),
                           param=ca.bone.param(),
                           cmt=ca.bone.cmt(),
                           model=ca.bone.model
                           )
                      )
}




ca.bone.derive.init <- function(camod) {
  ## Derive some initial values from
  ##   other initials and parameters.
  ## Returns: updated camod list after
  ##   calling copy.init()

  camod$init$A <- camod$init$B/10
  
  camod$init$TGFB <- camod$param$Pic0*1000
  camod$init$TGFBact <- camod$param$Pic0
  camod$init$OBfast <- camod$init$OB*camod$param$FracOBfast
  camod$init$OBslow <- camod$init$OB*(1-camod$param$FracOBfast)
  camod$init$M <- camod$param$k3*camod$init$RNK*camod$init$L/camod$param$k4
  camod$init$N <- camod$param$k1*camod$init$O*camod$init$L/camod$param$k2
  
  ca.bone.copy.init(camod)
}


ca.bone.copy.init <- function(camod) {

  ## Copy some initial compartment values
  ##   into the parameter list.
  ## Returns: updated camod list
  
  camod$param$Q0 <- camod$init$Q
  camod$param$OC0 <- camod$init$OC
  camod$param$RNK0 <- camod$init$RNK
  camod$param$RANKL0 <- camod$init$L
  camod$param$RNKL0 <- camod$init$L
  camod$param$RNK0 <- camod$init$RNK
  camod$param$OB0 <- camod$init$OB
  camod$param$ROB0 <- camod$init$ROB1
  camod$param$QboneInit <- camod$init$Qbone
  camod$param$OPG0 <- camod$init$O
  camod$param$RX20 <- camod$init$RX2
  camod$param$CREB0 <- camod$init$CREB
  camod$param$M0 <- camod$init$M
  camod$param$TGFBact0 <- camod$init$TGFBact
  camod$param$TGFB0 <- camod$init$TGFB
   
  camod
  
}


ca.bone.responses <- function(out,camod) {

  ## Derives some outcomes of interest from raw simulated output
  ##
  ## Arguments: data frame output from lsoda() call and model list
  ## Returns: updated simulation output with some markers derived
  ##          and some amounts of interest expressed as concentrations
  
  ## BSAP = bone-specific alkaline phosphatase
  ## sCTx = serum C-terminal telopeptide of type I collagen
  ## PTHconc = plasma pth concentration (pg/mL)
  ## PTHpM = plasma pth concentration (pM)
  ## CaConc = plasma calcium concentration (mM)
  ## calcitriol = plasma calcitriol (pM)
  
  out$PTHpM <- out$PTH/camod$param$V1
  out$PTHconc <- out$PTHpM*9.4
  out$BSAP <- with(out,OBfast+OBslow)
  out$sCTx <- out$OC
  out$CaConc <- out$P/camod$param$V1
  out$CalcitriolConc <- out$B/camod$param$V1
  
  out$UCALC <- out$UCAL
  
  out

}




