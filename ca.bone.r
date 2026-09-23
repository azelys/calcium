# updated based on fda's publication

library(deSolve)
library(lattice)
library(reshape)

source("ca.bone.lib.r")

camod <- ca.bone.load.model()
camod <- ca.bone.derive.init(camod)


## EVALUATION TIMES IN HOURS
## seq (start time, duration of simulation in days*24 hour, interval)
times <- seq(0,210*24,0.5)


## TERIPARATIDE DOSING EVENTS (TIMES IN HOURS)
#Description: seq (start time, duration of simulation in days*24 hour, interval)
teri.times <- seq(720,210*24,24)
teri.dose.mcg <- 100.00000
teri.dose <- teri.dose.mcg*1E6/9424.8


## New feature of the model to test/account for oral vitamin D supplement
vitd.times <- seq(720,210*24,24)
vitd.dose.mcg <- 0.5000000002
vitd.dose <- vitd.dose.mcg*1E6/417

## New feature of the model to compute and 24 hour urinary calcium excretion
events1 <- data.frame(
      var="TERISC", 
      time=teri.times,  
      value=teri.dose,  
      method="add"
    )
events2 <- data.frame(
      var="VDORAL", 
      time=vitd.times,  
      value=vitd.dose,  
      method="add"
    )
events3 <- data.frame(
      var="UCAL",
      time=times[times%%24==0],
      value=rep(0,length(times[times%%24==0])),
      method="rep"
    )

events <- rbind(events1,events2,events3)
events <- events[order(events$time),]


## ADD DOSING EVENTS TO EVALUATION TIMES
times <- sort(unique(c(times,events$time)))


## RUN THE MODEL

out <- as.data.frame(
                     lsoda(
                           unlist(camod$init[camod$cmt]), 
                           times, 
                           camod$model,
                           camod$param,
                           rtol=1e-10,
                           atol=1E-10,
                           ynames=F,
                           events=list(data=events)
                           )
                     )


## POST PROCESSING
out <- ca.bone.responses(out,camod)
out$UCALQD<-ifelse(out$time%%24==0,out$UCAL,0)


## Optional plotting routine
## REQUIRES loading of reshape and lattice libraries
print(xyplot(value~time/24|variable,
       data=melt(out,measure.vars=c("BSAP","sCTx","PTHpM","PTHconc",camod$cmt),id.vars="time"),
       type='l',par.strip.text=list(cex=0.8),
       scales=list(y=list(relation='free')),
       xlab="Time (days)",
       ylab="DV"
       )
)


