## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
#install.packages("FME")
library(FME)         # Package for MCMC simulation and model fitting
#install.packages("minpack.lm")
library(minpack.lm)  # Package for model fitting
#install.packages("reshape")
library(reshape)     # Package for melt function to reshape the table
#install.packages("truncnorm")
library(truncnorm)   # Package for the truncated normal distribution function 
#install.packages("EnvStats")
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
#install.packages("invgamma")
library(invgamma)    # Package for inverse gamma distribution function
#install.packages("foreach")
library(foreach)     # Package for parallel computing
#install.packages("doParalle")
library(doParallel)  # Package for parallel computing
#install.packages("bayesplot")
library(bayesplot)   # Package for MCMC traceplot
#rm(list=ls())
## Build mrgsolve-based PK Model
mod <- mcode ("pk", HumanPK.code)


############################################# Model Calibration with levenberg-marquart #############

#Plasma.DEHP <- read.csv(file="~/Dropbox/R/DEHP/2012_Plasma_DEHP.csv")
############ plasma #####################
Cal1 <- read.csv(file="~/Dropbox/R/DPHP/Cal1_plasma.csv")
Cal1_plasma_DPHP <- cbind.data.frame (Time = Cal1$time,
                                      Plasma_DPHP = Cal1$DPHP)
Cal1_plasma_MPHP <- cbind.data.frame (Time = Cal1$time,
                                      Plasma_MPHP = Cal1$MPHP)
Cal1_plasma_OH   <- cbind.data.frame (Time = Cal1$time,
                                      Plasma_OH = Cal1$OH)
Cal1_plasma_oxo  <- cbind.data.frame (Time = Cal1$time,
                                      Plasma_oxo = Cal1$oxo)



Cal2 <- read.csv(file="~/Dropbox/R/DPHP/Cal2_plasma.csv")
Cal2_plasma_DPHP <- cbind.data.frame (Time = Cal2$time,
                                      Plasma_DPHP = Cal2$DPHP)
Cal2_plasma_MPHP <- cbind.data.frame (Time = Cal2$time,
                                      Plasma_MPHP = Cal2$MPHP)
Cal2_plasma_OH   <- cbind.data.frame (Time = Cal2$time,
                                      Plasma_OH = Cal2$OH)
Cal2_plasma_oxo  <- cbind.data.frame (Time = Cal2$time,
                                      Plasma_oxo = Cal2$OXO)

BW_1         = 83
DOSEoral1    = 717 * BW_1 
BW_2          = 75
DOSEoral2     = 639 * BW_2  # ug Oral dose; 

######### urine ##################
Cal1a <- read.csv(file="~/Dropbox/R/DPHP/Cal1_urine.csv")
Cal1_urine_MPHP <- cbind.data.frame (Time = Cal1a$time,
                                     Urine_MPHP = Cal1a$MPHP/100*DOSEoral1)
Cal1_urine_OH <- cbind.data.frame (Time   = Cal1a$time,
                                   Urine_OH = Cal1a$OH/100*DOSEoral1)
Cal1_urine_oxo <- cbind.data.frame (Time  = Cal1a$time,
                                    Urine_oxo = Cal1a$oxo/100*DOSEoral1)


Cal2a <- read.csv(file="~/Dropbox/R/DPHP/Cal2_urine.csv")
Cal2_urine_MPHP <- cbind.data.frame (Time = Cal2a$time,
                                     Urine_MPHP = Cal2a$MPHP/100*DOSEoral2)
Cal2_urine_OH <- cbind.data.frame (Time = Cal2a$time,
                                   Urine_OH = Cal2a$OH/100*DOSEoral2)
Cal2_urine_oxo <- cbind.data.frame (Time = Cal2a$time,
                                    Urine_oxo = Cal2a$OXO/100*DOSEoral2)


# prediction function
Pred <- function(pars, pred=FALSE) {
  
  ##' Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  
  ### shared parameter
  tinterval   = 24
  TDoses      = 1
  
  #1
  # Exposure scenario #Inhalation_dermal_cal_1
  BW_1         = 83
  DOSEoral1    = 717 * BW_1   # ug Oral dose; 
  Cal_1    <- ev(ID=1, amt= DOSEoral1, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_cal_1  <- Cal_1 
  
  #2
  # Exposure scenario #Inhalation_dermal_cal_2
  BW_2          = 75
  DOSEoral2     = 639 * BW_2  # ug Oral dose; 
  Cal_2    <- ev(ID=1, amt= DOSEoral2, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_cal_2  <- Cal_2 
  
  
  # set up the exposure time
  tsamp=tgrid(0,48,0.01)  # 48 hours
  
  out.1 <- 
    
    mod %>% 
    param(pars, BW = 83)%>% 
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH, Urine_oxo) %>%
    update(atol = 1E-8,maxsteps = 100000000) %>%
    mrgsim_d(data = ex_cal_1, tgrid=tsamp)
  out.1<-cbind.data.frame (Time = out.1$time,
                           Plasma_oxo    = out.1$Plasma_oxo,
                           Plasma_DPHP   = out.1$Plasma_DPHP,
                           Plasma_MPHP   = out.1$Plasma_MPHP,
                           Plasma_OH     = out.1$Plasma_OH,
                           Urine_MPHP    = (out.1$Urine_MPHP),
                           Urine_OH      = (out.1$Urine_OH),
                           Urine_oxo     = (out.1$Urine_oxo))
  
  
  out.2 <-
    
    mod %>% 
    param(pars,BW = 75) %>%  
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000000) %>%
    mrgsim_d(data = ex_cal_2, tgrid=tsamp)
  out.2<-cbind.data.frame (Time = out.2$time,
                           Plasma_oxo    = out.2$Plasma_oxo,
                           Plasma_DPHP   = out.2$Plasma_DPHP,
                           Plasma_MPHP   = out.2$Plasma_MPHP,
                           Plasma_OH     = out.2$Plasma_OH,
                           Urine_MPHP    = (out.2$Urine_MPHP),
                           Urine_OH      = (out.2$Urine_OH),
                           Urine_oxo     = (out.2$Urine_oxo))
  
  
  return(list("out.1" = out.1,
              "out.2" = out.2))
  
}

## Cost fuction (from FME pckage) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out <-  Pred (pars)
  cost<- modCost(model=out$out.1, obs= Cal1_plasma_oxo,   x="Time")
  cost<- modCost(model=out$out.1, obs= Cal1_plasma_DPHP,  x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_plasma_MPHP,  x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_plasma_OH,    x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_urine_MPHP,   x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_urine_OH,     x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_urine_oxo,    x="Time", cost = cost)
  
  cost<- modCost(model=out$out.2, obs= Cal2_plasma_oxo,   x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_plasma_DPHP,  x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_plasma_MPHP,  x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_plasma_OH,    x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_urine_MPHP,   x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_urine_OH,     x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_urine_oxo,    x="Time", cost = cost)
  
  return(cost)
}



###################### Local sensitivity analysis #########################
## Choose the senstiive parameters in the model                          ##
## initial parmaeters######################################################

theta.int <- log(c(
  Vp                  =   3.2,    #: Volume of distribution parent compound L/kg bw
  Vm                  =   0.05,   # = Volume of distribution metabolite compound 
  ka                  =   0.27,     #= stomach absorption rate parent 1/h
  
  lambda_p            =   0.12,     #= lambda for parent 
  lambda_zm           =   0.17,    #= lambda for metabolite 
  lambda_u            =   2E-4,   #= lambda for metabolite urine 
  lambda_z_OH         =   0.27,       #= h-1
  lambda_z_oxo        =   0.16,      # = h-1
  
  lambda_u_OH         =   4.8e-3,     #= h-1
  lambda_u_oxo        =   4e-3,    # = h-1
  
  Fr_OH               =   0.5,    # = h-1
  Fr_oxo               =   0.17,  #  = h-1
  
  V_OH                =   0.5,   # = Volume of distribution parent compound L/kg bw
  V_oxo               =   0.5    #= Volume of distribution metabolite compound 
  
  
  
))

# Senstivity function (FME)
# Check the senstive parameters in the modelhttp://ec2-3-135-223-180.us-east-2.compute.amazonaws.com/graphics/plot_zoom_png?width=1200&height=868
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)
Sen=summary(SnsPlasma)
write.csv(Sen,file="~/Dropbox/NM/Class/DEHP/PK/DPHP/result/Sen.human.csv")

plot(summary(SnsPlasma))





######################### MCMC optimization ##########################
## set up senstivie or necessary parametesr as model input          ##
######################################################################
theta <- log(c(
  Vp                 =   1,    #: Volume of distribution parent compound L/kg bw
  Vm                  =   0.05,   # = Volume of distribution metabolite compound 
  ka                  =   3,     #= stomach absorption rate parent 1/h
  
  lambda_p            =   0.8,     #= lambda for parent 
  lambda_zm           =   1.7,    #= lambda for metabolite 
  lambda_u            =   5E-3,   #= lambda for metabolite urine 
  
  lambda_z_OH         =   0.7,       #= h-1
  lambda_z_oxo        =   0.16,      # = h-1
  
  lambda_u_OH         =   4e-2,     #= h-1
  lambda_u_oxo        =   1e-2,    # = h-1
  
  Fr_OH               =   0.5,    # = h-1
  Fr_oxo               =   0.17,  #  = h-1
  
  V_OH                =   0.5,   # = Volume of distribution parent compound L/kg bw
  V_oxo               =   0.5    #= Volume of distribution metabolite compound 
  
))

# Model fitting
#upper = c(Inf,Inf,Inf,Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
system.time(Fit<- modFit(f=MCcost, p=theta, method = "Nelder-Mead",
                         lower = c(-2,  -5,0.1,  -0.9,-1,-6.5,  -0.5,-2,  -10,  -7,    -5,  -5,   -5,-5),
                         upper = c( 1,0.01,  5,   0.1, 1.5,  -1,     2, 2,   -1,  -1,  0.01,0.01, 0.01, 0.01),
                         control = nls.lm.control(nprint=1)))

summary(Fit)
### para value
exp(Fit$par)

res=MCcost(Fit$par)$residuals$res
sum(res^2)


## Model calibration plot using ggplot2 
Sim.fit.A = Pred (Fit$par)$out.1         ## Time-course concentration profiles using estimated parameters under exposure senaior 1
Sim.fit.B = Pred (Fit$par)$out.2         ## Simulaiton of exposure scenaior 2


df.sim.A1  = cbind.data.frame (Time=Sim.fit.A$Time, Plasma_DPHP=Sim.fit.A$Plasma_DPHP)
df.sim.A2  = cbind.data.frame (Time=Sim.fit.A$Time, Plasma_MPHP=Sim.fit.A$Plasma_MPHP)
df.sim.A3  = cbind.data.frame (Time=Sim.fit.A$Time, Plasma_OH  =Sim.fit.A$Plasma_OH)
df.sim.A4  = cbind.data.frame (Time=Sim.fit.A$Time, Plasma_oxo =Sim.fit.A$Plasma_oxo)

df.sim.A6  = cbind.data.frame (Time=Sim.fit.A$Time, Urine_MPHP =Sim.fit.A$Urine_MPHP)
df.sim.A7  = cbind.data.frame (Time=Sim.fit.A$Time, Urine_OH   =Sim.fit.A$Urine_OH)
df.sim.A8  = cbind.data.frame (Time=Sim.fit.A$Time, Urine_oxo  =Sim.fit.A$Urine_oxo)

df.sim.B1  = cbind.data.frame (Time=Sim.fit.B$Time, Plasma_DPHP=Sim.fit.B$Plasma_DPHP)
df.sim.B2  = cbind.data.frame (Time=Sim.fit.B$Time, Plasma_MPHP=Sim.fit.B$Plasma_MPHP)
df.sim.B3  = cbind.data.frame (Time=Sim.fit.B$Time, Plasma_OH  =Sim.fit.B$Plasma_OH)
df.sim.B4  = cbind.data.frame (Time=Sim.fit.B$Time, Plasma_oxo =Sim.fit.B$Plasma_oxo)

df.sim.B6  = cbind.data.frame (Time=Sim.fit.B$Time, Urine_MPHP =Sim.fit.B$Urine_MPHP)
df.sim.B7  = cbind.data.frame (Time=Sim.fit.B$Time, Urine_OH   =Sim.fit.B$Urine_OH)
df.sim.B8  = cbind.data.frame (Time=Sim.fit.B$Time, Urine_oxo  =Sim.fit.B$Urine_oxo)



## Plot
plot.A1=
  ggplot() +
  geom_line (data = df.sim.A1,aes(Time, Plasma_DPHP), col="firebrick", lwd=2)+
  geom_point(data = Cal1_plasma_DPHP, aes(Time,Plasma_DPHP),size=2.5) + ylab("Concentration") 
plot.A1

plot.A2=
  ggplot() +
  geom_line (data = df.sim.A2,aes(Time,Plasma_MPHP), col="firebrick", lwd=2)+
  geom_point(data = Cal1_plasma_MPHP, aes(Time, Plasma_MPHP),size=2.5) + ylab("Concentration") 
plot.A2

plot.A3=
  ggplot() +
  geom_line (data = df.sim.A3,aes(Time,Plasma_OH), col="firebrick", lwd=2)+
  geom_point(data = Cal1_plasma_OH, aes(Time, Plasma_OH),size=2.5) + ylab("Concentration") 
plot.A3


plot.A4=
  ggplot() +
  geom_line (data = df.sim.A4,aes(Time,Plasma_oxo), col="firebrick", lwd=2)+
  geom_point(data = Cal1_plasma_oxo, aes(Time, Plasma_oxo),size=2.5) + ylab("Concentration") 
plot.A4


plot.B1=
  ggplot() +
  geom_line (data = df.sim.B1,aes(Time, Plasma_DPHP), col="firebrick", lwd=2)+
  geom_point(data = Cal2_plasma_DPHP, aes(Time,Plasma_DPHP),size=2.5) + ylab("Concentration") 
plot.B1

plot.B2=
  ggplot() +
  geom_line (data = df.sim.B2,aes(Time,Plasma_MPHP), col="firebrick", lwd=2)+
  geom_point(data = Cal2_plasma_MPHP, aes(Time, Plasma_MPHP),size=2.5) + ylab("Concentration") 
plot.B2

plot.B3=
  ggplot() +
  geom_line (data = df.sim.B3,aes(Time,Plasma_OH), col="firebrick", lwd=2)+
  geom_point(data = Cal2_plasma_OH, aes(Time, Plasma_OH),size=2.5) + ylab("Concentration") 
plot.B3


plot.B4=
  ggplot() +
  geom_line (data = df.sim.B4,aes(Time,Plasma_oxo), col="firebrick", lwd=2)+
  geom_point(data = Cal2_plasma_oxo, aes(Time, Plasma_oxo),size=2.5) + ylab("Concentration") 
plot.B4


###### URINE
plot.A1a=
  ggplot() +
  geom_line (data = df.sim.A6,aes(Time,Urine_MPHP), col="firebrick", lwd=2)+
  geom_point(data = Cal1_urine_MPHP,aes(Time, Urine_MPHP),size=2.5) + ylab("Concentration") 
plot.A1a

plot.A2a=
  ggplot() +
  geom_line (data = df.sim.A7,aes(Time,Urine_OH), col="firebrick", lwd=2)+
  geom_point(data = Cal1_urine_OH,aes(Time, Urine_OH),size=2.5) + ylab("Concentration") 
plot.A2a


plot.A3a=
  ggplot() +
  geom_line (data = df.sim.A8,aes(Time,Urine_oxo), col="firebrick", lwd=2)+
  geom_point(data = Cal1_urine_oxo,aes(Time, Urine_oxo),size=2.5) + ylab("Concentration") 
plot.A3a


plot.B1a=
  ggplot() +
  geom_line (data = df.sim.B6,aes(Time,Urine_MPHP), col="firebrick", lwd=2)+
  geom_point(data = Cal2_urine_MPHP,aes(Time, Urine_MPHP),size=2.5) + ylab("Concentration") 
plot.B1a

plot.B2a=
  ggplot() +
  geom_line (data = df.sim.B7,aes(Time,Urine_OH), col="firebrick", lwd=2)+
  geom_point(data = Cal2_urine_OH,aes(Time, Urine_OH),size=2.5) + ylab("Concentration") 
plot.B2a


plot.B3a=
  ggplot() +
  geom_line (data = df.sim.B8,aes(Time,Urine_oxo), col="firebrick", lwd=2)+
  geom_point(data = Cal2_urine_oxo,aes(Time, Urine_oxo),size=2.5) + ylab("Concentration") 
plot.B3a




############################################# Model Calibration with MCMC ###################################################
## EIGHT data sets was used in model evaluation                                                                               #                                                       #
#############################################################################################################################
############ plasma #####################
########### 1
############ plasma #####################
Opt1 <- read.csv(file="~/Dropbox/R/DPHP/Opt1_plasma.csv")
Opt1_plasma_DPHP <- cbind.data.frame (Time = Opt1$time,
                                      Plasma_DPHP = Opt1$DPHP)
Opt1_plasma_MPHP <- cbind.data.frame (Time = Opt1$time,
                                      Plasma_MPHP = Opt1$MPHP)
Opt1_plasma_OH   <- cbind.data.frame (Time = Opt1$time,
                                      Plasma_OH = Opt1$OH)
Opt1_plasma_oxo  <- cbind.data.frame (Time = Opt1$time,
                                      Plasma_oxo = Opt1$OXO)

Opt2 <- read.csv(file="~/Dropbox/R/DPHP/Opt2_plasma.csv")
Opt2_plasma_DPHP <- cbind.data.frame (Time = Opt2$time,
                                      Plasma_DPHP = Opt2$DPHP)
Opt2_plasma_MPHP <- cbind.data.frame (Time = Opt2$time,
                                      Plasma_MPHP = Opt2$MPHP)
Opt2_plasma_OH   <- cbind.data.frame (Time = Opt2$time,
                                      Plasma_OH = Opt2$OH)
Opt2_plasma_oxo  <- cbind.data.frame (Time = Opt2$time,
                                      Plasma_oxo = Opt2$OXO)

Opt3 <- read.csv(file="~/Dropbox/R/DPHP/Opt3_plasma.csv")
Opt3_plasma_DPHP <- cbind.data.frame (Time = Opt3$time,
                                      Plasma_DPHP = Opt3$DPHP)
Opt3_plasma_MPHP <- cbind.data.frame (Time = Opt3$time,
                                      Plasma_MPHP = Opt3$MPHP)
Opt3_plasma_OH   <- cbind.data.frame (Time = Opt3$time,
                                      Plasma_OH = Opt3$OH)
Opt3_plasma_oxo  <- cbind.data.frame (Time = Opt3$time,
                                      Plasma_oxo = Opt3$OXO)

Opt4 <- read.csv(file="~/Dropbox/R/DPHP/Opt4_plasma.csv")
Opt4_plasma_DPHP <- cbind.data.frame (Time = Opt4$time,
                                      Plasma_DPHP = Opt4$DPHP)
Opt4_plasma_MPHP <- cbind.data.frame (Time = Opt4$time,
                                      Plasma_MPHP = Opt4$MPHP)
Opt4_plasma_OH   <- cbind.data.frame (Time = Opt4$time,
                                      Plasma_OH = Opt4$OH)
Opt4_plasma_oxo  <- cbind.data.frame (Time = Opt4$time,
                                      Plasma_oxo = Opt4$OXO)

Opt5 <- read.csv(file="~/Dropbox/R/DPHP/Opt5_plasma.csv")
Opt5_plasma_DPHP <- cbind.data.frame (Time = Opt5$time,
                                      Plasma_DPHP = Opt5$DPHP)
Opt5_plasma_MPHP <- cbind.data.frame (Time = Opt5$time,
                                      Plasma_MPHP = Opt5$MPHP)
Opt5_plasma_OH   <- cbind.data.frame (Time = Opt5$time,
                                      Plasma_OH = Opt5$OH)
Opt5_plasma_oxo  <- cbind.data.frame (Time = Opt5$time,
                                      Plasma_oxo = Opt5$OXO)

Opt6 <- read.csv(file="~/Dropbox/R/DPHP/Opt6_plasma.csv")
Opt6_plasma_DPHP <- cbind.data.frame (Time = Opt6$time,
                                      Plasma_DPHP = Opt6$DPHP)
Opt6_plasma_MPHP <- cbind.data.frame (Time = Opt6$time,
                                      Plasma_MPHP = Opt6$MPHP)
Opt6_plasma_OH   <- cbind.data.frame (Time = Opt6$time,
                                      Plasma_OH = Opt6$OH)
Opt6_plasma_oxo  <- cbind.data.frame (Time = Opt6$time,
                                      Plasma_oxo = Opt6$OXO)




######### urine ##################
Opt1a <- read.csv(file="~/Dropbox/R/DPHP/Opt1_urine.csv")
Opt1_urine_MPHP <- cbind.data.frame (Time = Opt1a$time,
                                     Urine_MPHP = Opt1a$MPHP)
Opt1_urine_OH <- cbind.data.frame (Time   = Opt1a$time,
                                   Urine_OH = Opt1a$OH)
Opt1_urine_oxo <- cbind.data.frame (Time  = Opt1a$time,
                                    Urine_oxo = Opt1a$OXO)


Opt2a <- read.csv(file="~/Dropbox/R/DPHP/Opt2_urine.csv")
Opt2_urine_MPHP <- cbind.data.frame (Time = Opt2a$time,
                                     Urine_MPHP = Opt2a$MPHP)
Opt2_urine_OH <- cbind.data.frame (Time = Opt2a$time,
                                   Urine_OH = Opt2a$OH)
Opt2_urine_oxo <- cbind.data.frame (Time = Opt2a$time,
                                    Urine_oxo = Opt2a$OXO)


Opt3a <- read.csv(file="~/Dropbox/R/DPHP/Opt3_urine.csv")
Opt3_urine_MPHP <- cbind.data.frame (Time = Opt3a$time,
                                     Urine_MPHP = Opt3a$MPHP)
Opt3_urine_OH <- cbind.data.frame (Time   = Opt3a$time,
                                   Urine_OH = Opt3a$OH)
Opt3_urine_oxo <- cbind.data.frame (Time  = Opt3a$time,
                                    Urine_oxo = Opt3a$OXO)


Opt4a <- read.csv(file="~/Dropbox/R/DPHP/Opt4_urine.csv")
Opt4_urine_MPHP <- cbind.data.frame (Time = Opt4a$time,
                                     Urine_MPHP = Opt4a$MPHP)
Opt4_urine_OH <- cbind.data.frame (Time = Opt4a$time,
                                   Urine_OH = Opt4a$OH)
Opt4_urine_oxo <- cbind.data.frame (Time = Opt4a$time,
                                    Urine_oxo = Opt4a$OXO)


Opt5a <- read.csv(file="~/Dropbox/R/DPHP/Opt5_urine.csv")
Opt5_urine_MPHP <- cbind.data.frame (Time = Opt5a$time,
                                     Urine_MPHP = Opt5a$MPHP)
Opt5_urine_OH <- cbind.data.frame (Time = Opt5a$time,
                                   Urine_OH = Opt5a$OH)
Opt5_urine_oxo <- cbind.data.frame (Time = Opt5a$time,
                                    Urine_oxo = Opt5a$OXO)


Opt6a <- read.csv(file="~/Dropbox/R/DPHP/Opt6_urine.csv")
Opt6_urine_MPHP <- cbind.data.frame (Time = Opt6a$time,
                                     Urine_MPHP = Opt6a$MPHP)
Opt6_urine_OH <- cbind.data.frame (Time   = Opt6a$time,
                                   Urine_OH = Opt6a$OH)
Opt6_urine_oxo <- cbind.data.frame (Time  = Opt6a$time,
                                    Urine_oxo = Opt6a$OXO)



## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)

theta.MCMC <- log(c(
  Vp                  =   1.46246679 ,    #: Volume of distribution parent compound L/kg bw
  Vm                  =   0.78796963,   # = Volume of distribution metabolite compound 
  ka                  =   1.15901210,     #= stomach absorption rate parent 1/h
  
  lambda_p            =   1.02509013,     #= lambda for parent 
  lambda_zm           =   2.06176703,    #= lambda for metabolite 
  lambda_u            =   0.00199547,   #= lambda for metabolite urine 
  lambda_z_OH         =   0.79440136 ,       #= h-1
  lambda_z_oxo        =   0.19518004,      # = h-1
  
  lambda_u_OH         =   0.06556724,     #= h-1
  lambda_u_oxo        =   0.07657706,    # = h-1
  
  Fr_OH               =   0.28006279,    # = h-1
  Fr_oxo              =   0.68611923,  #  = h-1
  
  V_OH                =   0.66301511,   # = Volume of distribution parent compound L/kg bw
  V_oxo               =   0.86278857,    #= Volume of distribution metabolite compound 
  
  sig2                =   1, 
  
  sig_Vp                 = 1,               
  sig_Vm                 = 1,
  sig_ka                 = 1,
  
  sig_lambda_p           = 1,
  sig_lambda_zm          = 1,
  sig_lambda_u           = 1,
  sig_lambda_z_OH        =   1, 
  sig_lambda_z_oxo       =   1,
  
  sig_lambda_u_OH        =   1,     #= h-1
  sig_lambda_u_oxo       =   1,    # = h-1
  
  sig_Fr_OH              =   1,    # = h-1
  sig_Fr_oxo             =   1,  #  = h-1
  
  sig_V_OH               =   1,   # = Volume of distribution parent compound L/kg bw
  sig_V_oxo              =   1    #= 
  
  
))


which_sig <- grep("sig", names(theta.MCMC))

## Maximum likelihood estimation (MLE) fuction for MCMC
mcmc.fun <- function (pars, pred=FALSE){
  
  ## Get out of log domain
  pars.data <- lapply(pars [-which_sig],exp)
  names(pars.data) <- names(pars.data)
  
  
  ##### common parameter
  tinterval             = 24                            ## hr, Time interval
  TDoses                = 1                             ## Dose times
  
  # Exposure scenario for single low dose
  #1
  # Exposure scenario #Inhalation_dermal_cal_1
  BW_1          = 83
  DOSEoral1   = BW_1 * 717  # ug Oral dose; 
  
  Opt_1    <- ev(ID=1, amt= DOSEoral1, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_1  <- Opt_1 
  
  #2
  # Exposure scenario #Inhalation_dermal_Opt_2
  BW_2          = 75
  DOSEoral2   = BW_2 * 639   # ug Oral dose; 
  
  Opt_2    <- ev(ID=1, amt= DOSEoral2, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_2  <- Opt_2 
  
  #3
  # Exposure scenario #Inhalation_dermal_Opt_2
  BW_3          = 76
  DOSEoral3   = BW_3 * 781   # ug Oral dose; 
  
  Opt_3    <- ev(ID=1, amt= DOSEoral3, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_3  <- Opt_3 
  
  #4
  # Exposure scenario #Inhalation_dermal_Opt_2
  BW_4          = 74
  DOSEoral4   = BW_4 * 783   # ug Oral dose;  
  
  Opt_4    <- ev(ID=1, amt= DOSEoral4, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_4  <- Opt_4 
  
  #5
  # Exposure scenario #Inhalation_dermal_Opt_2
  BW_5          = 90
  DOSEoral5   = BW_5 * 775   # ug Oral dose;  
  
  Opt_5    <- ev(ID=1, amt= DOSEoral5, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_5  <- Opt_5 
  
  #5
  # Exposure scenario #Inhalation_dermal_Opt_2
  BW_6          = 108
  DOSEoral6   = BW_6 * 733   # ug Oral dose;  
  
  Opt_6    <- ev(ID=1, amt= DOSEoral6, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_6  <- Opt_6 
  
  ############
  
  
  # set up the simulation exposure time
  tsamp.B=tgrid(0,48,0.01)  # 48 hours
  
  ########### inhal+dermal 
  out.a <- 
    mod %>% 
    param(pars.data, BW = 83) %>%
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_1, tgrid=tsamp.B)
  
  out.a <- cbind.data.frame (Time =out.a$time,
                             Plasma_DPHP   = out.a$Plasma_DPHP,
                             Plasma_MPHP   = out.a$Plasma_MPHP,
                             Plasma_OH     = out.a$Plasma_OH,
                             Plasma_oxo    = out.a$Plasma_oxo,
                             Urine_MPHP    = (out.a$Urine_MPHP)/DOSEoral1*100,
                             Urine_OH      = (out.a$Urine_OH)/DOSEoral1*100,
                             Urine_oxo     = (out.a$Urine_oxo)/DOSEoral1*100)
  
  out.b <- 
    mod %>% 
    param(pars.data, BW = 75) %>%
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_2, tgrid=tsamp.B)
  
  out.b <- cbind.data.frame (Time =out.b$time,
                             Plasma_DPHP   = out.b$Plasma_DPHP,
                             Plasma_MPHP   = out.b$Plasma_MPHP,
                             Plasma_OH     = out.b$Plasma_OH,
                             Plasma_oxo    = out.b$Plasma_oxo,
                             Urine_MPHP    = (out.b$Urine_MPHP)/DOSEoral2*100,
                             Urine_OH      = (out.b$Urine_OH)/DOSEoral2*100,
                             Urine_oxo     = (out.b$Urine_oxo)/DOSEoral2*100)
  
  out.c <- 
    mod %>% 
    param(pars.data, BW = 76) %>%
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH,  Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_3, tgrid=tsamp.B)
  
  out.c <- cbind.data.frame (Time =out.c$time,
                             Plasma_DPHP   = out.c$Plasma_DPHP,
                             Plasma_MPHP   = out.c$Plasma_MPHP,
                             Plasma_OH     = out.c$Plasma_OH,
                             Plasma_oxo    = out.c$Plasma_oxo,
                             Urine_MPHP    = (out.c$Urine_MPHP)/DOSEoral3*100,
                             Urine_OH      = (out.c$Urine_OH)/DOSEoral3*100,
                             Urine_oxo     = (out.c$Urine_oxo)/DOSEoral3*100)
  
  
  out.d <- 
    mod %>% 
    param(pars.data, BW = 74) %>%
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_4, tgrid=tsamp.B)
  
  out.d <- cbind.data.frame (Time =out.d$time,
                             Plasma_DPHP   = out.d$Plasma_DPHP,
                             Plasma_MPHP   = out.d$Plasma_MPHP,
                             Plasma_OH     = out.d$Plasma_OH,
                             Plasma_oxo    = out.d$Plasma_oxo,
                             Urine_MPHP    = (out.d$Urine_MPHP)/DOSEoral4*100,
                             Urine_OH      = (out.d$Urine_OH)/DOSEoral4*100,
                             Urine_oxo     = (out.d$Urine_oxo)/DOSEoral4*100)
  
  out.e <- 
    mod %>% 
    param(pars.data, BW = 90) %>%
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_5, tgrid=tsamp.B)
  
  out.e <- cbind.data.frame (Time =out.e$time,
                             Plasma_DPHP   = out.e$Plasma_DPHP,
                             Plasma_MPHP   = out.e$Plasma_MPHP,
                             Plasma_OH     = out.e$Plasma_OH,
                             Plasma_oxo    = out.e$Plasma_oxo,
                             Urine_MPHP    = (out.e$Urine_MPHP)/DOSEoral5*100,
                             Urine_OH      = (out.e$Urine_OH)/DOSEoral5*100,
                             Urine_oxo     = (out.e$Urine_oxo)/DOSEoral5*100)
  
  out.f <- 
    mod %>% 
    param(pars.data, BW = 108) %>%
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_6, tgrid=tsamp.B)
  
  out.f <- cbind.data.frame (Time   = out.f$time,
                             Plasma_DPHP   = out.f$Plasma_DPHP,
                             Plasma_MPHP   = out.f$Plasma_MPHP,
                             Plasma_OH     = out.f$Plasma_OH,
                             Plasma_oxo    = out.f$Plasma_oxo,
                             Urine_MPHP    = (out.f$Urine_MPHP)/DOSEoral6*100,
                             Urine_OH      = (out.f$Urine_OH)/DOSEoral6*100,
                             Urine_oxo     = (out.f$Urine_oxo)/DOSEoral6*100)
  
  
  
  if (pred) return (list("out.a" = out.a,     ## Exposure scenario a
                         "out.b" = out.b,     ## Exposure scenario b 
                         "out.c" = out.c,     ## Exposure scenario c 
                         "out.d" = out.d,
                         "out.e" = out.e,
                         "out.f" = out.f))     ## Exposure scenario ee 
  
  ## Data mutate with prediction
  out.ap  = out.a  [which(out.a$Time  %in% Opt1_plasma_DPHP$Time),]
  out.bp  = out.b  [which(out.b$Time  %in% Opt2_plasma_DPHP$Time),]
  out.cp  = out.c  [which(out.c$Time  %in% Opt3_plasma_DPHP$Time),]
  out.dp  = out.d  [which(out.d$Time  %in% Opt4_plasma_DPHP$Time),]
  out.ep  = out.e  [which(out.e$Time  %in% Opt5_plasma_DPHP$Time),]
  out.fp  = out.f  [which(out.f$Time  %in% Opt6_plasma_DPHP$Time),]
  
  
  out.au  = out.a  [which(out.a$Time  %in% Opt1a$time),]
  out.bu  = out.b  [which(out.b$Time  %in% Opt2a$time),]
  out.cu  = out.c  [which(out.c$Time  %in% Opt3a$time),]
  out.du  = out.d  [which(out.d$Time  %in% Opt4a$time),]
  out.eu  = out.e  [which(out.e$Time  %in% Opt5a$time),]
  out.fu  = out.f  [which(out.f$Time  %in% Opt6a$time),]
  
  
  ## log.Predition 
  log.yhat.a.Plasma.DPHP               <-log(out.ap$Plasma_DPHP)
  log.yhat.a.Plasma.MPHP               <-log(out.ap$Plasma_MPHP)
  log.yhat.a.Plasma.OH                 <-log(out.ap$Plasma_OH)
  log.yhat.a.Plasma.oxo                <-log(out.ap$Plasma_oxo)
  log.yhat.a.Urine.MPHP                <-log(out.au$Urine_MPHP)
  log.yhat.a.Urine.OH                  <-log(out.au$Urine_OH)
  log.yhat.a.Urine.oxo                 <-log(out.au$Urine_oxo)
  
  log.yhat.b.Plasma.DPHP               <-log(out.bp$Plasma_DPHP)
  log.yhat.b.Plasma.MPHP               <-log(out.bp$Plasma_MPHP)
  log.yhat.b.Plasma.OH                 <-log(out.bp$Plasma_OH)
  log.yhat.b.Plasma.oxo                <-log(out.bp$Plasma_oxo)
  log.yhat.b.Urine.MPHP                <-log(out.bu$Urine_MPHP)
  log.yhat.b.Urine.OH                  <-log(out.bu$Urine_OH)
  log.yhat.b.Urine.oxo                 <-log(out.bu$Urine_oxo)
  
  log.yhat.c.Plasma.DPHP               <-log(out.cp$Plasma_DPHP)
  log.yhat.c.Plasma.MPHP               <-log(out.cp$Plasma_MPHP)
  log.yhat.c.Plasma.OH                 <-log(out.cp$Plasma_OH)
  log.yhat.c.Plasma.oxo                <-log(out.cp$Plasma_oxo)
  log.yhat.c.Urine.MPHP                <-log(out.cu$Urine_MPHP)
  log.yhat.c.Urine.OH                  <-log(out.cu$Urine_OH)
  log.yhat.c.Urine.oxo                 <-log(out.cu$Urine_oxo)
  
  log.yhat.d.Plasma.DPHP               <-log(out.dp$Plasma_DPHP)
  log.yhat.d.Plasma.MPHP               <-log(out.dp$Plasma_MPHP)
  log.yhat.d.Plasma.OH                 <-log(out.dp$Plasma_OH)
  log.yhat.d.Plasma.oxo                <-log(out.dp$Plasma_oxo)
  log.yhat.d.Urine.MPHP                <-log(out.du$Urine_MPHP)
  log.yhat.d.Urine.OH                  <-log(out.du$Urine_OH)
  log.yhat.d.Urine.oxo                 <-log(out.du$Urine_oxo)
  
  log.yhat.e.Plasma.DPHP               <-log(out.ep$Plasma_DPHP)
  log.yhat.e.Plasma.MPHP               <-log(out.ep$Plasma_MPHP)
  log.yhat.e.Plasma.OH                 <-log(out.ep$Plasma_OH)
  log.yhat.e.Plasma.oxo                <-log(out.ep$Plasma_oxo)
  log.yhat.e.Urine.MPHP                <-log(out.eu$Urine_MPHP)
  log.yhat.e.Urine.OH                  <-log(out.eu$Urine_OH)
  log.yhat.e.Urine.oxo                 <-log(out.eu$Urine_oxo)
  
  log.yhat.f.Plasma.DPHP               <-log(out.fp$Plasma_DPHP)
  log.yhat.f.Plasma.MPHP               <-log(out.fp$Plasma_MPHP)
  log.yhat.f.Plasma.OH                 <-log(out.fp$Plasma_OH)
  log.yhat.f.Plasma.oxo                <-log(out.fp$Plasma_oxo)
  log.yhat.f.Urine.MPHP                <-log(out.fu$Urine_MPHP)
  log.yhat.f.Urine.OH                  <-log(out.fu$Urine_OH)
  log.yhat.f.Urine.oxo                 <-log(out.fu$Urine_oxo)
  
  ## log.Observed data
  log.y.a.Plasma.DPHP               <-log(Opt1_plasma_DPHP$Plasma_DPHP)
  log.y.a.Plasma.MPHP               <-log(Opt1_plasma_MPHP$Plasma_MPHP)
  log.y.a.Plasma.OH                 <-log(Opt1_plasma_OH$Plasma_OH)
  log.y.a.Plasma.oxo                <-log(Opt1_plasma_oxo$Plasma_oxo)
  log.y.a.Urine.MPHP                <-log(Opt1_urine_MPHP$Urine_MPHP)
  log.y.a.Urine.OH                  <-log(Opt1_urine_OH$Urine_OH)
  log.y.a.Urine.oxo                 <-log(Opt1_urine_oxo$Urine_oxo)
  
  log.y.b.Plasma.DPHP               <-log(Opt2_plasma_DPHP$Plasma_DPHP)
  log.y.b.Plasma.MPHP               <-log(Opt2_plasma_MPHP$Plasma_MPHP)
  log.y.b.Plasma.OH                 <-log(Opt2_plasma_OH$Plasma_OH)
  log.y.b.Plasma.oxo                <-log(Opt2_plasma_oxo$Plasma_oxo)
  log.y.b.Urine.MPHP                <-log(Opt2_urine_MPHP$Urine_MPHP)
  log.y.b.Urine.OH                  <-log(Opt2_urine_OH$Urine_OH)
  log.y.b.Urine.oxo                 <-log(Opt2_urine_oxo$Urine_oxo)
  
  log.y.c.Plasma.DPHP               <-log(Opt3_plasma_DPHP$Plasma_DPHP)
  log.y.c.Plasma.MPHP               <-log(Opt3_plasma_MPHP$Plasma_MPHP)
  log.y.c.Plasma.OH                 <-log(Opt3_plasma_OH$Plasma_OH)
  log.y.c.Plasma.oxo                <-log(Opt3_plasma_oxo$Plasma_oxo)
  log.y.c.Urine.MPHP                <-log(Opt3_urine_MPHP$Urine_MPHP)
  log.y.c.Urine.OH                  <-log(Opt3_urine_OH$Urine_OH)
  log.y.c.Urine.oxo                 <-log(Opt3_urine_oxo$Urine_oxo)
  
  log.y.d.Plasma.DPHP               <-log(Opt4_plasma_DPHP$Plasma_DPHP)
  log.y.d.Plasma.MPHP               <-log(Opt4_plasma_MPHP$Plasma_MPHP)
  log.y.d.Plasma.OH                 <-log(Opt4_plasma_OH$Plasma_OH)
  log.y.d.Plasma.oxo                <-log(Opt4_plasma_oxo$Plasma_oxo)
  log.y.d.Urine.MPHP                <-log(Opt4_urine_MPHP$Urine_MPHP)
  log.y.d.Urine.OH                  <-log(Opt4_urine_OH$Urine_OH)
  log.y.d.Urine.oxo                 <-log(Opt4_urine_oxo$Urine_oxo)
  
  log.y.e.Plasma.DPHP               <-log(Opt5_plasma_DPHP$Plasma_DPHP)
  log.y.e.Plasma.MPHP               <-log(Opt5_plasma_MPHP$Plasma_MPHP)
  log.y.e.Plasma.OH                 <-log(Opt5_plasma_OH$Plasma_OH)
  log.y.e.Plasma.oxo                <-log(Opt5_plasma_oxo$Plasma_oxo)
  log.y.e.Urine.MPHP                <-log(Opt5_urine_MPHP$Urine_MPHP)
  log.y.e.Urine.OH                  <-log(Opt5_urine_OH$Urine_OH)
  log.y.e.Urine.oxo                 <-log(Opt5_urine_oxo$Urine_oxo)
  
  log.y.f.Plasma.DPHP               <-log(Opt6_plasma_DPHP$Plasma_DPHP)
  log.y.f.Plasma.MPHP               <-log(Opt6_plasma_MPHP$Plasma_MPHP)
  log.y.f.Plasma.OH                 <-log(Opt6_plasma_OH$Plasma_OH)
  log.y.f.Plasma.oxo                <-log(Opt6_plasma_oxo$Plasma_oxo)
  log.y.f.Urine.MPHP                <-log(Opt6_urine_MPHP$Urine_MPHP)
  log.y.f.Urine.OH                  <-log(Opt6_urine_OH$Urine_OH)
  log.y.f.Urine.oxo                 <-log(Opt6_urine_oxo$Urine_oxo)
  
  
  log.yhat        <- c(log.yhat.a.Plasma.DPHP, log.yhat.a.Plasma.MPHP, log.yhat.a.Plasma.OH,log.yhat.a.Plasma.oxo,
                       log.yhat.a.Urine.MPHP,  log.yhat.a.Urine.OH,    log.yhat.a.Urine.oxo, 
                       
                       log.yhat.b.Plasma.DPHP, log.yhat.b.Plasma.MPHP, log.yhat.b.Plasma.OH,log.yhat.b.Plasma.oxo,
                       log.yhat.b.Urine.MPHP,  log.yhat.b.Urine.OH,    log.yhat.b.Urine.oxo,  
                       
                       log.yhat.c.Plasma.DPHP, log.yhat.c.Plasma.MPHP, log.yhat.c.Plasma.OH,log.yhat.c.Plasma.oxo,
                       log.yhat.c.Urine.MPHP,  log.yhat.c.Urine.OH,    log.yhat.c.Urine.oxo,  
                       
                       log.yhat.d.Plasma.DPHP, log.yhat.d.Plasma.MPHP, log.yhat.d.Plasma.OH,log.yhat.d.Plasma.oxo,
                       log.yhat.d.Urine.MPHP,  log.yhat.d.Urine.OH,    log.yhat.d.Urine.oxo,
                       
                       log.yhat.e.Plasma.DPHP, log.yhat.e.Plasma.MPHP, log.yhat.e.Plasma.OH,log.yhat.e.Plasma.oxo,
                       log.yhat.e.Urine.MPHP,  log.yhat.e.Urine.OH,    log.yhat.e.Urine.oxo, 
                       
                       log.yhat.f.Plasma.DPHP, log.yhat.f.Plasma.MPHP, log.yhat.f.Plasma.OH,log.yhat.f.Plasma.oxo,
                       log.yhat.f.Urine.MPHP,  log.yhat.f.Urine.OH,    log.yhat.f.Urine.oxo)   
  
  log.y           <- c(log.y.a.Plasma.DPHP, log.y.a.Plasma.MPHP, log.y.a.Plasma.OH,log.y.a.Plasma.oxo,
                       log.y.a.Urine.MPHP,  log.y.a.Urine.OH,    log.y.a.Urine.oxo, 
                       
                       log.y.b.Plasma.DPHP, log.y.b.Plasma.MPHP, log.y.b.Plasma.OH,log.y.b.Plasma.oxo,
                       log.y.b.Urine.MPHP,  log.y.b.Urine.OH,    log.y.b.Urine.oxo,
                       
                       log.y.c.Plasma.DPHP, log.y.c.Plasma.MPHP, log.y.c.Plasma.OH,log.y.c.Plasma.oxo,
                       log.y.c.Urine.MPHP,  log.y.c.Urine.OH,    log.y.c.Urine.oxo,  
                       
                       log.y.d.Plasma.DPHP, log.y.d.Plasma.MPHP, log.y.d.Plasma.OH,log.y.d.Plasma.oxo,
                       log.y.d.Urine.MPHP,  log.y.d.Urine.OH,    log.y.d.Urine.oxo,
                       
                       log.y.e.Plasma.DPHP, log.y.e.Plasma.MPHP, log.y.e.Plasma.OH,log.y.e.Plasma.oxo,
                       log.y.e.Urine.MPHP,  log.y.e.Urine.OH,    log.y.e.Urine.oxo,
                       
                       log.y.f.Plasma.DPHP, log.y.f.Plasma.MPHP, log.y.f.Plasma.OH,log.y.f.Plasma.oxo,
                       log.y.f.Urine.MPHP,  log.y.f.Urine.OH,    log.y.f.Urine.oxo)   
  
  #log.yhat        <- c(log.yhat.a.DPHP, log.yhat.a.MPHP,     
  #    log.yhat.b.DPHP, log.yhat.b.MPHP,       
  #  log.yhat.c.DPHP, log.yhat.c.MPHP,       
  #  log.yhat.d.DPHP, log.yhat.d.MPHP, 
  #  log.yhat.e.DPHP, log.yhat.e.MPHP,      
  #  log.yhat.f.DPHP, log.yhat.f.MPHP)   
  
  #log.y           <- c(log.y.a.DPHP, log.y.a.MPHP,     
  # log.y.b.DPHP, log.y.b.MPHP,        
  # log.y.c.DPHP, log.y.c.MPHP,      
  #log.y.d.DPHP, log.y.d.MPHP, 
  # log.y.e.DPHP, log.y.e.MPHP,      
  #log.y.f.DPHP, log.y.f.MPHP)   
  
  # The method of Maximum likelihood
  sig2            <- as.numeric((exp(pars[which_sig][1])))
  log_likelihood  <- -2*sum ((dnorm (log.y,
                                     mean =log.yhat,
                                     sd = sqrt(sig2),
                                     log = TRUE)))
  return(log_likelihood)
  
}

## Define the Prior distributions: either normal or lognormal distribution
## nomral distribution
Prior <- function(pars) {
  
  ## Population level
  # The likelihood for population mean (parameters)
  pars.data = exp(pars[-which_sig])
  sig  <- as.numeric (exp(pars[which_sig][2:15]))                 # Coefficient of variation from population variance; sigmal0
  sig2 <- as.numeric (exp(pars[which_sig][1]))                    # error variances from model esidual
  
  mean           = exp(theta.MCMC[-which_sig])
  CV             = 1                                            # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
  sd             = mean*CV
  
  # Calculate likelihoods of each parameters; P(u|M,S)
  prior_pars     = dtruncnorm(pars.data, 
                              a = qnorm(0.025, mean = mean, sd = sd), 
                              b = qnorm(0.975, mean = mean, sd = sd), 
                              mean = mean, sd = sd ) 
  
  # The likelihood for population variance; P(sigmal^2|sigmal0^2)
  CU             = 1                                              # Coefficient of uncertainty (CU) (Hack et al., 2006)
  CV.sig         = exp(theta.MCMC[which_sig])[2:15]               # Singmal0
  alpha          = (2+1)/(CU^2)                                   # Shape parametrer of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  beta           = (alpha-1)*CV.sig^2                             # Scale parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  
  # Calculate likelihoods of model error (sig2) and population variance (sig) parameters
  prior_sig      = dinvgamma (sig, shape = alpha , scale = beta)  # Prior distribution for population vraicne; sigma2
  prior_sig2     = dunif (sig2, min = 0.01, max = 3.3)            # Error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014)   
  
  ## individual level; P(theta|u,sigmal^2)
  mean_i         = prior_pars
  sd_i           = sqrt(prior_sig)
  prior_pars_i   = dtruncnorm (prior_pars, 
                               a = qnorm(0.025, mean = mean_i, sd = sd_i), 
                               b = qnorm(0.975, mean = mean_i, sd = sd_i), 
                               mean = mean_i, sd = sd_i) 
  
  # log-transformed (log-likelihoods of each parameters)
  log.pri.pars   = log (prior_pars)
  log.pri.sig    = log (prior_sig)
  log.pri.pars.i = log (prior_pars_i)
  log.pri.sig2   = log (prior_sig2)
  
  # Maximau likelihood estimation (MLE): negative log-likelihood function, (-2 times sum of log-likelihoods)
  MLE =  -2*sum(log.pri.pars, log.pri.sig , log.pri.pars.i,log.pri.sig2)  
  
  return(MLE)
}



#################### MCMC simulation with parallel computing ############################


#install.packages(c('doMC', 'foreach')) 
library(foreach) 
library(doMC)
#install.packages("doSNOW")
#install.packages("doParallel") 
library("doParallel")
#install.packages("doMPI")
library(parallel)
#install.packages("drake") 
library(drake)
mrgsolve::loadso
loadso(mod)

detectCores()                                ## check the cores
cl<- makeCluster(detectCores())              ## use all cores in our system     
registerDoParallel(cl)                       ## registers a cluster of all the cores on our system

# start time
strt<-Sys.time()

# parallel
system.time(
  MCMC <- foreach( i = 1:4, .packages = c('mrgsolve','magrittr','FME','truncnorm','EnvStats','invgamma','dplyr')) %dopar% {
    mod <- mcode ("pk", HumanPK.code)
    loadso(mod)
    modMCMC(f             = mcmc.fun, 
            p             = theta.MCMC, 
            niter         = 300000,           ## iteration number 
            jump          = 0.01,             ## jump function generation new parameters distribution using covariate matrix
            lower = c(-2,  -5,-0.5,  -0.9,-1,    -8,  -0.5,-1.9,  -5, -2.9,    -5,  -5,   -5,-5, -Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf),
            upper = c(1,  0.1,   2,   0.1, 1.5,  -5,    0.5, 2,   -1, -0.5,  0.01,-0.1,  0.1, 1, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
            prior         = Prior,            ## prior function
            updatecov     = 50,               ## Adaptative Metropolis
            var0          = NULL,             ## initial model variance;
            wvar0         = 0.01,             ## "Weight" for the initial model variance
            ntrydr        = 2,                ## Delayed Rejection
            burninlength  = 150000,           ## number of initial iterations to be removed from output.
            outputlength  = 5000)            ## number of output iterations           
    
  }
)


#end time
print(Sys.time()-strt)

stopCluster(cl)   


## Performance four chains to check the convergences
MC.H.1 = as.mcmc (MCMC[[1]]$pars)     # first  chain
MC.H.2 = as.mcmc (MCMC[[2]]$pars)     # second chain
MC.H.3 = as.mcmc (MCMC[[3]]$pars)     # third  chain
MC.H.4 = as.mcmc (MCMC[[4]]$pars)     # fourth chain
combinedchains = mcmc.list(MC.H.1,MC.H.2,MC.H.3,MC.H.4) ## combine all chains
gelman.diag (combinedchains)          # Gelman convergence diagnosis
heidel.diag (combinedchains)          # covergence diagnosis/Heidelberger and Welch's convergence diagnostic
gelman.plot (combinedchains)          # gelman plot

# Save the poseior parameters (95% CI)
quan.Human = exp(summary(MC.H.1)$quantiles)  


# Trace plot using bayesplot
## Covergences plot
#color_scheme_set("gray")
mcmc_trace (
  combinedchains,
  pars =names(theta.MCMC),
  size = 0.5,
  facet_args = list(nrow = 2)) +
  ggplot2::scale_color_brewer()


## output the MCMC results
write.csv(quan.Human,file="~/PK/DPHP/result/Human.summary_pos.csv")
write.csv(MC.H.1,file="~/PK/DPHP/result/Human.pos.csv")
saveRDS(MCMC[[1]],file ='~/PK/DPHP/result/Human.MCMC.rds')
saveRDS(combinedchains,file='~/PK/DPHP/result/Human.comb.rds')


## Plot using MCMC parameters

MCMC.B = mcmc.fun (par = MCMC[[1]]$bestpar, pred=TRUE)$out.B
MCMC.C = mcmc.fun (par = MCMC[[1]]$bestpar, pred=TRUE)$out.C
df.B=cbind.data.frame (Time  = MCMC.B$Time, 
                       CA2    = MCMC.B$CA2)
df.C=cbind.data.frame (Time  = MCMC.C$Time, 
                       AU3    = MCMC.C$AU3)

Plot_MCMC =
  ggplot() +
  geom_line(data  = df.B,aes(Time,CA2), col="firebrick", lwd=2)+
  geom_line(data  = df.C,aes(Time,AU3), col="green", lwd=2) +
  geom_point(data = Human.2,aes(Time,CA2),col="firebrick",size=2) + ylab("Conc") +
  geom_point(data = Human.3,aes(Time,AU3),col="green",size=2) 
Plot_MCMC 


plot.B1=
  ggplot() +
  geom_line(data  = df.B,aes(Time,CA2), col="firebrick", lwd=2)+
  geom_point(data = Human.2,aes(Time,CA2),col="firebrick",size=2) + ylab("Conc")

plot.B1

ggsave("Plot_MCMC.tiff",scale = 1,
       plot = Plot_MCMC ,
       path = "~/PK/DPHP/",
       width = 25, height = 20, units = "cm",dpi=320)

