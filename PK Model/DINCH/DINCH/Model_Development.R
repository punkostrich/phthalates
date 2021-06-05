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

######## DINCH ################
############################################# Model Calibration with levenberg-marquart #############

#Plasma.DEHP <- read.csv(file="~/Dropbox/R/DEHP/2012_Plasma_DEHP.csv")

Cal1 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DINCH/Cal1.csv")
Cal1_MHNCH <- cbind.data.frame (Time = Cal1$time,
                                MHNCH   = Cal1$OH.MINCH)
Cal1_MINCH <- cbind.data.frame (Time = Cal1$time,
                                MINCH   = Cal1$MINCH)
Cal1_oxo <- cbind.data.frame (Time = Cal1$time,
                                oxo   = Cal1$oxo.MINCH)
Cal1_cx <- cbind.data.frame (Time = Cal1$time,
                                cx   = Cal1$cx.MINCH)
Cal1_CHDA <- cbind.data.frame (Time = Cal1$time,
                                CHDA   = Cal1$CHDA)

Cal2 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DINCH/Cal2.csv")
Cal2_MHNCH <- cbind.data.frame (Time = Cal2$time,
                                MHNCH   = Cal2$OH.MINCH)
Cal2_MINCH <- cbind.data.frame (Time = Cal2$time,
                                MINCH   = Cal2$MINCH)
Cal2_oxo <- cbind.data.frame (Time = Cal2$time,
                              oxo   = Cal2$oxo.MINCH)
Cal2_cx <- cbind.data.frame (Time = Cal2$time,
                             cx   = Cal2$cx.MINCH)
Cal2_CHDA <- cbind.data.frame (Time = Cal2$time,
                               CHDA   = Cal2$CHDA)




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
  BW_1          = 66.2
  DOSEoral1   = BW_1 * 1000   # ug Oral dose; 
  
  Cal_1    <- ev(ID=1, amt= DOSEoral1, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_cal_1  <- Cal_1 
  
  #2
  # Exposure scenario #Inhalation_dermal_cal_2
  BW_2          = 57
  DOSEoral2   = BW_2 * 1000   # ug Oral dose; 
  
  Cal_2    <- ev(ID=1, amt= DOSEoral2, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_cal_2  <- Cal_2 
  
 
  # set up the exposure time
  tsamp=tgrid(0,72,0.1)  # 48 hours
  
  out.1 <- 
    
    mod %>% 
    param(pars, BW = 66.2)%>% 
    Req(Urine_MHNCH, Urine_MINCH, Urine_cx, Urine_oxo, Urine_CHDA) %>%
    update(atol = 1E-15,maxsteps = 100000000) %>%
    mrgsim_d(data = ex_cal_1, tgrid=tsamp)
  out.1<-cbind.data.frame (Time = out.1$time,
                           MHNCH   = (out.1$Urine_MHNCH)/BW_1,
                           MINCH   = (out.1$Urine_MINCH)/BW_1,
                           cx   = (out.1$Urine_cx)/BW_1,
                           oxo   = (out.1$Urine_oxo)/BW_1,
                           CHDA   = (out.1$Urine_CHDA)/BW_1)
  
  
  out.2 <-
    
    mod %>% 
    param(pars,BW = 57) %>%  
    Req(Urine_MHNCH, Urine_MINCH, Urine_cx, Urine_oxo, Urine_CHDA)%>%
    update(atol = 1E-15,maxsteps = 1000000000) %>%
    mrgsim_d(data = ex_cal_2, tgrid=tsamp)
  out.2<-cbind.data.frame (Time = out.2$time,
                           MHNCH   = (out.2$Urine_MHNCH)/BW_2,
                           MINCH   = (out.2$Urine_MINCH)/BW_2,
                           cx   = (out.2$Urine_cx)/BW_2,
                           oxo   = (out.2$Urine_oxo)/BW_2,
                           CHDA   = (out.2$Urine_CHDA)/BW_2)
  
  
  return(list("out.1" = out.1,
              "out.2" = out.2))
  
}

## Cost fuction (from FME pckage) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out <-  Pred (pars)
  cost<- modCost(model=out$out.1, obs= Cal1_MHNCH, weight='std', x="Time")
  cost<- modCost(model=out$out.1, obs= Cal1_MINCH, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_cx, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_oxo, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_CHDA, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_MHNCH, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_MINCH, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_cx, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_oxo, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Cal2_CHDA, weight='std', x="Time", cost = cost)
  return(cost)
}



###################### Local sensitivity analysis #########################
## Choose the senstiive parameters in the model                          ##
## initial parmaeters######################################################

theta.int <- log(c(
  Vp                  =   5.3,     # Volume of distribution parent compound L/kg bw
  Vm                  =   0.05,    #: Volume of distribution metabolite compound 
  
  ka                  =   0.2,     #: stomach absorption rate parent 1/h
  
  BW                  =   82.3,    #: kg, Bodyweight (EPA Factors Handbook, 2011)
  
  lambda_p            =   0.08,    #: lambda for parent (h-1)
  
  
  lambda_zm      =   0.2,     #: h-1
  lambda_z_MHNCH      =   1,       #: h-1
  lambda_z_cx         =   1,       #: h-1
  lambda_z_oxo        =   1,       #: h-1
  lambda_z_CHDA       =   1,       #: h-1
  
  lambda_u       =   0.1,     #: h-1
  lambda_u_MHNCH      =   0.5,     #: h-1
  lambda_u_cx         =   0.2,     #: h-1
  lambda_u_oxo        =   0.2,     #: h-1
  lambda_u_CHDA       =   1,       #: h-1
  
  Fr_MINCH            =   0.5,     #: h-1
  Fr_MHNCH            =   0.5,     #: h-1
  Fr_cx               =   0.17,    #: h-1
  Fr_oxo              =   0.13    #: h-1
  
))

# Senstivity function (FME)
# Check the senstive parameters in the modelhttp://ec2-3-135-223-180.us-east-2.compute.amazonaws.com/graphics/plot_zoom_png?width=1200&height=868
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)
Sen=summary(SnsPlasma)
write.csv(Sen,file="C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/Sen.human.csv")

plot(summary(SnsPlasma))





######################### MCMC optimization ##########################
## set up senstivie or necessary parametesr as model input          ##
######################################################################
theta <- log(c(
  #Vp                  =   5.3,     # Volume of distribution parent compound L/kg bw
  #Vm                  =   0.05,    #: Volume of distribution metabolite compound 
  
  ka                  =   0.2,     #: stomach absorption rate parent 1/h
  
  #BW                  =   82.3,    #: kg, Bodyweight (EPA Factors Handbook, 2011)
  
  lambda_p            =   0.08,    #: lambda for parent (h-1)
  
  
  lambda_zm      =   0.2,     #: h-1
  #lambda_z_MHNCH      =   1,       #: h-1
  #lambda_z_cx         =   1,       #: h-1
  #lambda_z_oxo        =   1,       #: h-1
  #lambda_z_CHDA       =   1,       #: h-1
  
  lambda_u       =   0.1,     #: h-1
  lambda_u_MHNCH      =   0.5,     #: h-1
  lambda_u_cx         =   0.2,     #: h-1
  lambda_u_oxo        =   0.2,     #: h-1
  #lambda_u_CHDA       =   1,       #: h-1
  #Fr_MINCH            =   0.5,     #: h-1
  Fr_MHNCH            =   0.5,     #: h-1
  Fr_cx               =   0.17,    #: h-1
  Fr_oxo              =   0.13    #: h-1
))

# Model fitting
#upper = c(Inf,Inf,Inf,Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
system.time(Fit<- modFit(f=MCcost, p=theta, method = "Nelder-Mead",lower = c(-2,-5,-5,-5,-5,-5,-5,-5,-5,-5),
                         upper = c(5,3,Inf,Inf,Inf,Inf,Inf,0.01,0.01,0.01),
                         control = nls.lm.control(nprint=1)))

summary(Fit)
### para value
exp(Fit$par)

res=MCcost(Fit$par)$residuals$res
sum(res^2)


## Model calibration plot using ggplot2 
Sim.fit.A = Pred (Fit$par)$out.1         ## Time-course concentration profiles using estimated parameters under exposure senaior 1
Sim.fit.B = Pred (Fit$par)$out.2         ## Simulaiton of exposure scenaior 2


df.sim.A1  = cbind.data.frame (Time=Sim.fit.A$Time, MHNCH=Sim.fit.A$MHNCH)
df.sim.A2  = cbind.data.frame (Time=Sim.fit.A$Time, MINCH=Sim.fit.A$MINCH)
df.sim.A3  = cbind.data.frame (Time=Sim.fit.A$Time, cx=Sim.fit.A$cx)
df.sim.A4  = cbind.data.frame (Time=Sim.fit.A$Time, oxo=Sim.fit.A$oxo)
df.sim.A5  = cbind.data.frame (Time=Sim.fit.A$Time, CHDA=Sim.fit.A$CHDA)
df.sim.B1  = cbind.data.frame (Time=Sim.fit.B$Time, MHNCH=Sim.fit.B$MHNCH)
df.sim.B2  = cbind.data.frame (Time=Sim.fit.B$Time, MINCH=Sim.fit.B$MINCH)
df.sim.B3  = cbind.data.frame (Time=Sim.fit.B$Time, cx=Sim.fit.B$cx)
df.sim.B4  = cbind.data.frame (Time=Sim.fit.B$Time, oxo=Sim.fit.B$oxo)
df.sim.B5  = cbind.data.frame (Time=Sim.fit.B$Time, CHDA=Sim.fit.B$CHDA)


## Plot
plot.A1=
  ggplot() +
  geom_line (data = df.sim.A1,aes(Time,MHNCH), col="firebrick", lwd=2)+
  geom_point(data = Cal1_MHNCH, aes(Time, MHNCH),size=2.5) + ylab("Concentration") 
plot.A1


plot.A2=
  ggplot() +
  geom_line (data = df.sim.A2,aes(Time,MINCH), col="firebrick", lwd=2)+
  geom_point(data = Cal1_MINCH, aes(Time, MINCH),size=2.5) + ylab("Concentration") 
plot.A2


plot.A3=
  ggplot() +
  geom_line (data = df.sim.A3,aes(Time,cx), col="firebrick", lwd=2)+
  geom_point(data = Cal1_cx,aes(Time, cx),size=2.5) + ylab("Concentration") 
plot.A3


plot.A4=
  ggplot() +
  geom_line (data = df.sim.A4,aes(Time,oxo), col="firebrick", lwd=2)+
  geom_point(data = Cal1_oxo,aes(Time, oxo),size=2.5) + ylab("Concentration") 
plot.A4


plot.A5=
  ggplot() +
  geom_line (data = df.sim.A5,aes(Time,CHDA), col="firebrick", lwd=2)+
  geom_point(data = Cal1_CHDA,aes(Time, CHDA),size=2.5) + ylab("Concentration") 
plot.A5


plot.B1=
  ggplot() +
  geom_line (data = df.sim.B1,aes(Time,MHNCH), col="firebrick", lwd=2)+
  geom_point(data = Cal2_MHNCH, aes(Time, MHNCH),size=2.5) + ylab("Concentration") 
plot.B1


plot.B2=
  ggplot() +
  geom_line (data = df.sim.B2,aes(Time,MINCH), col="firebrick", lwd=2)+
  geom_point(data = Cal2_MINCH, aes(Time, MINCH),size=2.5) + ylab("Concentration") 
plot.B2


plot.B3=
  ggplot() +
  geom_line (data = df.sim.B3,aes(Time,cx), col="firebrick", lwd=2)+
  geom_point(data = Cal2_cx,aes(Time, cx),size=2.5) + ylab("Concentration") 
plot.B3


plot.B4=
  ggplot() +
  geom_line (data = df.sim.B4,aes(Time,oxo), col="firebrick", lwd=2)+
  geom_point(data = Cal2_oxo,aes(Time, oxo),size=2.5) + ylab("Concentration") 
plot.B4


plot.B5=
  ggplot() +
  geom_line (data = df.sim.B5,aes(Time,CHDA), col="firebrick", lwd=2)+
  geom_point(data = Cal2_CHDA,aes(Time, CHDA),size=2.5) + ylab("Concentration") 
plot.B5




plot.A7=
  ggplot() +
  geom_line (data = df.sim.A1,aes(Time,MHNCH), col="firebrick", lwd=2)+
  geom_point(data = Cal1_MHNCH, aes(Time, MHNCH),col="firebrick",size=2.5) + ylab("Concentration") +
  geom_line (data = df.sim.A2,aes(Time,MINCH), col="green", lwd=2)+
  geom_point(data = Cal1_MINCH, aes(Time, MINCH),col="green",size=2.5) + 
  geom_line (data = df.sim.A3,aes(Time,cx), col="blue", lwd=2)+
  geom_point(data = Cal1_cx,aes(Time, cx),col="blue",size=2.5) +
  geom_line (data = df.sim.A4,aes(Time,oxo), col="black", lwd=2)+
  geom_point(data = Cal1_oxo,aes(Time, oxo),col="black",size=2.5) + 
  geom_line (data = df.sim.A5,aes(Time,CHDA), col="tan", lwd=2)+
  geom_point(data = Cal1_CHDA,aes(Time, CHDA),col="tan",size=2.5) + 
  geom_line (data = df.sim.B1,aes(Time,MHNCH), col="orange", lwd=2)+
  geom_point(data = Cal2_MHNCH, aes(Time, MHNCH),col="orange",size=2.5) +  
  geom_line (data = df.sim.B2,aes(Time,MINCH), col="red", lwd=2)+
  geom_point(data = Cal2_MINCH, aes(Time, MINCH),col="red",size=2.5) + 
  geom_line (data = df.sim.B3,aes(Time,cx), col="steelblue", lwd=2)+
  geom_point(data = Cal2_cx,aes(Time, cx),col="steelblue",size=2.5) +
  geom_line (data = df.sim.B4,aes(Time,oxo), col="green4", lwd=2)+
  geom_point(data = Cal2_oxo,aes(Time, oxo),col="green4",size=2.5) + 
  geom_line (data = df.sim.B5,aes(Time,CHDA), col="white", lwd=2)+
  geom_point(data = Cal2_CHDA,aes(Time, CHDA),col="white",size=2.5) 

plot.A7
plot.A1
plot.A2
plot.A3
plot.A4
plot.A5
plot.B1
plot.B2
plot.B3
plot.B4
plot.B5





############################################# Model Calibration with MCMC ###################################################
## EIGHT data sets was used in model evaluation                                                                               #                                                       #
#############################################################################################################################

# input callibrated data set;
#Human.2 <- read.csv(file="~/Dropbox/R/DEHP/2004_Plasma.csv")
Opt1 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DINCH/Opt1.csv")
Opt1 <- Opt1[Opt1$time != 0, ]
Opt1 <- Opt1[Opt1$OH.MINCH != 0, ]
Opt1 <- Opt1[Opt1$MINCH != 0, ]
Opt1 <- Opt1[Opt1$oxo.MINCH != 0, ]
Opt1 <- Opt1[Opt1$cx.MINCH != 0, ]
Opt1 <- Opt1[Opt1$CHDA != 0, ]
Opt1 <- na.omit(Opt1)
Opt1_MHNCH <- cbind.data.frame (Time = Opt1$time,
                                MHNCH   = Opt1$OH.MINCH)
Opt1_MINCH <- cbind.data.frame (Time = Opt1$time,
                                MINCH   = Opt1$MINCH)
Opt1_oxo <- cbind.data.frame (Time = Opt1$time,
                              oxo   = Opt1$oxo.MINCH)
Opt1_cx <- cbind.data.frame (Time = Opt1$time,
                             cx   = Opt1$cx.MINCH)


Opt1_CHDA <- cbind.data.frame (Time = Opt1$time,
                               CHDA   = Opt1$CHDA)


#Opt1_MHNCH <- na.omit(Opt1_MHNCH)
#Opt1_MINCH <- na.omit(Opt1_MINCH)
#Opt1_cx <- na.omit(Opt1_cx)
#Opt1_oxo <- na.omit(Opt1_oxo)
#Opt1_CHDA <- na.omit(Opt1_CHDA)


Opt2 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DINCH/Opt2.csv")
Opt2 <- Opt2[Opt2$time != 0, ]
Opt2 <- Opt2[Opt2$OH.MINCH != 0, ]
Opt2 <- Opt2[Opt2$MINCH != 0, ]
Opt2 <- Opt2[Opt2$oxo.MINCH != 0, ]
Opt2 <- Opt2[Opt2$cx.MINCH != 0, ]
Opt2 <- Opt2[Opt2$CHDA != 0, ]
Opt2 <- na.omit(Opt2)

Opt2_MHNCH <- cbind.data.frame (Time = Opt2$time,
                                MHNCH   = Opt2$OH.MINCH)
Opt2_MINCH <- cbind.data.frame (Time = Opt2$time,
                                MINCH   = Opt2$MINCH)
Opt2_oxo <- cbind.data.frame (Time = Opt2$time,
                              oxo   = Opt2$oxo.MINCH)
Opt2_cx <- cbind.data.frame (Time = Opt2$time,
                             cx   = Opt2$cx.MINCH)
Opt2_CHDA <- cbind.data.frame (Time = Opt2$time,
                               CHDA   = Opt2$CHDA)

#Opt2_MHNCH <- na.omit(Opt2_MHNCH)
#Opt2_MINCH <- na.omit(Opt2_MINCH)
#Opt2_cx <- na.omit(Opt2_cx)
#Opt2_oxo <- na.omit(Opt2_oxo)
#Opt2_CHDA <- na.omit(Opt2_CHDA)


Opt3 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DINCH/Opt3.csv")
Opt3 <- Opt3[Opt3$time != 0, ]
Opt3 <- Opt3[Opt3$OH.MINCH != 0, ]
Opt3 <- Opt3[Opt3$MINCH != 0, ]
Opt3 <- Opt3[Opt3$oxo.MINCH != 0, ]
Opt3 <- Opt3[Opt3$cx.MINCH != 0, ]
Opt3 <- na.omit(Opt3)

Opt3_MHNCH <- cbind.data.frame (Time = Opt3$time,
                                MHNCH   = Opt3$OH.MINCH)
Opt3_MINCH <- cbind.data.frame (Time = Opt3$time,
                                MINCH   = Opt3$MINCH)
Opt3_oxo <- cbind.data.frame (Time = Opt3$time,
                              oxo   = Opt3$oxo.MINCH)
Opt3_cx <- cbind.data.frame (Time = Opt3$time,
                             cx   = Opt3$cx.MINCH)


#Opt3_MHNCH <- Opt3_MHNCH[Opt3_MHNCH$MHNCH != 0, ] 
#Opt3_MINCH <- Opt3_MINCH[Opt3_MINCH$MINCH != 0,]
#Opt3_cx    <- Opt3_cx[Opt3_cx$cx != 0,]
#Opt3_oxo   <- Opt3_cx[Opt3_oxo$oxo != 0,]

#Opt3_MHNCH <- na.omit(Opt3_MHNCH)
#Opt3_MINCH <- na.omit(Opt3_MINCH)
#Opt3_cx <- na.omit(Opt3_cx)
#Opt3_oxo <- na.omit(Opt3_oxo)


Opt4 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DINCH/Opt4.csv")
Opt4 <- Opt4[Opt4$time != 0, ]
Opt4 <- Opt4[Opt4$OH.MINCH != 0, ]
Opt4 <- Opt4[Opt4$MINCH != 0, ]
Opt4 <- Opt4[Opt4$cx.MINCH != 0, ]
Opt4 <- na.omit(Opt4)

Opt4_MHNCH <- cbind.data.frame (Time = Opt4$time,
                                MHNCH   = Opt4$OH.MINCH)
Opt4_MINCH <- cbind.data.frame (Time = Opt4$time,
                                MINCH   = Opt4$MINCH)
Opt4_cx <- cbind.data.frame (Time = Opt4$time,
                             cx   = Opt4$cx.MINCH)



#Opt4_MHNCH <- Opt4_MHNCH[Opt4_MHNCH$MHNCH != 0, ] 
#Opt4_MINCH <- Opt4_MINCH[Opt4_MINCH$MINCH != 0,]
#Opt4_cx    <- Opt4_cx[Opt4_cx$cx != 0,]

#Opt4_MHNCH <- na.omit(Opt4_MHNCH)
#Opt4_MINCH <- na.omit(Opt4_MINCH)
#Opt4_cx <- na.omit(Opt4_cx)



## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)

theta.MCMC <- log(c(
  Vp                  =   5.3,     # Volume of distribution parent compound L/kg bw
  Vm                  =   0.5,    #: Volume of distribution metabolite compound 
  
  ka                  =   0.19655912,     #: stomach absorption rate parent 1/h
  
  lambda_p            =   0.33092431,    #: lambda for parent (h-1)
  lambda_zm           =   2.42063938,     #: h-1
  lambda_u            =   0.03410414,     #: h-1
  
  lambda_z_MHNCH      =   1,       #: h-1
  lambda_z_cx         =   1,       #: h-1
  lambda_z_oxo        =   1,       #: h-1
  lambda_z_CHDA       =   1,       #: h-1
  
  lambda_u_MHNCH      =   0.25295638,     #: h-1
  lambda_u_cx         =   0.11143978,     #: h-1
  lambda_u_oxo        =   0.69180619,     #: h-1
  lambda_u_CHDA       =   1,       #: h-1
  
  Fr_MINCH            =   0.5,     #: h-1
  Fr_MHNCH            =   0.45026139,     #: h-1
  Fr_cx               =   0.27608643,    #: h-1
  Fr_oxo              =   0.15636782,    #: h-1
  
  sig2                =   1,      ## Model residuals; mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of most parametesr is 0.3, some parameters is 0.5 due to possible variation)
  sig_Vp                = 1,                ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_Vm                = 1,
  sig_ka                = 1,

  sig_lambda_p          = 1,
  sig_lambda_zm         = 1,
  sig_lambda_u          = 1,
  
  sig_lambda_z_MHNCH      =   1,       #: h-1
  sig_lambda_z_cx         =   1,       #: h-1
  sig_lambda_z_oxo        =   1,       #: h-1
  sig_lambda_z_CHDA       =   1,       #: h-1
  
  sig_lambda_u_MHNCH      =   1,     #: h-1
  sig_lambda_u_cx         =   1,     #: h-1
  sig_lambda_u_oxo        =   1,     #: h-1
  sig_lambda_u_CHDA       =   1,       #: h-1
  
  sig_Fr_MINCH            =   1,     #: h-1
  sig_Fr_MHNCH            =   1,     #: h-1
  sig_Fr_cx               =   1,    #: h-1
  sig_Fr_oxo              =   1    #: h-1

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
  BW_1          = 73
  DOSEoral1   = BW_1 * 1000   # ug Oral dose; 
  
  Opt_1    <- ev(ID=1, amt= DOSEoral1, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_1  <- Opt_1 
  
  #2
  # Exposure scenario #Inhalation_dermal_Opt_2
  BW_2          = 106.7
  DOSEoral2   = BW_2 * 1000   # ug Oral dose; 
  
  Opt_2    <- ev(ID=1, amt= DOSEoral2, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_2  <- Opt_2 
  
  #3
  # Exposure scenario #Inhalation_dermal_Opt_2
  DOSEoral3   = 50 * 1000   # ug Oral dose; 
  
  Opt_3    <- ev(ID=1, amt= DOSEoral3, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_3  <- Opt_3 
  
  
  #4
  # Exposure scenario #Inhalation_dermal_Opt_2
  DOSEoral4   = 50 * 1000   # ug Oral dose; 
  
  Opt_4    <- ev(ID=1, amt= DOSEoral4, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_4  <- Opt_4 
  
  
  ############
 
  
  # set up the simulation exposure time
  tsamp.B=tgrid(0,72,0.1)  # 48 hours
  
  ########### inhal+dermal 
  out.a <- 
    mod %>% 
    param(pars.data, BW = 73) %>%
    Req(Urine_MHNCH, Urine_MINCH, Urine_cx, Urine_oxo, Urine_CHDA)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_1, tgrid=tsamp.B)
  
  out.a <- cbind.data.frame (Time =out.a$time,
                             MHNCH   = (out.a$Urine_MHNCH)/BW_1,
                             MINCH   = (out.a$Urine_MINCH)/BW_1,
                             cx   = (out.a$Urine_cx)/BW_1,
                             oxo   = (out.a$Urine_oxo)/BW_1,
                             CHDA   = (out.a$Urine_CHDA)/BW_1)
  
  
  out.b <- 
    mod %>% 
    param(pars.data, BW = 106.7) %>%
    Req(Urine_MHNCH, Urine_MINCH, Urine_cx, Urine_oxo, Urine_CHDA)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_2, tgrid=tsamp.B)
  
  out.b <- cbind.data.frame (Time =out.b$time,
                             MHNCH   = (out.b$Urine_MHNCH)/BW_2,
                             MINCH   = (out.b$Urine_MINCH)/BW_2,
                             cx   = (out.b$Urine_cx)/BW_2,
                             oxo   = (out.b$Urine_oxo)/BW_2,
                             CHDA   = (out.b$Urine_CHDA)/BW_2)
  
  
  out.c <- 
    mod %>% 
    param(pars.data) %>%
    Req(Urine_MHNCH, Urine_MINCH, Urine_cx, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_3, tgrid=tsamp.B)
  
  out.c <- cbind.data.frame (Time =out.c$time,
                             MHNCH   = out.c$Urine_MHNCH,
                             MINCH   = out.c$Urine_MINCH,
                             cx   = out.c$Urine_cx,
                             oxo   = out.c$Urine_oxo)
  
  
  out.d <- 
    mod %>% 
    param(pars.data) %>%
    Req(Urine_MHNCH, Urine_MINCH, Urine_cx)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_4, tgrid=tsamp.B)
  
  out.d <- cbind.data.frame (Time =out.d$time,
                             MHNCH   = out.d$Urine_MHNCH,
                             MINCH   = out.d$Urine_MINCH,
                             cx   = out.d$Urine_cx)
  
  
  
  if (pred) return (list("out.a" = out.a,     ## Exposure scenario a
                         "out.b" = out.b,     ## Exposure scenario b 
                         "out.c" = out.c,     ## Exposure scenario c 
                         "out.d" = out.d))     ## Exposure scenario ee 
  
  ## Data mutate with prediction
  out.a  = out.a  [which(out.a$Time  %in% Opt1$time),]
  out.b  = out.b  [which(out.b$Time  %in% Opt2$time),]
  out.c  = out.c  [which(out.c$Time  %in% Opt3$time),]
  out.d  = out.d  [which(out.d$Time  %in% Opt4$time),]
  
  
  #out.a.MHNCH  = out.a  [which(out.a$Time  %in% Opt1_MHNCH$Time),]
 # out.a.MINCH  = out.a  [which(out.a$Time  %in% Opt1_MINCH$Time),]
 # out.a.cx     = out.a  [which(out.a$Time  %in% Opt1_cx$Time),]
 # out.a.oxo    = out.a  [which(out.a$Time  %in% Opt1_oxo$Time),]
 # out.a.CHDA   = out.a  [which(out.a$Time  %in% Opt1_CHDA$Time),]
  
 # out.b.MHNCH  = out.b  [which(out.b$Time  %in% Opt2_MHNCH$Time),]
 # out.b.MINCH  = out.b  [which(out.b$Time  %in% Opt2_MINCH$Time),]
 # out.b.cx     = out.b  [which(out.b$Time  %in% Opt2_cx$Time),]
 # out.b.oxo    = out.b  [which(out.b$Time  %in% Opt2_oxo$Time),]
 # out.b.CHDA   = out.b  [which(out.b$Time  %in% Opt2_CHDA$time),]
  
 # out.c.MHNCH  = out.c  [which(out.c$Time  %in% Opt3_MHNCH$Time),]
 # out.c.MINCH  = out.c  [which(out.c$Time  %in% Opt3_MINCH$Time),]
#  out.c.cx     = out.c  [which(out.c$Time  %in% Opt3_cx$Time),]
 # out.c.oxo    = out.c  [which(out.c$Time  %in% Opt3_oxo$Time),]
 
 # out.d.MHNCH  = out.d  [which(out.d$Time  %in% Opt4_MHNCH$Time),]
 # out.d.MINCH  = out.d  [which(out.d$Time  %in% Opt4_MINCH$Time),]
 # out.d.cx     = out.d  [which(out.d$Time  %in% Opt4_cx$Time),]

 
  ## log.Predition 
  #log.yhat.a.MHNCH      <-log(out.a.MHNCH$MHNCH)
  #log.yhat.a.MINCH      <-log(out.a.MINCH$MINCH)
  #log.yhat.a.cx         <-log(out.a.cx$cx)
  #log.yhat.a.oxo        <-log(out.a.oxo$oxo)
  #log.yhat.a.CHDA       <-log(out.a.CHDA$CHDA)
  
  #log.yhat.b.MHNCH      <-log(out.b.MHNCH$MHNCH)
  #log.yhat.b.MINCH      <-log(out.b.MINCH$MINCH)
 # log.yhat.b.cx         <-log(out.b.cx$cx)
  #log.yhat.b.oxo        <-log(out.b.oxo$oxo)
  #log.yhat.b.CHDA       <-log(out.b.CHDA$CHDA)
  
 # log.yhat.c.MHNCH      <-log(out.c.MHNCH$MHNCH)
 # log.yhat.c.MINCH      <-log(out.c.MINCH$MINCH)
 # log.yhat.c.cx         <-log(out.c.cx$cx)
 # log.yhat.c.oxo        <-log(out.c.oxo$oxo)
  
  #log.yhat.d.MHNCH      <-log(out.d.MHNCH$MHNCH)
 # log.yhat.d.MINCH      <-log(out.d.MINCH$MINCH)
 # log.yhat.d.cx         <-log(out.d.cx$cx)
  
  log.yhat.a.MHNCH      <-log(out.a$MHNCH)
  log.yhat.a.MINCH      <-log(out.a$MINCH)
  log.yhat.a.cx         <-log(out.a$cx)
  log.yhat.a.oxo        <-log(out.a$oxo)
  log.yhat.a.CHDA       <-log(out.a$CHDA)
  
  log.yhat.b.MHNCH      <-log(out.b$MHNCH)
  log.yhat.b.MINCH      <-log(out.b$MINCH)
  log.yhat.b.cx         <-log(out.b$cx)
  log.yhat.b.oxo        <-log(out.b$oxo)
  log.yhat.b.CHDA       <-log(out.b$CHDA)
  
  log.yhat.c.MHNCH      <-log(out.c$MHNCH)
  log.yhat.c.MINCH      <-log(out.c$MINCH)
  log.yhat.c.cx         <-log(out.c$cx)
  log.yhat.c.oxo        <-log(out.c$oxo)
  
  log.yhat.d.MHNCH      <-log(out.d$MHNCH)
  log.yhat.d.MINCH      <-log(out.d$MINCH)
  log.yhat.d.cx         <-log(out.d$cx)
  
 

  
  ## log.Observed data
  log.y.a.MHNCH         <-log(Opt1_MHNCH$MHNCH)
  log.y.a.MINCH         <-log(Opt1_MINCH$MINCH)
  log.y.a.cx            <-log(Opt1_cx$cx)
  log.y.a.oxo           <-log(Opt1_oxo$oxo)
  log.y.a.CHDA          <-log(Opt1_CHDA$CHDA)
  
  log.y.b.MHNCH         <-log(Opt2_MHNCH$MHNCH)
  log.y.b.MINCH         <-log(Opt2_MINCH$MINCH)
  log.y.b.cx            <-log(Opt2_cx$cx)
  log.y.b.oxo           <-log(Opt2_oxo$oxo)
  log.y.b.CHDA          <-log(Opt2_CHDA$CHDA)
  
  log.y.c.MHNCH         <-log(Opt3_MHNCH$MHNCH)
  log.y.c.MINCH         <-log(Opt3_MINCH$MINCH)
  log.y.c.cx            <-log(Opt3_cx$cx)
  log.y.c.oxo           <-log(Opt3_oxo$oxo)
  
  log.y.d.MHNCH         <-log(Opt4_MHNCH$MHNCH)
  log.y.d.MINCH         <-log(Opt4_MINCH$MINCH)
  log.y.d.cx            <-log(Opt4_cx$cx)
  
  
  
  log.yhat        <- c(log.yhat.a.MHNCH, log.yhat.a.MINCH, log.yhat.a.cx, log.yhat.a.oxo,log.yhat.a.CHDA,      
                       log.yhat.b.MHNCH, log.yhat.b.MINCH, log.yhat.b.cx, log.yhat.b.oxo,log.yhat.b.CHDA,       
                       log.yhat.c.MHNCH, log.yhat.c.MINCH, log.yhat.c.cx, log.yhat.c.oxo,        
                       log.yhat.d.MHNCH, log.yhat.d.MINCH, log.yhat.d.cx)   
  
  log.y           <- c(log.y.a.MHNCH, log.y.a.MINCH, log.y.a.cx, log.y.a.oxo,log.y.a.CHDA,      
                       log.y.b.MHNCH, log.y.b.MINCH, log.y.b.cx, log.y.b.oxo,log.y.b.CHDA,       
                       log.y.c.MHNCH, log.y.c.MINCH, log.y.c.cx, log.y.c.oxo,        
                       log.y.d.MHNCH, log.y.d.MINCH, log.y.d.cx)   
  
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
  sig  <- as.numeric (exp(pars[which_sig][2:19]))                 # Coefficient of variation from population variance; sigmal0
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
  CV.sig         = exp(theta.MCMC[which_sig])[2:19]               # Singmal0
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
            niter         = 400000,           ## iteration number 
            jump          = 0.01,             ## jump function generation new parameters distribution using covariate matrix
            lower = c(-1,-3,-3,-4,-4, -8,  -1, -1, -1, -1,  -9,-9,-9, -1,   -1,  -7,  -7,  -7,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf),
            upper = c(+2,+1,+2,+2,+2,0.1, 0.1,0.1,0.1,0.1,  +1,+1,+1, +1, 0.01,0.01,0.01,0.01, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
            prior         = Prior,            ## prior function
            updatecov     = 50,               ## Adaptative Metropolis
            var0          = NULL,             ## initial model variance;
            wvar0         = 0.01,             ## "Weight" for the initial model variance
            ntrydr        = 2,                ## Delayed Rejection
            burninlength  = 200000,           ## number of initial iterations to be removed from output.
            outputlength  = 10000)            ## number of output iterations           
    
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
write.csv(quan.Human,file="C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/Human.summary_pos.csv")
write.csv(MC.H.1,file="C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/Human.pos.csv")
saveRDS(MCMC[[1]],file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/Human.MCMC.rds')
saveRDS(combinedchains,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/Human.comb.rds')


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
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DINCH/PK/DINCH/",
       width = 25, height = 20, units = "cm",dpi=320)

