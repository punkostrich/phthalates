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

############### DEHA #############
############################################# Model Calibration with levenberg-marquart #############

#Plasma.DEHP <- read.csv(file="~/Dropbox/R/DEHP/2012_Plasma_DEHP.csv")

Cal1 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHA/Cal1.csv")
Cal1_OH <- cbind.data.frame (Time = Cal1$time_oh,
                                OH   = Cal1$mass_oh)
Cal1_oxo <- cbind.data.frame (Time = Cal1$time_oxo,
                                oxo   = Cal1$mass_oxo)
Cal1_cx <- cbind.data.frame (Time = Cal1$time_cx,
                                cx   = Cal1$mass_cx)



# prediction function
Pred <- function(pars, pred=FALSE) {
  
  ##' Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  
  ### shared parameter
  tinterval   = 24
  TDoses      = 1
  
  #1
  DOSEoral1   = 10 * 1000   # ug Oral dose; 
  
  Cal_1    <- ev(ID=1, amt= DOSEoral1, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_cal_1  <- Cal_1 
  

 
  # set up the exposure time
  tsamp=tgrid(0,48,0.01)  # 48 hours
  
  out.1 <- 
    
    mod %>% 
    param(pars)%>% 
    Req(Urine_OH, Urine_cx, Urine_oxo) %>%
    update(atol = 1E-15,maxsteps = 100000000) %>%
    mrgsim_d(data = ex_cal_1, tgrid=tsamp)
  out.1<-cbind.data.frame (Time = out.1$time,
                           OH   = (out.1$Urine_OH),
                           cx   = (out.1$Urine_cx),
                           oxo   = (out.1$Urine_oxo))
  
  
  
  return(list("out.1" = out.1))
  
}

## Cost fuction (from FME pckage) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out <-  Pred (pars)
  cost<- modCost(model=out$out.1, obs= Cal1_OH, weight='std', x="Time")
  cost<- modCost(model=out$out.1, obs= Cal1_cx, weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Cal1_oxo, weight='std', x="Time", cost = cost)
  return(cost)
}



###################### Local sensitivity analysis #########################
## Choose the senstiive parameters in the model                          ##
## initial parmaeters######################################################

theta.int <- log(c(
  Vp                  =   5.3,    # : Volume of distribution parent compound L/kg bw
  Vm                  =  0.5,    #: Volume of distribution metabolite compound 
  
  ka                  =   0.2,     #: stomach absorption rate parent 1/h
  
  BW                  =   82.3,    #: kg, Bodyweight (EPA Factors Handbook, 2011)
  
  lambda_p            =   0.08,    #: lambda for parent (h-1)
  
  
  lambda_zm           =   0.2,       #: h-1
  lambda_z_OH         =   0.4,       #: h-1
  lambda_z_cx         =   0.4,       #: h-1
  lambda_z_oxo        =   0.4,       #: h-1
  
  
  lambda_u_OH         =   2.8E-4,    # : h-1
  lambda_u_cx         =   8E-4,     #: h-1
  lambda_u_oxo        =   2E-4,     #: h-1
  
  
  Fr_OH               =   0.5,    # : h-1
  Fr_oxo              =   0.5    # : h-1
  
))

# Senstivity function (FME)
# Check the senstive parameters in the modelhttp://ec2-3-135-223-180.us-east-2.compute.amazonaws.com/graphics/plot_zoom_png?width=1200&height=868
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)
Sen=summary(SnsPlasma)
write.csv(Sen,file="C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/Sen.human.csv")

plot(summary(SnsPlasma))





######################### MCMC optimization ##########################
## set up senstivie or necessary parametesr as model input          ##
######################################################################
theta <- log(c(
  #Vp                  =   5.3,    # : Volume of distribution parent compound L/kg bw
  #Vm                  =  0.5,    #: Volume of distribution metabolite compound 
  
  ka                  =   0.2,     #: stomach absorption rate parent 1/h
  
  #BW                  =   82.3,    #: kg, Bodyweight (EPA Factors Handbook, 2011)
  
  lambda_p            =   0.08,    #: lambda for parent (h-1)
  
  
  lambda_zm           =   0.2,       #: h-1
  lambda_z_OH         =   0.4,       #: h-1
  lambda_z_cx         =   0.4,       #: h-1
  lambda_z_oxo        =   0.4,       #: h-1
  
  
  lambda_u_OH         =   2.8E-4,    # : h-1
  lambda_u_cx         =   8E-4,     #: h-1
  lambda_u_oxo        =   2E-4,     #: h-1
  
  
  Fr_OH               =   0.5,    # : h-1
  Fr_oxo               =   0.5    # : h-1
))

# Model fitting
#upper = c(Inf,Inf,Inf,Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
system.time(Fit<- modFit(f=MCcost, p=theta, method = "Nelder-Mead",lower = c(-5,-5,-5,-1,-1,-1,-10,-10,-10,-5,-5),
                         upper = c(3,3,3,3,5,5,3,3,3,0.01,0.01),
                         control = nls.lm.control(nprint=1)))

summary(Fit)
### para value
exp(Fit$par)

res=MCcost(Fit$par)$residuals$res
sum(res^2)


## Model calibration plot using ggplot2 
Sim.fit.A = Pred (Fit$par)$out.1         ## Time-course concentration profiles using estimated parameters under exposure senaior 1

df.sim.A1  = cbind.data.frame (Time=Sim.fit.A$Time, OH=Sim.fit.A$OH)
df.sim.A2  = cbind.data.frame (Time=Sim.fit.A$Time, cx=Sim.fit.A$cx)
df.sim.A3  = cbind.data.frame (Time=Sim.fit.A$Time, oxo=Sim.fit.A$oxo)

## Plot
plot.A1=
  ggplot() +
  geom_line (data = df.sim.A1,aes(Time,OH), col="firebrick", lwd=2)+
  geom_point(data = Cal1_OH, aes(Time, OH),size=2.5) + ylab("Concentration") 
plot.A1


plot.A2=
  ggplot() +
  geom_line (data = df.sim.A2,aes(Time,cx), col="firebrick", lwd=2)+
  geom_point(data = Cal1_cx,aes(Time, cx),size=2.5) + ylab("Concentration") 
plot.A2


plot.A3=
  ggplot() +
  geom_line (data = df.sim.A3,aes(Time,oxo), col="firebrick", lwd=2)+
  geom_point(data = Cal1_oxo,aes(Time, oxo),size=2.5) + ylab("Concentration") 
plot.A3



plot.A1
plot.A2
plot.A3






############################################# Model Calibration with MCMC ###################################################
## EIGHT data sets was used in model evaluation                                                                               #                                                       #
#############################################################################################################################

# input callibrated data set;
#Human.2 <- read.csv(file="~/Dropbox/R/DEHP/2004_Plasma.csv")
Opt1 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHA/Opt1A.csv")
Opt1_OH <- cbind.data.frame (Time = Opt1$time,
                             OH   = Opt1$mass_oh)
Opt1_oxo <- cbind.data.frame (Time = Opt1$time,
                              oxo   = Opt1$mass_oxo)
Opt1_cx <- cbind.data.frame (Time = Opt1$time,
                             cx   = Opt1$mass_cx)
Opt1_OH <- na.omit(Opt1_OH)
Opt1_oxo <- na.omit(Opt1_oxo)
Opt1_cx <- na.omit(Opt1_cx)


Opt2 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHA/Opt2A.csv")
Opt2_OH <- cbind.data.frame (Time = Opt2$time,
                             OH   = Opt2$mass_oh)
Opt2_oxo <- cbind.data.frame (Time = Opt2$time,
                              oxo   = Opt2$mass_oxo)
Opt2_cx <- cbind.data.frame (Time = Opt2$time,
                             cx   = Opt2$mass_cx)
Opt2_OH <- na.omit(Opt2_OH)
Opt2_oxo <- na.omit(Opt2_oxo)
Opt2_cx <- na.omit(Opt2_cx)


Opt3 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHA/Opt3A.csv")
Opt3_OH <- cbind.data.frame (Time = Opt3$time,
                             OH   = Opt3$mass_oh)
Opt3_oxo <- cbind.data.frame (Time = Opt3$time,
                              oxo   = Opt3$mass_oxo)
Opt3_cx <- cbind.data.frame (Time = Opt3$time,
                             cx   = Opt3$mass_cx)
Opt3_OH <- na.omit(Opt3_OH)
Opt3_oxo <- na.omit(Opt3_oxo)
Opt3_cx <- na.omit(Opt3_cx)

Opt3_OH  <- Opt3_OH[Opt3_OH$Time != 0, ]
Opt3_oxo <- Opt3_oxo[Opt3_oxo$Time != 0, ]
Opt3_cx  <- Opt3_cx[Opt3_cx$Time != 0, ]


Opt4 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHA/Opt4A.csv")
Opt4_OH <- cbind.data.frame (Time = Opt4$time,
                             OH   = Opt4$mass_oh)
Opt4_oxo <- cbind.data.frame (Time = Opt4$time,
                              oxo   = Opt4$mass_oxo)
Opt4_cx <- cbind.data.frame (Time = Opt4$time,
                             cx   = Opt4$mass_cx)
Opt4_OH <- na.omit(Opt4_OH)
Opt4_oxo <- na.omit(Opt4_oxo)
Opt4_cx <- na.omit(Opt4_cx)



## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)

theta.MCMC <- log(c(
  Vp                 =   5.3,    # : Volume of distribution parent compound L/kg bw
  Vm                 =  0.5,    #: Volume of distribution metabolite compound 
  
  ka                  =   4.206204444 ,     #: stomach absorption rate parent 1/h
  
  
  lambda_p            =   7.729299285,    #: lambda for parent (h-1)
  
  
  lambda_zm           =   2.900184044,       #: h-1
  lambda_z_OH         =   0.400726425,       #: h-1
  lambda_z_cx         =   0.396226689,       #: h-1
  lambda_z_oxo        =   0.398307769,       #: h-1
  
  
  lambda_u_OH         =   0.001094498,    # : h-1
  lambda_u_cx         =   0.001157938,     #: h-1
  lambda_u_oxo        =   0.001577265,     #: h-1
  
  
  Fr_OH               =   0.280284167,    # : h-1
  Fr_oxo              =   0.522078175 ,    # : h-1
  
  
  sig2                =   1,      ## Model residuals; mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of most parametesr is 0.3, some parameters is 0.5 due to possible variation)
  sig_Vp                = 1,                ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_Vm                = 1,
  
  sig_ka                 = 1, 
  sig_lambda_p           = 1,
  
  sig_lambda_zm          = 1,
  sig_lambda_z_OH        = 1,
  sig_lambda_z_cx        = 1,
  sig_lambda_z_oxo       = 1,
  
  sig_lambda_u_OH        = 1,
  sig_lambda_u_cx        = 1,
  sig_lambda_u_oxo       = 1,
  
  sig_Fr_OH              = 2,
  sig_Fr_oxo             = 1
  
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
  DOSEoral1   = 10 * 1000   # ug Oral dose; 
  
  Opt_1    <- ev(ID=1, amt= DOSEoral1, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_1  <- Opt_1 
  
  #2
  # Exposure scenario #Inhalation_dermal_Opt_2
  DOSEoral2   = 10 * 1000   # ug Oral dose; 
  
  Opt_2    <- ev(ID=1, amt= DOSEoral2, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_2  <- Opt_2 
  
  #3
  # Exposure scenario #Inhalation_dermal_Opt_2
  DOSEoral3   = 10 * 1000   # ug Oral dose; 
  
  Opt_3    <- ev(ID=1, amt= DOSEoral3, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_3  <- Opt_3 
  
  
  #4
  # Exposure scenario #Inhalation_dermal_Opt_2
  DOSEoral4   = 10 * 1000   # ug Oral dose; 
  
  Opt_4    <- ev(ID=1, amt= DOSEoral4, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_Opt_4  <- Opt_4 
  
  
  ############
 
  
  # set up the simulation exposure time
  tsamp.B=tgrid(0,48,0.01)  # 48 hours
  
  ########### inhal+dermal 
  out.a <- 
    mod %>% 
    param(pars.data) %>%
    Req(Urine_OH, Urine_cx, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_1, tgrid=tsamp.B)
  
  out.a <- cbind.data.frame (Time = out.a$time,
                             OH   = (out.a$Urine_OH),
                             cx   = (out.a$Urine_cx),
                             oxo   = (out.a$Urine_oxo))
  
  
  out.b <- 
    mod %>% 
    param(pars.data) %>%
    Req(Urine_OH, Urine_cx, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_2, tgrid=tsamp.B)
  
  out.b <- cbind.data.frame (Time = out.b$time,
                             OH   = (out.b$Urine_OH),
                             cx   = (out.b$Urine_cx),
                             oxo  = (out.b$Urine_oxo))
  
  
  out.c <- 
    mod %>% 
    param(pars.data) %>%
    Req(Urine_OH, Urine_cx, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_3, tgrid=tsamp.B)
  
  out.c <- cbind.data.frame (Time = out.c$time,
                             OH   = (out.c$Urine_OH),
                             cx   = (out.c$Urine_cx),
                             oxo  = (out.c$Urine_oxo))
  
  
  out.d <- 
    mod %>% 
    param(pars.data) %>%
    Req(Urine_OH, Urine_cx, Urine_oxo)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_Opt_4, tgrid=tsamp.B)
  
  out.d <- cbind.data.frame (Time = out.d$time,
                             OH   = (out.d$Urine_OH),
                             cx   = (out.d$Urine_cx),
                             oxo  = (out.d$Urine_oxo))
  
  
  
  if (pred) return (list("out.a" = out.a,     ## Exposure scenario a
                         "out.b" = out.b,     ## Exposure scenario b 
                         "out.c" = out.c,     ## Exposure scenario c 
                         "out.d" = out.d))     ## Exposure scenario ee 
  
  ## Data mutate with prediction
  out.a  = out.a  [which(out.a$Time  %in% Opt1$time),]
  out.b  = out.b  [which(out.b$Time  %in% Opt2$time),]
  out.c  = out.c  [which(out.c$Time  %in% Opt3$time),]
  out.d  = out.d  [which(out.d$Time  %in% Opt4$time),]

  
  
  ## log.Predition 
  log.yhat.a.OH         <-log(out.a$OH)
  log.yhat.a.cx         <-log(out.a$cx)
  log.yhat.a.oxo        <-log(out.a$oxo)

  log.yhat.b.OH         <-log(out.b$OH)
  log.yhat.b.cx         <-log(out.b$cx)
  log.yhat.b.oxo        <-log(out.b$oxo)

  log.yhat.c.OH         <-log(out.c$OH)
  log.yhat.c.cx         <-log(out.c$cx)
  log.yhat.c.oxo        <-log(out.c$oxo)
  
  log.yhat.d.OH         <-log(out.d$OH)
  log.yhat.d.cx         <-log(out.d$cx)
  log.yhat.d.oxo        <-log(out.d$oxo)
 

  
  ## log.Observed data
  log.y.a.OH            <-log(Opt1_OH$OH)
  log.y.a.cx            <-log(Opt1_cx$cx)
  log.y.a.oxo           <-log(Opt1_oxo$oxo)
  
  log.y.b.OH            <-log(Opt2_OH$OH)
  log.y.b.cx            <-log(Opt2_cx$cx)
  log.y.b.oxo           <-log(Opt2_oxo$oxo)
  
  log.y.c.OH            <-log(Opt3_OH$OH)
  log.y.c.cx            <-log(Opt3_cx$cx)
  log.y.c.oxo           <-log(Opt3_oxo$oxo)
  
  log.y.d.OH            <-log(Opt4_OH$OH)
  log.y.d.cx            <-log(Opt4_cx$cx)
  log.y.d.oxo           <-log(Opt4_oxo$oxo)
  
  
  
  log.yhat        <- c(log.yhat.a.OH, log.yhat.a.cx, log.yhat.a.oxo,
                       log.yhat.b.OH, log.yhat.b.cx, log.yhat.b.oxo,  
                       log.yhat.c.OH, log.yhat.c.cx, log.yhat.c.oxo,        
                       log.yhat.d.OH, log.yhat.d.cx, log.yhat.d.oxo)   
  
  log.y           <- c(log.y.a.OH, log.y.a.cx, log.y.a.oxo,    
                       log.y.b.OH, log.y.b.cx, log.y.b.oxo,       
                       log.y.c.OH, log.y.c.cx, log.y.c.oxo,        
                       log.y.d.OH, log.y.d.cx, log.y.d.oxo)   
  
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
  sig  <- as.numeric (exp(pars[which_sig][2:14]))                 # Coefficient of variation from population variance; sigmal0
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
  CV.sig         = exp(theta.MCMC[which_sig])[2:14]               # Singmal0
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
            niter         = 200000,           ## iteration number 
            jump          = 0.01,             ## jump function generation new parameters distribution using covariate matrix
            lower = c(-1,-3,-4,  -4,  -4,  -4,  -4,  -4,  -10, -10, -10,   -5,  -5,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf),
            upper = c(+2,+1,1.5,2.1, 1.5,-0.1,-0.1,-0.1, 0.01,0.01,0.01, -0.5,-0.1,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
            prior         = Prior,            ## prior function
            updatecov     = 50,               ## Adaptative Metropolis
            var0          = NULL,             ## initial model variance;
            wvar0         = 0.01,             ## "Weight" for the initial model variance
            ntrydr        = 2,                ## Delayed Rejection
            burninlength  = 100000,           ## number of initial iterations to be removed from output.
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
write.csv(quan.Human,file="C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/Human.summary_pos.csv")
write.csv(MC.H.1,file="C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/Human.pos.csv")
saveRDS(MCMC[[1]],file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/Human.MCMC.rds')
saveRDS(combinedchains,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/Human.comb.rds')


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
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHA/PK/DEHP/",
       width = 25, height = 20, units = "cm",dpi=320)

