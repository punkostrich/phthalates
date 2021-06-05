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

## Build mrgsolve-based PK Model
mod <- mcode ("pk", HumanPK.code)


############################################# Model Calibration with levenberg-marquart #############

#Plasma.DEHP <- read.csv(file="~/Dropbox/R/DEHP/2012_Plasma_DEHP.csv")
Plasma.DEHP <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2012_Plasma_DEHP.csv")
names(Plasma.DEHP)=c("Time", "CA_DEHP")

Plasma.MEHP <- read.csv(file="C:/Users/Punkostrich/Dropbox/R//DEHP/2012_Plasma_MEHP.csv")
names(Plasma.MEHP)=c("Time", "CA")

Urine.MEHP  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2012_Urine_MEHP.csv")
names(Urine.MEHP)=c("Time", "AU")

Inhal  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2018_Inhalation.csv")
names(Inhal)=c("Time", "AU")


# prediction function
Pred <- function(pars, pred=FALSE) {
  
  ##' Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  
  # Exposure scenario #1
  BW          = 82
  tinterval   = 24
  TDoses      = 1
  
  
  PDOSEoral1  = 52758                       # ug Oral dose; 
  DOSEoral1   = PDOSEoral1*1
  
  oral1 <- ev(ID=1, amt= DOSEoral1, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex1 <- oral1
  
  # Exposure scenario #2
  
  PDOSEinhal1  = 97.96                       # ug inhal dose; 
  DOSEinhal1   = PDOSEinhal1*1
  
  inhal1 <- ev(ID=1, amt= DOSEinhal1, ii=tinterval, addl=TDoses-1, cmt="AC", tinf = 3, replicate = FALSE)
  ex2 <- inhal1
  
  
  # set up the exposure time
  tsamp=tgrid(0,48,0.1)  # 48 hours
  
  out.1 <- 
    
    mod %>% 
    param(pars, BW = 82)%>% 
    Req(Plasma, Plasma_MEHP, Urine) %>%
    update(atol = 1E-15,maxsteps = 100000000) %>%
    mrgsim_d(data = ex1, tgrid=tsamp)
  out.1<-cbind.data.frame (Time = out.1$time,
                           CA_DEHP  = out.1$Plasma,
                           CA   = out.1$Plasma_MEHP,
                           AU   = (out.1$Urine/DOSEoral1)*100)
  
  
  out.2 <-
    
    mod %>% 
    param(pars,BW = 60) %>%  
    Req(Urine)%>%
    update(atol = 1E-15,maxsteps = 1000000000) %>%
    mrgsim_d(data = ex2, tgrid=tsamp)
  out.2<-cbind.data.frame (Time = out.2$time,
                           AU   = (out.2$Urine/DOSEinhal1)*100)
  
  
  return(list("out.1" = out.1,
              "out.2" = out.2))
  
}

## Cost fuction (from FME pckage) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out <-  Pred (pars)
  cost<- modCost(model=out$out.1, obs= Plasma.DEHP,  weight='std', x="Time")
  cost<- modCost(model=out$out.1, obs= Plasma.MEHP,weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.1, obs= Urine.MEHP,weight='std', x="Time", cost = cost)
  cost<- modCost(model=out$out.2, obs= Inhal,weight='std', x="Time", cost = cost)
  return(cost)
}



###################### Local sensitivity analysis #########################
## Choose the senstiive parameters in the model                          ##
## initial parmaeters######################################################

theta.int <- log(c(
  Vp                  =   23.5,    # Volume of distribution parent compound 
  Vm                  =   23.5,    #: Volume of distribution metabolite compound 
  
  ka                  =   1.6,     #: stomach absorption rate parent 1/h
  #ks   		            =   5,     #: skin Absorption rate to plasma (1/h)
  
  #F                   =   1,       #: transfer ratio from parent to metabolite 
  
  lambda_p            =   0.2,     #: lambda for parent 
  lambda_zm           =   0.2,    #: lambda for metabolite 
  lambda_u            =   0.004   #: lambda for metabolite urine 
  
  #MW                 =   391,     #: g/mol, DEHP molecular mass 
  #Free                =   0.0001,  #: Free fraction of DEHP in plasma 
  #FreeM               =   0.007   #: Free fraction of MEHP in plasma 
  
  #Kvoid               =   0.06974  #: (L/hr), Daily urine volume rate (L/hr)
  
))

# Senstivity function (FME)
# Check the senstive parameters in the modelhttp://ec2-3-135-223-180.us-east-2.compute.amazonaws.com/graphics/plot_zoom_png?width=1200&height=868
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)
Sen=summary(SnsPlasma)
write.csv(Sen,file="Sen.human.csv")

plot(summary(SnsPlasma))





######################### MCMC optimization ##########################
## set up senstivie or necessary parametesr as model input          ##
######################################################################
theta <- log(c(
  Vp                  =   23.5,    # Volume of distribution parent compound 
  Vm                  =   23.5,    #: Volume of distribution metabolite compound 
  
  ka                  =   1.6,     #: stomach absorption rate parent 1/h
  #ks   		            =   0.6,     #: skin Absorption rate to plasma (1/h)
  
  #F                   =   1,       #: transfer ratio from parent to metabolite 
  
  lambda_p            =   0.2,     #: lambda for parent 
  lambda_zm           =   0.2,    #: lambda for metabolite 
  lambda_u            =   0.004   #: lambda for metabolite urine 
  
  #MW                 =   391,     #: g/mol, DEHP molecular mass 
  #Free                =   0.0001,  #: Free fraction of DEHP in plasma 
  #FreeM               =   0.007,   #: Free fraction of MEHP in plasma 
  
  #Kvoid               =   0.06974  #: (L/hr), Daily urine volume rate (L/hr)
  
))

# Model fitting
#upper = c(Inf,Inf,Inf,Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
system.time(Fit<- modFit(f=MCcost, p=theta, method = "Nelder-Mead",lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf),
                         upper = c(6,4,1,Inf,Inf,Inf),
                         control = nls.lm.control(nprint=1)))

summary(Fit)
### para value
exp(Fit$par)

res=MCcost(Fit$par)$residuals$res
sum(res^2)


## Model calibration plot using ggplot2 
Sim.fit.A = Pred (Fit$par)$out.1         ## Time-course concentration profiles using estimated parameters under exposure senaior A
Sim.fit.B = Pred (Fit$par)$out.2         ## Simulaiton of exposure scenaior B

df.sim.A1 = cbind.data.frame (Time=Sim.fit.A$Time, CA_DEHP=Sim.fit.A$CA_DEHP)
df.sim.A2 = cbind.data.frame (Time=Sim.fit.A$Time, CA=Sim.fit.A$CA)
df.sim.A3 = cbind.data.frame (Time=Sim.fit.A$Time, AU=Sim.fit.A$AU)
df.sim.A4 = cbind.data.frame (Time=Sim.fit.B$Time, AU=Sim.fit.B$AU)

## Plot
plot.A1=
  ggplot() +
  geom_line (data = df.sim.A1,aes(Time,CA_DEHP), col="firebrick", lwd=2)+
  geom_point(data = Plasma.DEHP, aes(Time, CA_DEHP),size=2.5) + ylab("Concentration") 


plot.A2=
  ggplot() +
  geom_line (data = df.sim.A2,aes(Time,CA), col="firebrick", lwd=2)+
  geom_point(data = Plasma.MEHP, aes(Time, CA)) + ylab("Concentration") 


plot.A3=
  ggplot() +
  geom_line (data = df.sim.A3,aes(Time,AU), col="firebrick", lwd=2)+
  geom_point(data = Urine.MEHP ,aes(Time, AU),size=2.5) + ylab("Concentration") 


plot.A4=
  ggplot() +
  geom_line (data = df.sim.A4,aes(Time,AU), col="firebrick", lwd=2)+
  geom_point(data = Inhal,aes(Time, AU),size=2.5) + ylab("Concentration") 


plot.A5=
  ggplot() +
  geom_line (data = df.sim.A1,aes(Time,CA_DEHP), col="firebrick", lwd=2)+
  geom_point(data = Plasma.DEHP, aes(Time, CA_DEHP),col="firebrick",size=2.5) + ylab("Concentration") +
  geom_line (data = df.sim.A2,aes(Time,CA), col="blue", lwd=2)+
  geom_point(data = Plasma.MEHP, aes(Time, CA), col="blue",) + 
  geom_line (data = df.sim.A3,aes(Time,AU), col="green", lwd=2)+
  geom_point(data = Urine.MEHP ,aes(Time, AU),col="green",size=2.5) +
  geom_line (data = df.sim.A4,aes(Time,AU), col="yellow", lwd=2)+
  geom_point(data = Inhal,aes(Time, AU),size=2.5,col="yellow") 

plot.A1
plot.A2
plot.A3
plot.A4
plot.A5



############################################# Model Calibration with MCMC ###################################################
## One data sets was used in model evaluation                                                                               #                                                       #
#############################################################################################################################

# input callibrated data set; (Koch et al., 2004)
#Human.2 <- read.csv(file="~/Dropbox/R/DEHP/2004_Plasma.csv")
Human.2 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2004_Plasma.csv")
names(Human.2)=c("Time", "CA2")
DOSEoral_2    = 48500 

Human.3 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2004_Urine.csv")
names(Human.3)=c("Time", "AU3")
DOSEoral_3    = 48500 

## Fixed the physiological parameters;
## Input a new initial parameters
## Population mean and model error (sig2)

theta.MCMC <- log(c(
  Vp                  =   7.44537489,    # Volume of distribution parent compound 
  Vm                  =   0.19722783,    #: Volume of distribution metabolite compound 
  
  ka                  =   2.71825288,     #: stomach absorption rate parent 1/h
  
  lambda_p            =   0.14658602,     #: lambda for parent 
  lambda_zm           =   1.49905196,    #: lambda for metabolite 
  lambda_u            =   0.03292043,   #: lambda for metabolite urine 
  
 
  
  sig2                =   1,      ## Model residuals; mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of most parametesr is 0.3, some parameters is 0.5 due to possible variation)
  sig_Vp                = 1,                ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_Vm                = 1,
  sig_ka                = 1,
  
  sig_lambda_p          = 1,
  sig_lambda_zm         = 1,
  sig_lambda_u          = 1
 
))


which_sig <- grep("sig", names(theta.MCMC))

## Maximum likelihood estimation (MLE) fuction for MCMC
mcmc.fun <- function (pars, pred=FALSE){
  
  ## Get out of log domain
  pars.data <- lapply(pars [-which_sig],exp)
  names(pars.data) <- names(pars.data)
  
  # Exposure scenario for single low dose
  
  BW                    = 75                            ## kg, body weight
  tinterval             = 24                            ## hr, Time interval
  TDoses                = 1                             ## Dose times
  DOSEoral_2004         = 48500                          ## ug
  
  
  oral_2004 <- ev(ID=1, amt= DOSEoral_2004, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.B <- oral_2004
  
  
  # set up the simulation exposure time
  tsamp.B=tgrid(0,48,0.1)  # 48 hours
  
  out.B <- 
    mod %>% 
    param(pars.data, BW = 75) %>%
    Req(Plasma_MEHP)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex.B, tgrid=tsamp.B)
  
  out.B <- cbind.data.frame (  Time =out.B$time,
                               CA2  =out.B$Plasma_MEHP)
  
  
  out.C <- 
    mod %>% 
    param(pars.data, BW = 75) %>%
    Req(Urine)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex.B, tgrid=tsamp.B)
  
  out.C <- cbind.data.frame (  Time =out.C$time,
                               AU3  =out.C$Urine)
  
  if (pred) return (list("out.B" = out.B,      ## Exposure scenario B 
                         "out.C" = out.C))     ## Exposure scenario C 
  
  ## Data mutate with prediction
  out.B = out.B [which(out.B$Time %in% Human.2$Time),]
  out.C = out.C [which(out.C$Time %in% Human.3$Time),]
  
  ## log.Predition 
  log.yhat.B      <-log(out.B$CA2)
  log.yhat.C      <-log(out.C$AU3)
  
  ## log.Observed data
  log.y.B         <-log(Human.2$CA2)
  log.y.C         <-log(Human.3$AU3)
  
  log.yhat        <- c(log.yhat.B,log.yhat.C)
  log.y           <- c(log.y.B,log.y.C)
  
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
  sig  <- as.numeric (exp(pars[which_sig][2:7]))                 # Coefficient of variation from population variance; sigmal0
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
  CV.sig         = exp(theta.MCMC[which_sig])[2:7]               # Singmal0
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
            niter         = 50000,           ## iteration number 
            jump          = 0.01,             ## jump function generation new parameters distribution using covariate matrix
            lower = c(0.01,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf),
            prior         = Prior,            ## prior function
            updatecov     = 50,               ## Adaptative Metropolis
            var0          = NULL,             ## initial model variance;
            wvar0         = 0.01,             ## "Weight" for the initial model variance
            ntrydr        = 2,                ## Delayed Rejection
            burninlength  = 20000,           ## number of initial iterations to be removed from output.
            outputlength  = 2000)            ## number of output iterations           
    
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
write.csv(quan.Human,file="C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/Human.summary_pos.csv")
write.csv(MC.H.1,file="C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/Human.pos.csv")
saveRDS(MCMC[[1]],file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/Human.MCMC.rds')
saveRDS(combinedchains,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/Human.comb.rds')


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
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/",
       width = 25, height = 20, units = "cm",dpi=320)

