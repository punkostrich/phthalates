## loading R packages
library(magrittr)   # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(reshape2)   # melt function to reshape the table
library(tidyverse)  # Needed for the pipe %>% operator

## Loading human mrgModelfitting code
Human.MCMC  <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/Human.MCMC.rds")


## Population mean and model error (sig2)

theta.Human <- log(c(
  Vp                  =   5.3,     # Volume of distribution parent compound L/kg bw
  Vm                  =   0.05,    #: Volume of distribution metabolite compound 
  
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



which_sig <- grep("sig", names(theta.Human))

### save theta names #############
theta.names <- c(
  Vp                  =   5.3,     # Volume of distribution parent compound L/kg bw
  Vm                  =   0.05,    #: Volume of distribution metabolite compound 
  
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
  
)

saveRDS(names(theta.names),file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/theta.names.RDS')

##################################

########################################## Table 3 ################################ 
# Mean (2.5%, 97.5%) of prior and posterior distribution for DEHP                 #
###################################################################################

## Prior parameters
## Population mean (u)

## M value: Mean of population mean u
## S vaelu: SD of population mean u: set the CV equal to 50%

mean.human  = exp(theta.Human  [-which_sig]) # M-value for human
human.sd    = mean.human*0.5                 # S-value for human

## 50%, 2.5% and 97.5% for prior population mean u
a.human  = qnorm(0.025, mean = mean.human, sd = human.sd)
b.human  = qnorm(0.975, mean = mean.human, sd = human.sd)
#c.human  = qnorm(0.5,mean = mean.human, sd = human.sd)

## Median (2.5%, 97.5%) for Posterior parameters of population mean
## Estimate the quantile of posterior parameters
quan.Human  = summary(as.mcmc(Human.MCMC$pars))$quantiles  

exp(quan.Human)


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

saveRDS(names(theta.MCMC),file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/theta.names.DINCH.RDS')
