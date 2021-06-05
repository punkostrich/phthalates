## loading R packages
library(magrittr)   # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(reshape2)   # melt function to reshape the table
library(tidyverse)  # Needed for the pipe %>% operator

## Loading human, rat, mouse, monkey MCMC data, from line 557 of the mrgModelfitting code
Human.MCMC  <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/Human.MCMC.rds")


## Population mean and model error (sig2)

theta.Human <- log(c(
  #Vp_b                  =   140,    # Volume of distribution parent compound 
  #Vm_b                  =   40,    #: Volume of distribution metabolite compound 
  #ka                  =   0.23,     #: stomach absorption rate parent 1/h
  lambda_p            =   0.5,     #: lambda for parent 
  lambda_zm           =   1,    #: lambda for metabolite 
  lambda_u            =   0.02,   #: lambda for metabolite urine 
  
  sig2                =   0.5,      ## Model residuals; mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of most parametesr is 0.3, some parameters is 0.5 due to possible variation)
  #sig_Vp                = 1.5,                ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
 # sig_Vm                = 0.5,
  #sig_ka                = 0.5,
  sig_lambda_p          = 0.5,
  sig_lambda_zm         = 0.5,
  sig_lambda_u          = 0.5
 
))



which_sig <- grep("sig", names(theta.Human))

### save theta names #############
theta.names <- c(
  #Vp_b                  =   140,    # Volume of distribution parent compound 
 # Vm_b                  =   40,    #: Volume of distribution metabolite compound 
  #ka                  =   0.23,     #: stomach absorption rate parent 1/h
  lambda_p            =   0.5,     #: lambda for parent 
  lambda_zm           =   1,    #: lambda for metabolite 
  lambda_u            =   0.02,   #: lambda for metabolite urine 
  
  sig2                =   0.5,      ## Model residuals; mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of most parametesr is 0.3, some parameters is 0.5 due to possible variation)
  #sig_Vp                = 1.5,                ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
 # sig_Vm                = 0.5,
 # sig_ka                = 0.5,
  sig_lambda_p          = 0.5,
  sig_lambda_zm         = 0.5,
  sig_lambda_u          = 0.5
  
)

saveRDS(names(theta.names),file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/theta.names.RDS')

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



theta.names.DEP <- c(
  Vp                  =   140,    # Volume of distribution parent compound 
  Vm                  =   40,    #: Volume of distribution metabolite compound 
  ka                  =   0.23,     #: stomach absorption rate parent 1/h
  lambda_p            =   0.5,     #: lambda for parent 
  lambda_zm           =   1,    #: lambda for metabolite 
  lambda_u            =   0.02,   #: lambda for metabolite urine 
  
  sig2                =   0.5,      ## Model residuals; mostly between 0.3 and 0.5 (corresponding to coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of most parametesr is 0.3, some parameters is 0.5 due to possible variation)
  sig_Vp                = 1.5,                ## Default value of 0.3 and 0.5 was used to represent a moderate level of variation (Hack et al., 2006; Chiu et al., 2009)
  sig_Vm                = 0.5,
  sig_ka                = 0.5,
  sig_lambda_p          = 0.5,
  sig_lambda_zm         = 0.5,
  sig_lambda_u          = 0.5
  
)

saveRDS(names(theta.names.DEP),file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/theta.names.DEP.RDS')
