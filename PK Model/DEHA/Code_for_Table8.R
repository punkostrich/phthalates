## loading R packages
library(magrittr)   # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(reshape2)   # melt function to reshape the table
library(tidyverse)  # Needed for the pipe %>% operator
library(FME)
## Loading human mrgModelfitting code
Human.MCMC  <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/Human.MCMC.rds")


## Population mean and model error (sig2)

theta.Human <- log(c(
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
  
  
  sig2                =   1,      
  
  
  sig_Vp                = 1,               
  sig_Vm                = 1,
  sig_ka                = 1,
  sig_lambda_p           = 1,
  
  sig_lambda_zm      = 1,
  sig_lambda_z_OH        = 1,
  sig_lambda_z_cx        = 1,
  sig_lambda_z_oxo       = 1,
  
  sig_lambda_u_OH        = 1,
  sig_lambda_u_cx        = 1,
  sig_lambda_u_oxo       = 1,
  
  sig_Fr_OH              = 1,
  sig_Fr_oxo             = 1
  
))



which_sig <- grep("sig", names(theta.Human))

### save theta names #############
theta.names <- c(
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
  
  
  sig2                =   1,    
  
  
  sig_Vp                = 1,             
  sig_Vm                = 1,
  sig_ka                = 1,
  sig_lambda_p           = 1,
  
  sig_lambda_zm      = 1,
  sig_lambda_z_OH        = 1,
  sig_lambda_z_cx        = 1,
  sig_lambda_z_oxo       = 1,
  
  sig_lambda_u_OH        = 1,
  sig_lambda_u_cx        = 1,
  sig_lambda_u_oxo       = 1,
  
  sig_Fr_OH              = 1,
  sig_Fr_oxo             = 1
  
)

saveRDS(names(theta.names),file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/theta.names.RDS')

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


### save theta names #############
theta.names <- c(
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
  
 
  sig2                =   1,     
  sig_Vp                = 1,               
  sig_Vm                = 1,
  sig_ka                = 1,
  sig_lambda_p           = 1,
  
  sig_lambda_zm      = 1,
  sig_lambda_z_OH        = 1,
  sig_lambda_z_cx        = 1,
  sig_lambda_z_oxo       = 1,
  
  sig_lambda_u_OH        = 1,
  sig_lambda_u_cx        = 1,
  sig_lambda_u_oxo       = 1,
  
  sig_Fr_OH              = 1,
  sig_Fr_oxo             = 1
  
  
  #sig_Free              = 0.5,
  #sig_FreeM             = 0.5,
  #sig_Kvoid             = 0.5
  
)

saveRDS(names(theta.names),file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/theta.names.DEHA.RDS')

