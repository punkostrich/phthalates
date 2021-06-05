########################################## Fig. 6 ################################ 
# circular plot of sensitivity analysis                                           #
###################################################################################
## loading R packages
library(magrittr)   # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)      # Needed for the pipe %>% operator
library(mrgsolve)   # Needed to run the main PBPK code
library(reshape)    # melt function to reshape the table
library(ggplot2)    # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(grid)       # for plotting the figure
library(lattice)    # for plotting the figure

rm(list=ls())

## Input PBPK model
DEP.PK.code     <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/humanPK.RDS')
DEHA.PK.code    <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/humanPK.RDS')
DEHP.PK.code    <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/humanPK.RDS')
DINCH.PK.code   <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/humanPK.RDS')
DPHP.PK.code    <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/humanPK.RDS')

# load pk model
mod.DEP         <- mcode ("DEP_PK",   DEP.PK.code)
mod.DEHA        <- mcode ("DEHA_PK",  DEHA.PK.code)
mod.DEHP        <- mcode ("DEHP_PK",  DEHP.PK.code)
mod.DINCH       <- mcode ("DINCH_PK", DINCH.PK.code)
mod.DPHP        <- mcode ("DPHP_PK",  DPHP.PK.code)

## Loading MCMC data
DEP.MCMC         <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/Human.MCMC.rds")
DEHA.MCMC        <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/Human.MCMC.rds")
DEHP.MCMC        <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/Human.MCMC.rds")
DINCH.MCMC       <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/Human.MCMC.rds")
DPHP.MCMC        <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/Human.MCMC.rds")

## loading the theta names
DEP.names        <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/theta.names.DEP.RDS")
DEHA.names       <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/theta.names.DEHA.RDS")
DEHP.names       <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/theta.names.DEHP.RDS")
DINCH.names      <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/theta.names.DINCH.RDS")
DPHP.names       <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/theta.names.DPHP.RDS")

which_sig.DEP          <- grep("sig", DEP.names)
which_sig.DEHA         <- grep("sig", DEHA.names)
which_sig.DEHP         <- grep("sig", DEHP.names)
which_sig.DINCH        <- grep("sig", DINCH.names)
which_sig.DPHP         <- grep("sig", DPHP.names)


## Sensitivity analysis 3 CHEMICAL DP, DEHP, DPHP
Pred <- function (pars.DEP, pars.DEHP){
  
  which_sig.DEP          <- grep("sig", DEP.names)
  which_sig.DEHP         <- grep("sig", DEHP.names)
  
  ## Get out of log domain
  pars.DEP    <- lapply(pars.DEP  [-which_sig.DEP],exp) # pars.mouse is the entire parameter set; pars.mouse [-which_sig] means to keep parameters with "sig" only, and then do exp transformation, then reassign to pars.mouse
  pars.DEHP   <- lapply(pars.DEHP [-which_sig.DEHP],exp)
  
  
  ## Repeat dose exposure scenario: 
   
  BW                = 82.3                                ## human body weight
  
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
  
  PDOSEoral.DEP     = 1                                ## ug/kg/d oral dose
  PDOSEoral.DEHP    = 0.1                                ## ug/kg/d
  #PDOSEoral.DPHP    = 0.1                                ## ug/kg/d
   
  DOSEoral.DEP     = PDOSEoral.DEP   * BW                                ## ug/kg/d
  DOSEoral.DEHP    = PDOSEoral.DEHP  * BW                                 ## ug/kg/d
  #DOSEoral.DPHP    = PDOSEoral.DPHP  * BW                                ## ug/kg/d
  
  PDOSEother.DEP     = 1                                ## ug/kg/d other exposure dose
  PDOSEother.DEHP    = 1                                ## ug/kg/d
  #PDOSEother.DPHP    = 1                                ## ug/kg/d
  
  DOSEother.DEP     = PDOSEother.DEP   * BW                                ## ug/kg/d
  DOSEother.DEHP    = PDOSEother.DEHP  * BW                                 ## ug/kg/d
  #DOSEother.DPHP    = PDOSEother.DPHP  * BW                                ## ug/kg/d

  
  ############ DOSE
  ex.oral.DEP          <- ev(ID=1, amt= DOSEoral.DEP, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.other.DEP         <- ev(ID=1, amt= DOSEother.DEP, ii=tinterval, addl=TDoses-1, cmt="AC", tinf = 24, replicate = FALSE)
  events.DEP           <- ex.oral.DEP + ex.other.DEP
  #events.DEP           <- ex.oral.DEP 
  
  ex.oral.DEHP         <- ev(ID=1, amt= DOSEoral.DEHP, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.other.DEHP        <- ev(ID=1, amt= DOSEother.DEHP, ii=tinterval, addl=TDoses-1, cmt="AC", tinf = 24, replicate = FALSE)
  events.DEHP          <- ex.oral.DEHP + ex.other.DEHP
  #events.DEHP           <- ex.oral.DEHP
  
  #ex.oral.DPHP         <- ev(ID=1, amt= DOSEoral.DPHP, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
 # ex.other.DPHP        <- ev(ID=1, amt= DOSEother.DPHP, ii=tinterval, addl=TDoses-1, cmt="AC", tinf = 24, replicate = FALSE)
  #events.DPHP          <- ex.oral.DPHP + ex.other.DPHP
 # events.DPHP           <- ex.oral.DPHP
  
  ## set up the exposure time
  tsamp           = tgrid(0,tinterval*(TDoses-1)+24*7,1)          ## 2 weeks and simulated for 24 * 7 hours after dosing at one-hour resolution
  
    
  ## Get a prediction
   out.DEP <- 
    mod.DEP  %>%
    param(pars.DEP) %>%
    Req(AUC_CA,AUC_CAM)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = events.DEP, tgrid=tsamp) %>%
    filter(time!=0)     # The code can produce time-dependent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  outdf.DEP   = cbind.data.frame (Time       = out.DEP$time, 
                                  AUC_CA     = out.DEP$AUC_CA,
                                  AUC_CAM    = out.DEP$AUC_CAM) 
  
  out.DEHP <- 
    mod.DEHP  %>%
    param(pars.DEHP) %>%
    Req(AUC_CA,AUC_CAM)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = events.DEHP, tgrid=tsamp) %>%
    filter(time!=0)     # The code can produce time-DEHPendent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  outdf.DEHP   = cbind.data.frame (Time       = out.DEHP$time, 
                                  AUC_CA      = out.DEHP$AUC_CA,
                                  AUC_CAM     = out.DEHP$AUC_CAM) 
   
   return (list("outdf.DEP"   = outdf.DEP,
                "outdf.DEHP"  = outdf.DEHP))
  
}

############## DINCH 
Pred.DINCH <- function (pars.DINCH){
  
  
  which_sig.DINCH        <- grep("sig", DINCH.names)
  
  pars.DINCH  <- lapply(pars.DINCH [-which_sig.DINCH],exp)
  
  
  ## Repeat dose exposure scenario: 
  
  BW                = 82.3                                ## human body weight
  
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
  
  PDOSEoral.DINCH   = 1                                ## ug/kg/d
  DOSEoral.DINCH   = PDOSEoral.DINCH * BW                                 ## ug/kg/d
  
  PDOSEother.DINCH   = 1                                ## ug/kg/d
  DOSEother.DINCH   = PDOSEother.DINCH * BW                                 ## ug/kg/d

 
  ############ DOSE
  
  ex.oral.DINCH        <- ev(ID=1, amt= DOSEoral.DINCH, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.other.DINCH       <- ev(ID=1, amt= DOSEother.DINCH, ii=tinterval, addl=TDoses-1, cmt="AC_DINCH", tinf = 24, replicate = FALSE)
  events.DINCH         <- ex.oral.DINCH + ex.other.DINCH
  #events.DINCH           <- ex.oral.DINCH
  
  
  ## set up the exposure time
  tsamp           = tgrid(0,tinterval*(TDoses-1)+24*7,1)          ## 2 weeks and simulated for 24 * 7 hours after dosing at one-hour resolution
  
  
  ## Get a prediction
  out.DINCH <- 
    mod.DINCH  %>%
    param(pars.DINCH) %>%
    Req(AUC_CA,AUC_CAM)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = events.DINCH, tgrid=tsamp) %>%
    filter(time!=0)     # The code can produce time-DINCHendent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  outdf.DINCH   = cbind.data.frame (Time      = out.DINCH$time, 
                                    AUC_CA      = out.DINCH$AUC_CA,
                                    AUC_CAM     = out.DINCH$AUC_CAM) 
  
  return (list("outdf.DINCH" = outdf.DINCH))
  
}

########################## DEHA
Pred.DEHA <- function (pars.DEHA){

  which_sig.DEHA         <- grep("sig", DEHA.names)
  
  pars.DEHA   <- lapply(pars.DEHA  [-which_sig.DEHA],exp)
  
  
  ## Repeat dose exposure scenario: 
  
  BW                = 82.3                                ## human body weight
  
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
  
  PDOSEoral.DEHA    = 1                                ## ug/kg/d
  DOSEoral.DEHA    = PDOSEoral.DEHA  * BW                                ## ug/kg/d
  
  PDOSEother.DEHA    = 1                                ## ug/kg/d
  DOSEother.DEHA    = PDOSEother.DEHA  * BW                                ## ug/kg/d
  
 
  
  ex.oral.DEHA         <- ev(ID=1, amt= DOSEoral.DEHA, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.other.DEHA        <- ev(ID=1, amt= DOSEother.DEHA, ii=tinterval, addl=TDoses-1, cmt="AC_DEHA", tinf = 24, replicate = FALSE)
  events.DEHA          <- ex.oral.DEHA + ex.other.DEHA
  #events.DEHA           <- ex.oral.DEHA
  
  ## set up the exposure time
  tsamp           = tgrid(0,tinterval*(TDoses-1)+24*7,1)          ## 2 weeks and simulated for 24 * 7 hours after dosing at one-hour resolution
  
  out.DEHA <- 
    mod.DEHA  %>%
    param(pars.DEHA) %>%
    Req(AUC_CA,AUC_CAM)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = events.DEHA, tgrid=tsamp) %>%
    filter(time!=0)     # The code can produce time-DEHAendent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  outdf.DEHA   = cbind.data.frame (Time       = out.DEHA$time, 
                                   AUC_CA     = out.DEHA$AUC_CA,
                                   AUC_CAM    = out.DEHA$AUC_CAM) 
  
  return (list("outdf.DEHA"  = outdf.DEHA))
  
}



############## DPHP 
Pred.DPHP <- function (pars.DPHP){
  
  
  which_sig.DPHP        <- grep("sig", DPHP.names)
  
  pars.DPHP  <- lapply(pars.DPHP [-which_sig.DPHP],exp)
  
  
  ## Repeat dose exposure scenario: 
  
  BW                = 82.3                                ## human body weight
  
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
  
  PDOSEoral.DPHP   = 1                                ## ug/kg/d
  DOSEoral.DPHP   = PDOSEoral.DPHP * BW                                 ## ug/kg/d
  
  PDOSEother.DPHP   = 1                                ## ug/kg/d
  DOSEother.DPHP   = PDOSEother.DPHP * BW                                 ## ug/kg/d
  
  
  ############ DOSE
  
  ex.oral.DPHP        <- ev(ID=1, amt= DOSEoral.DPHP, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.other.DPHP       <- ev(ID=1, amt= DOSEother.DPHP, ii=tinterval, addl=TDoses-1, cmt="AC_DPHP", tinf = 24, replicate = FALSE)
  events.DPHP         <- ex.oral.DPHP + ex.other.DPHP
  #events.DPHP           <- ex.oral.DPHP
  
  
  ## set up the exposure time
  tsamp           = tgrid(0,tinterval*(TDoses-1)+24*7,1)          ## 2 weeks and simulated for 24 * 7 hours after dosing at one-hour resolution
  
  
  ## Get a prediction
  out.DPHP <- 
    mod.DPHP  %>%
    param(pars.DPHP) %>%
    Req(AUC_CA,AUC_CAM)%>%
    update(atol = 1E-5,rtol=1E-5,maxsteps=2000) %>%
    mrgsim_d(data = events.DPHP, tgrid=tsamp) %>%
    filter(time!=0)     # The code can produce time-DPHPendent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  outdf.DPHP   = cbind.data.frame (Time      = out.DPHP$time, 
                                    AUC_CA      = out.DPHP$AUC_CA,
                                    AUC_CAM     = out.DPHP$AUC_CAM) 
  
  return (list("outdf.DPHP" = outdf.DPHP))
  
}


####################### END of MODEL ####################

pars.DEP   = DEP.MCMC$bestpar
pars.DEHA  = DEHA.MCMC$bestpar
pars.DEHP  = DEHP.MCMC$bestpar
pars.DINCH = DINCH.MCMC$bestpar
pars.DPHP  = DPHP.MCMC$bestpar

R = Pred(pars.DEP, pars.DEHP)
R.DINCH = Pred.DINCH(pars.DINCH)
R.DEHA = Pred.DEHA(pars.DEHA)
R.DPHP = Pred.DPHP(pars.DPHP)

################## Create the matrix for normalized sensitivity coefficient data ##################
NSC_DEP    = matrix(nrow=6,ncol=2)
NSC_DEHA   = matrix(nrow=13,ncol=2)
NSC_DEHP   = matrix(nrow=6,ncol=2)
NSC_DINCH  = matrix(nrow=18,ncol=2)
NSC_DPHP   = matrix(nrow=14,ncol=2)


for (i in 1:6) {
  
  pars.DEP.new        <- log(c(exp(pars.DEP[i])*1.01,exp(pars.DEP[-i])))
  pars.DEHP.new       <- log(c(exp(pars.DEHP[i])*1.01,exp(pars.DEHP[-i])))
  # Each cycle, generate a new value of parameter i, and delete parameter i, so that you can proceed to the next parameter i+1
  
  Rnew                <- Pred(pars.DEP.new, pars.DEHP.new)
  delta.P.DEP         <- 100 #exp(pars.DEP[i])/(exp(pars.DEP[i])) * 1#00
  delta.P.DEHP        <- 100 #exp(pars.DEHP[i])/(exp(pars.DEHP[i])) * 1#00
  
  
  ## Estimated the AUC
  ## DEP
  DEP.AUC.CA.new      =  Rnew$outdf.DEP %>% filter (Time == 14) %>% select (AUC_CA)
  DEP.AUC.CA.ori      =  R$outdf.DEP    %>% filter (Time == 14) %>% select (AUC_CA)
  DEP.AUC.CAM.new     =  Rnew$outdf.DEP %>% filter (Time == 14) %>% select (AUC_CAM)
  DEP.AUC.CAM.ori     =  R$outdf.DEP    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DEP    =  DEP.AUC.CA.new  - DEP.AUC.CA.ori
  delta.AUC.CAM.DEP   =  DEP.AUC.CAM.new - DEP.AUC.CAM.ori

  ## DEHP
  DEHP.AUC.CA.new      =  Rnew$outdf.DEHP %>% filter (Time == 14) %>% select (AUC_CA)
  DEHP.AUC.CA.ori      =  R$outdf.DEHP    %>% filter (Time == 14) %>% select (AUC_CA)
  DEHP.AUC.CAM.new     =  Rnew$outdf.DEHP %>% filter (Time == 14) %>% select (AUC_CAM)
  DEHP.AUC.CAM.ori     =  R$outdf.DEHP    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DEHP    =  DEHP.AUC.CA.new  - DEHP.AUC.CA.ori
  delta.AUC.CAM.DEHP   =  DEHP.AUC.CAM.new - DEHP.AUC.CAM.ori
  
 
  NSC_DEP    [i, 1]   <- as.numeric((delta.AUC.CA.DEP/DEP.AUC.CA.ori)     * delta.P.DEP)
  NSC_DEP    [i, 2]   <- as.numeric((delta.AUC.CAM.DEP/DEP.AUC.CAM.ori)   * delta.P.DEP)
  NSC_DEHP   [i, 1]   <- as.numeric((delta.AUC.CA.DEHP/DEHP.AUC.CA.ori)   * delta.P.DEHP)
  NSC_DEHP   [i, 2]   <- as.numeric((delta.AUC.CAM.DEHP/DEHP.AUC.CAM.ori) * delta.P.DEHP)
  
  
  cat("iteration = ", i , "\n") 

}

print(pars.DEHA.new)

for (i in 1:13) {
  #DEHA.old       <- log((exp(pars.DEHA[i])*1))  
 # DEHA.new       <- log((exp(pars.DEHA[i])*1.01))  

  #pars.DEHA.new  [i]   <- DEHA.new
  pars.DEHA.new       <- log(c(exp(pars.DEHA[i])*1.01,exp(pars.DEHA[-i]))) 

  Rnew.DEHA           <- Pred.DEHA(pars.DEHA.new)
  delta.P.DEHA        <- exp(pars.DEHA[i])/(exp(pars.DEHA[i])) * 100 # Basically, the ratio is 100 for each parameter.

  ## Estimated the AUC
  
  ## DEHA
  DEHA.AUC.CA.new      =  Rnew.DEHA$outdf.DEHA %>% filter (Time == 14) %>% select (AUC_CA)
  DEHA.AUC.CA.ori      =  R.DEHA$outdf.DEHA %>% filter (Time == 14) %>% select (AUC_CA)
  DEHA.AUC.CAM.new     =  Rnew.DEHA$outdf.DEHA %>% filter (Time == 14) %>% select (AUC_CAM)
  DEHA.AUC.CAM.ori     =  R.DEHA$outdf.DEHA %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DEHA    =  DEHA.AUC.CA.new - DEHA.AUC.CA.ori
  delta.AUC.CAM.DEHA   =  DEHA.AUC.CAM.new - DEHA.AUC.CAM.ori
  
  NSC_DEHA   [i, 1]   <- as.numeric((delta.AUC.CA.DEHA/DEHA.AUC.CA.ori) * delta.P.DEHA)
  NSC_DEHA   [i, 2]   <- as.numeric((delta.AUC.CAM.DEHA/DEHA.AUC.CAM.ori) * delta.P.DEHA)
 
 # pars.DEHA.new  [i]   <- DEHA.old 
  
  cat("iteration = ", i , "\n") 
}


for (i in 1:18) {
  pars.DINCH.new      <- log(c(exp(pars.DINCH[i])*1.01,exp(pars.DINCH[-i])))
 
  Rnew                <- Pred.DINCH(pars.DINCH.new)
  delta.P.DINCH       <- exp(pars.DINCH[i])/(exp(pars.DINCH[i])) * 100
  
  
  ## Estimated the AUC
  ## DINCH
  DINCH.AUC.CA.new      =  Rnew$outdf.DINCH %>% filter (Time == 14) %>% select (AUC_CA)
  DINCH.AUC.CA.ori      =  R.DINCH$outdf.DINCH    %>% filter (Time == 14) %>% select (AUC_CA)
  DINCH.AUC.CAM.new     =  Rnew$outdf.DINCH %>% filter (Time == 14) %>% select (AUC_CAM)
  DINCH.AUC.CAM.ori     =  R.DINCH$outdf.DINCH    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DINCH    =  DINCH.AUC.CA.new - DINCH.AUC.CA.ori
  delta.AUC.CAM.DINCH   =  DINCH.AUC.CAM.new - DINCH.AUC.CAM.ori
  
  NSC_DINCH  [i, 1]   <- as.numeric((delta.AUC.CA.DINCH/DINCH.AUC.CA.ori)   * delta.P.DINCH)
  NSC_DINCH  [i, 2]   <- as.numeric((delta.AUC.CAM.DINCH/DINCH.AUC.CAM.ori) * delta.P.DINCH)
  
  cat("iteration = ", i , "\n") 
  
}

###### DPHP #####################

for (i in 1:14) {
  pars.DPHP.new      <- log(c(exp(pars.DPHP[i])*1.01,exp(pars.DPHP[-i])))
  
  Rnew                <- Pred.DPHP(pars.DPHP.new)
  delta.P.DPHP       <- exp(pars.DPHP[i])/(exp(pars.DPHP[i])) * 100
  
  
  ## Estimated the AUC
  ## DPHP
  DPHP.AUC.CA.new      =  Rnew$outdf.DPHP %>% filter (Time == 14) %>% select (AUC_CA)
  DPHP.AUC.CA.ori      =  R.DPHP$outdf.DPHP    %>% filter (Time == 14) %>% select (AUC_CA)
  DPHP.AUC.CAM.new     =  Rnew$outdf.DPHP %>% filter (Time == 14) %>% select (AUC_CAM)
  DPHP.AUC.CAM.ori     =  R.DPHP$outdf.DPHP    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DPHP    =  DPHP.AUC.CA.new  - DPHP.AUC.CA.ori
  delta.AUC.CAM.DPHP   =  DPHP.AUC.CAM.new - DPHP.AUC.CAM.ori
  
  NSC_DPHP  [i, 1]   <- as.numeric((delta.AUC.CA.DPHP/DPHP.AUC.CA.ori)   * delta.P.DPHP)
  NSC_DPHP  [i, 2]   <- as.numeric((delta.AUC.CAM.DPHP/DPHP.AUC.CAM.ori) * delta.P.DPHP)
  
  
  cat("iteration = ", i , "\n") 
}


colnames (NSC_DEP)   = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DEP)   = DEP.names[1:6]
colnames (NSC_DEHA)  = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DEHA)  = DEHA.names[1:13]
colnames (NSC_DEHP)  = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DEHP)  = DEHP.names[1:6]
colnames (NSC_DINCH) = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DINCH) = DINCH.names[1:18]
colnames (NSC_DPHP)  = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DPHP)  = DPHP.names[1:14]


##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################

Circle.plot <- function (melt.data){ # melt.data is an argument of Circle.plot function.

# Set a number of 'empty bar' to add at the end of each group
  empty_bar= 3
  to_add = data.frame(matrix(NA, empty_bar*nlevels(as.factor(melt.data$group)), ncol(melt.data)) )
  colnames(to_add) = colnames(melt.data)
  to_add$group=rep(levels(as.factor(melt.data$group)), each=empty_bar)
  melt.data=rbind(melt.data, to_add)
  melt.data=melt.data %>% arrange(group)
  melt.data$id=seq(1, nrow(melt.data)) # id is the number of rows. In total, there were 6+6+13+6+18 rows.

# Get the name and the y position of each label
 label_data=melt.data
 number_of_bar=nrow(label_data) # in total, there were 27 rows, 27 parameters * 4 
 angle<- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
 label_data$hjust<-ifelse( angle < -90, 1, 0)
 label_data$angle<-ifelse(angle < -90, angle+180, angle)

 # prepare a data frame for base lines
 base_data=melt.data %>%
 group_by(group) %>%
 summarize(start=min(id), end=max(id) - empty_bar) %>%
 rowwise() %>%
 mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
 grid_data = base_data
 grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
 grid_data$start = grid_data$start - 1
 grid_data = grid_data[-1,]

# Make the plot
 windowsFonts(Times=windowsFont("Times New Roman"))
 
 p.cir.plot <- 
 ggplot(melt.data, aes(x=as.factor(id), y = abs(value*100), fill = group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
 geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=1) +
 scale_fill_brewer(palette = "Dark2") +
   

# Add a val=80/60/40/20 lines. I do it at the beginning to make sure barplots are OVER it.
 geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "#6E6E6E", alpha=1, size=0.8 , inherit.aes = FALSE ) +
 geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "#6E6E6E", alpha=1, size=0.8 , inherit.aes = FALSE ) +
 geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "#6E6E6E", alpha=1, size=0.8 , inherit.aes = FALSE ) +
 geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "#6E6E6E", alpha=1, size=0.8 , inherit.aes = FALSE ) +
  
# Add text showing the value of each 80/60/40/20 lines
 annotate("text", x = rep(max(melt.data$id),4), y = c(20, 40, 60, 80), label = c("20%", "40%", "60%", "80%") , 
          color="black", size=10 , angle=0, fontface="bold", hjust=1) +
 #geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.9) +
 ylim(-100,120) +
 theme_minimal() +
 theme(
    legend.position         = "none",
    text                    = element_text (family = "Times"),
    panel.background        = element_blank (),
    plot.background         = element_blank (),
    axis.text               = element_blank(),
    axis.title              = element_blank(),
    panel.grid              = element_blank(),
    plot.margin             = unit(rep(-1,4), "cm")
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=abs(value*100)+4, label=par, hjust=hjust), 
            color="black", fontface="bold.italic",
            alpha = 1, size = 16, angle= label_data$angle, inherit.aes = FALSE,parse = TRUE) +
  
   
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -4, xend = end, yend = -4), colour = "black", alpha=1, size=1.2 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -12, label=group), hjust=c(1,1,0.8,0,0), colour = "black", alpha=0.8, size= 11, fontface="bold", inherit.aes = FALSE)
  
return (p.cir.plot)
}

###################### Fig. 8a; AUC of plasma_parent ####################
melt.DEP.CA        = melt(NSC_DEP[,1])
melt.DEP.CA$group  = c("DEP")
melt.DEP.CA$par    = rownames(melt.DEP.CA)

melt.DEHA.CA        = melt(NSC_DEHA[,1])
melt.DEHA.CA$group  = c("DEHA")
melt.DEHA.CA$par    = rownames(melt.DEHA.CA)

melt.DEHP.CA        = melt(NSC_DEHP[,1])
melt.DEHP.CA$group  = c("DEHP")
melt.DEHP.CA$par    = rownames(melt.DEHP.CA)

melt.DINCH.CA        = melt(NSC_DINCH[,1])
melt.DINCH.CA$group  = c("DINCH")
melt.DINCH.CA$par    = rownames(melt.DINCH.CA)

melt.DPHP.CA        = melt(NSC_DPHP[,1])
melt.DPHP.CA$group  = c("DPHP")
melt.DPHP.CA$par    = rownames(melt.DPHP.CA)


melt.data.CA         = rbind (melt.DEP.CA, melt.DEHA.CA, melt.DEHP.CA, melt.DINCH.CA, melt.DPHP.CA)


library(tidyverse)
melt.data.CA  = melt.data.CA  %>%
  arrange(value) %>% 
  mutate(par=gsub("Vp", "italic(V[p])", par),
         par=gsub("Vm", "italic(V[m])", par),
         par=gsub("ka", "italic(k[a])", par),
         par=gsub("lambda_p", "italic(lambda[p])", par),
         par=gsub("lambda_zm", "italic(lambda[zm])", par),
         par=gsub("lambda_u", "italic(lambda[u])", par))


p                = Circle.plot (melt.data.CA%>%filter(abs(value)>0.01))
p


####### Save plot #######
ggsave("Sens_parentB_newAbcd.tiff",scale = 1,
       plot = p,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot",
       width = 40, height = 40, units = "cm",dpi=320)


#dev.off()


melt.DEP.CAM        = melt(NSC_DEP[,2])
melt.DEP.CAM$group  = c("MEP")
melt.DEP.CAM$par    = rownames(melt.DEP.CAM)

melt.DEHA.CAM        = melt(NSC_DEHA[,2])
melt.DEHA.CAM$group  = c("MEHA")
melt.DEHA.CAM$par    = rownames(melt.DEHA.CAM)

melt.DEHP.CAM        = melt(NSC_DEHP[,2])
melt.DEHP.CAM$group  = c("MEHP")
melt.DEHP.CAM$par    = rownames(melt.DEHP.CAM)

melt.DINCH.CAM        = melt(NSC_DINCH[,2])
melt.DINCH.CAM$group  = c("MINCH")
melt.DINCH.CAM$par    = rownames(melt.DINCH.CAM)

melt.DPHP.CAM        = melt(NSC_DPHP[,2])
melt.DPHP.CAM$group  = c("MPHP")
melt.DPHP.CAM$par    = rownames(melt.DPHP.CAM)


melt.data.CAM         = rbind (melt.DEP.CAM, melt.DEHA.CAM, melt.DEHP.CAM, melt.DINCH.CAM, melt.DPHP.CAM)

melt.data.CAM  = melt.data.CAM  %>%
  arrange(value) %>% 
  mutate(par=gsub("Vp", "italic(V[p])", par),
         par=gsub("Vm", "italic(V[m])", par),
         par=gsub("ka", "italic(k[a])", par),
         par=gsub("lambda_p", "italic(lambda[p])", par),
         par=gsub("lambda_zm", "italic(lambda[zm])", par),
         par=gsub("lambda_u", "italic(lambda[u])", par))



p2                = Circle.plot (melt.data.CAM%>%filter(abs(value)>0.01))
p2

ggsave("Sens_metaB_newAb.tiff",scale = 1,
       plot = p2,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot",
       width = 40, height = 40, units = "cm",dpi=320)

