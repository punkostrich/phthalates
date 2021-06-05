## loading R packages
library(magrittr)   # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)      # Needed for the pipe %>% operator
library(mrgsolve)   # Needed to run the main PBPK code
library(reshape)    # melt function to reshape the table
library(ggplot2)    # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(grid)       # for plotting the figure
library(lattice)    # for plotting the figure

dev.off()
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


product.DEP.names         <- readRDS(file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/product.DEP.RDS')
product.DEHA.names        <- readRDS(file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/product.DEHA.RDS')
product.DEHP.names        <- readRDS(file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/product.DEHP.RDS')
product.DINCH.names       <- readRDS(file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/product.DINCH.RDS')
product.DPHP.names        <- readRDS(file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/product.DPHP.RDS')

## Sensitivity analysis 3 CHEMICAL DEHP
Pred.DEP <- function (product.pars.DEP){
  
  pars.DEP       =    DEP.MCMC$bestpar
  
  which_sig.DEP          <- grep("sig", DEP.names)
 
  ## Get out of log domain
  pars.DEP    <- lapply(pars.DEP  [-which_sig.DEP],exp) # pars.mouse is the entire parameter set; pars.mouse [-which_sig] means to keep parameters with "sig" only, and then do exp transformation, then reassign to pars.mouse

  ## Repeat dose exposure scenario: 
   
  BW                = 82.3                                ## human body weight
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
  
  PDOSEoral.DEP     = product.pars.DEP[2]                                ## ug/kg/d oral do## ug/kg/d
  DOSEoral.DEP      = PDOSEoral.DEP   * BW                                ## ug/kg/d
 
  PDOSEother.DEP    = product.pars.DEP[1] + product.pars.DEP[3]                                ## ug/kg/d other exposure dose
  DOSEother.DEP     = PDOSEother.DEP   * BW                                ## ug/kg/d

  ############ DOSE
  ex.oral.DEP          <- ev(ID=1, amt= DOSEoral.DEP, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.other.DEP         <- ev(ID=1, amt= DOSEother.DEP, ii=tinterval, addl=TDoses-1, cmt="AC", tinf = 24, replicate = FALSE)
  events.DEP           <- ex.oral.DEP + ex.other.DEP
  
 
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
  
 
   return (list("outdf.DEP"   = outdf.DEP))
  
}

################# DEHP ########################
Pred.DEHP <- function (product.pars.DEHP){
  
  pars.DEHP       =    DEHP.MCMC$bestpar
  
  which_sig.DEHP          <- grep("sig", DEHP.names)

  ## Get out of log domain
  pars.DEHP    <- lapply(pars.DEHP  [-which_sig.DEHP],exp) # pars.mouse is the entire parameter set; pars.mouse [-which_sig] means to keep parameters with "sig" only, and then do exp transformation, then reassign to pars.mouse
  
  ## Repeat dose exposure scenario: 
  
  BW                = 82.3                                ## human body weight
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
  
  PDOSEoral.DEHP     = product.pars.DEHP[2]                                ## ug/kg/d oral do## ug/kg/d
  PDOSEother.DEHP    = product.pars.DEHP[1] + product.pars.DEHP[3] 
  
  DOSEoral.DEHP      = PDOSEoral.DEHP   * BW                                ## ug/kg/d
  DOSEother.DEHP     = PDOSEother.DEHP   * BW                                ## ug/kg/d
  
  ############ DOSE
  ex.oral.DEHP          <- ev(ID=1, amt= DOSEoral.DEHP, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.other.DEHP         <- ev(ID=1, amt= DOSEother.DEHP, ii=tinterval, addl=TDoses-1, cmt="AC", tinf = 24, replicate = FALSE)
  events.DEHP           <- ex.oral.DEHP + ex.other.DEHP
  

  ## set up the exposure time
  tsamp           = tgrid(0,tinterval*(TDoses-1)+24*7,1)          ## 2 weeks and simulated for 24 * 7 hours after dosing at one-hour resolution

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
  
  return (list("outdf.DEHP"  = outdf.DEHP))
  
}




############## DINCH 
Pred.DINCH <- function (product.pars.DINCH){
  
  pars.DINCH       =    DINCH.MCMC$bestpar
  which_sig.DINCH        <- grep("sig", DINCH.names)
  
  pars.DINCH  <- lapply(pars.DINCH [-which_sig.DINCH],exp)
  
  
  ## Repeat dose exposure scenario: 
  
  BW                = 82.3                                ## human body weight
  
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
 
  PDOSEoral.DINCH    = product.pars.DINCH[2]                                ## ug/kg/d oral do## ug/kg/d
  PDOSEother.DINCH   = product.pars.DINCH[1] + product.pars.DINCH[3] 
  
  DOSEoral.DINCH     = PDOSEoral.DINCH * BW                                 ## ug/kg/d
  DOSEother.DINCH    = PDOSEother.DINCH * BW                                 ## ug/kg/d

 
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
Pred.DEHA <- function (product.pars.DEHA){

  pars.DEHA       =    DEHA.MCMC$bestpar
  which_sig.DEHA         <- grep("sig", DEHA.names)
  
  pars.DEHA   <- lapply(pars.DEHA  [-which_sig.DEHA],exp)
  
  
  ## Repeat dose exposure scenario: 
  
  BW                = 82.3                                ## human body weight
  
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
  
  PDOSEoral.DEHA    = product.pars.DEHA[2]                                ## ug/kg/d oral do## ug/kg/d
  PDOSEother.DEHA   = product.pars.DEHA[1] + product.pars.DEHA[3] 
  
  DOSEoral.DEHA    = PDOSEoral.DEHA  * BW                                ## ug/kg/d
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
Pred.DPHP <- function (product.pars.DPHP){
  
  pars.DPHP      =    DPHP.MCMC$bestpar
  which_sig.DPHP        <- grep("sig", DPHP.names)
  
  pars.DPHP  <- lapply(pars.DPHP [-which_sig.DPHP],exp)
  
  
  ## Repeat dose exposure scenario: 
  
  BW                = 82.3                                ## human body weight
  
  tinterval         = 24                                  ## Time interval
  TDoses            = 14                                  ## The number of dosing for two weeks
  
  PDOSEoral.DPHP    = product.pars.DPHP[2]                                ## ug/kg/d oral do## ug/kg/d
  PDOSEother.DPHP   = product.pars.DPHP[1] + product.pars.DPHP[3] 
                              
  DOSEoral.DPHP     = PDOSEoral.DPHP * BW                                 ## ug/kg/d
  DOSEother.DPHP    = PDOSEother.DPHP * BW                                 ## ug/kg/d
  
  
  ############ DOSE
  
  ex.oral.DPHP        <- ev(ID=1, amt= DOSEoral.DPHP, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex.other.DPHP       <- ev(ID=1, amt= DOSEother.DPHP, ii=tinterval, addl=TDoses-1, cmt="AC_DPHP", tinf = 24, replicate = FALSE)
  events.DPHP         <- ex.oral.DPHP + ex.other.DPHP
  
  
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


####################### END of PK MODEL #############################
#############                                   #####################
###################### START OF EXPOSURE MODEL ######################
#### DEP
Pred.Product.DEP <- function (pars){
  
 # perc   = pars[1]
  C_pcp  = pars[1]
  
  sn = 1000000
  
  ##### percentage in plastic
  life = 10     # plastic service life
  dens = 2   #kg/m3
  dens_DEP <- 1.12 * 1E3 #kg/m3
  Cpcp_ori <- 0.002
  
  ######## household size ########
  household = 2.65
  
  A   = 1
  y0  = min(19113, 19113 * Cpcp_ori / dens_DEP / (Cpcp_ori/dens_DEP + (1-Cpcp_ori)/1E3) * 3.7)
  
  
  ## convective mass transfer coefficient log-normal
  h    =  3.6 
  ## indoor ventilation rate Q normal
  Q    = 221 
  ## inhalation rate normal
  InhR = 15 
  ## boday weight normal
  BW   = 80.8         ## kg
  ## tsp log-normal
  TSP  = 20             
  ## exposure duration percentage assume 90% indoor
  ED   = 0.9  ## -
  ## particle-air partition coefficient normal
  Kp   = 4.6E-5 
  ## Calculate DEHP conc in air 
  Cair = h * y0 * A/(h * A + (1 + Kp * TSP) * Q)    ## ug/m3
  ## inhalation exposure
  IE   = (Cair * InhR * ED + Cair * Kp * TSP * InhR * ED)/BW   ## ug/day/kg
  
  
  ######################################## 2 ####################################
  ################################ Exposure scenario #Oral ######################
  
  #dust oral
  IngR  =  50  ## mg/d
  Kdust <- Kp/8.32                ## m3/ug
  OE = (Cair * Kdust * IngR)/BW  
  
  
  ############################## 3 ###################################
  ##################### Exposure scenario #Dermal ####################
  kpg = 1.56
  SA =  2.1 
  Fr = 0.24
  
  cosmatic = 5.6
  
  SE <- (Cair * kpg * SA * Fr * ED)/BW + (C_pcp * cosmatic * 0.005 *0.5)/BW  ## ug/day/kg
  
  
  exposure = IE + OE + SE
  
  return(c( IE, OE, SE))
  
}

####### DEHA
Pred.Product.DEHA <- function (pars){
  
  volume = pars[1]
  perc   = pars[2]
  Cp_ori = pars[3]
  
  sn = 1000000
  
  ##### percentage in plastic
  life = 10     # plastic service life
  dens = 2   #kg/m3
  dens_DEHA   <- 0.922 * 1E3 #kg/m3

  
  ######## household size ########
  household = 2.65
  
  A_p = volume * life / perc / dens / 328 /1000000
  A   = volume * life / perc / dens / 328 /1000000 * household 
  y0  = min(10.4,10.4 * perc / dens_DEHA / (perc/dens_DEHA + (1-perc)/1.5E3) * 3.7) 
  
  
  ## convective mass transfer coefficient log-normal
  h =  3.6 
  ## indoor ventilation rate Q normal
  Q = 221 
  ## inhalation rate normal
  InhR = 15 
  ## boday weight normal
  BW = 80.8         ## kg
  ## tsp log-normal
  TSP = 20             
  ## exposure duration percentage assume 90% indoor
  ED = 0.9  ## -
  ## particle-air partition coefficient normal
  Kp = 0.02 
  ## Calculate DEHP conc in air 
  Cair = h * y0 * A/(h * A + (1 + Kp * TSP) * Q)    ## ug/m3
  ## inhalation exposure
  IE = (Cair * InhR * ED + Cair * Kp * TSP * InhR * ED)/BW   ## ug/day/kg
  
  
  ######################################## 2 ####################################
  ################################ Exposure scenario #Oral ######################
  #diet oral
  ##kpf of DEHP
  kpf = 1410
  
  ## concentration in food contacting material Cp (%) log-normal
  #Cp_ori = 0.0005 
  Cp    = Cp_ori * 1E6   ## ug/g
  ## meat consumption (g/kg-day) normal
  meat  = 1.8 
  ## cheese consmption normal
  dairy = 16 
  ## fat consumption normal
  fat   = 1.2   
  Diet  = Cp/kpf * (meat*BW + dairy + fat*BW)      ##ug/day
  
  #dust oral
  IngR  = 50  ## mg/d
  Kdust <- Kp/8.32  
  OE = (Cair * Kdust * IngR)/BW  + Diet/BW ## ug/day/kg
  
  
  ############################## 3 ###################################
  ##################### Exposure scenario #Dermal ####################
  kpg = 1.77
  SA =  2.1 
  Fr  = 0.24 
  SE = (Cair * kpg * SA * Fr * ED)/BW  ## ug/day/kg
  
  exposure = IE + OE + SE
  
  return(c( IE, OE, SE))
  
}


############### DEHP 
Pred.Product.DEHP <- function (pars){
  
  volume = pars[1]
  perc   = pars[2]
  Cp_ori = pars[3]
  
  sn = 1000000
  
  ##### percentage in plastic
  
  #perc = 0.3 
  life = 10     # plastic service life
  dens = 2   #kg/m3
  dens_DEHP   <- 0.99 * 1E3 #kg/m3
  ######## household size ########
  household = 2.65
  
  A_p = volume * life / perc / dens / 328 /1000000
  A = volume * life / perc / dens / 328 /1000000 * household 
  y0  = min(2.5, 2.5 * perc / dens_DEHP / (perc/dens_DEHP + (1-perc)/1.5E3) * 3.7) ## ug/m3,
  
  
  ## convective mass transfer coefficient log-normal
  h =  3.6 
  ## indoor ventilation rate Q normal
  Q = 221 
  ## inhalation rate normal
  InhR = 15 
  ## boday weight normal
  BW = 80.8         ## kg
  ## tsp log-normal
  TSP = 20             
  ## exposure duration percentage assume 90% indoor
  ED = 0.9  ## -
  ## particle-air partition coefficient normal
  Kp = 0.06 
  ## Calculate DEHP conc in air 
  Cair = h * y0 * A/(h * A + (1 + Kp * TSP) * Q)    ## ug/m3
  ## inhalation exposure
  IE = (Cair * InhR * ED + Cair * Kp * TSP * InhR * ED)/BW   ## ug/day/kg
  
  
  ######################################## 2 ####################################
  ################################ Exposure scenario #Oral ######################
  #diet oral
  ##kpf of DEHP
  kpf = 5126
  
  ## concentration in food contacting material Cp (%) log-normal
  #Cp_ori = 0.0005 
  Cp = Cp_ori * 1E6   ## ug/g
  ## meat consumption (g/kg-day) normal
  meat = 1.8 
  ## cheese consmption normal
  dairy = 16 
  ## fat consumption normal
  fat = 1.2   
  Diet = Cp/kpf * (meat*BW + dairy + fat*BW)      ##ug/day
  
  #dust oral
  IngR  = 50  ## mg/d
  Kdust <- Kp/8.32  
  OE = (Cair * Kdust * IngR)/BW  + Diet/BW ## ug/day/kg
  
  
  ############################## 3 ###################################
  ##################### Exposure scenario #Dermal ####################
  kpg = 4.68
  SA = 2.1 
  Fr  = 0.24 
  SE = (Cair * kpg * SA * Fr * ED)/BW  ## ug/day/kg
  
  exposure = IE + OE + SE
  
  return(c( IE, OE, SE))
  
}


################# DINCH ############
Pred.Product.DINCH <- function (pars){
  
  volume = pars[1]
  perc   = pars[2]
  Cp_ori = pars[3]
  dens_DINCH   <- 0.954 * 1E3 #kg/m3
  
  sn = 1000000
  
  ##### percentage in plastic
  life = 10     # plastic service life
  dens = 2   #kg/m3
  
  ######## household size ########
  household = 2.65
  
  A_p = volume * life / perc / dens / 328 /1000000
  A   = volume * life / perc / dens / 328 /1000000 * household 
  y0  = min(1.6, 1.6 * perc / dens_DINCH / (perc/dens_DINCH + (1-perc)/1.5E3) * 3.7) ## ug/m3,
  
  
  ## convective mass transfer coefficient log-normal
  h    =  3.6 
  ## indoor ventilation rate Q normal
  Q    = 221 
  ## inhalation rate normal
  InhR = 15 
  ## boday weight normal
  BW   = 80.8         ## kg
  ## tsp log-normal
  TSP  = 20             
  ## exposure duration percentage assume 90% indoor
  ED   = 0.9  ## -
  ## particle-air partition coefficient normal
  Kp   = 0.09 
  ## Calculate DEHP conc in air 
  Cair = h * y0 * A/(h * A + (1 + Kp * TSP) * Q)    ## ug/m3
  ## inhalation exposure
  IE   = (Cair * InhR * ED + Cair * Kp * TSP * InhR * ED)/BW   ## ug/day/kg
  
  
  ######################################## 2 ####################################
  ################################ Exposure scenario #Oral ######################
  #diet oral
  ##kpf of DEHP
  kpf = 213981
  
  ## concentration in food contacting material Cp (%) log-normal
  #Cp_ori = 0.0005 
  Cp    = Cp_ori * 1E6   ## ug/g
  ## meat consumption (g/kg-day) normal
  meat  = 1.8 
  ## cheese consmption normal
  dairy = 16 
  ## fat consumption normal
  fat   = 1.2   
  Diet  = Cp/kpf * (meat*BW + dairy + fat*BW)      ##ug/day
  
  #dust oral
  IngR  =  50  ## mg/d
  Kdust <- Kp/8.32  
  OE = (Cair * Kdust * IngR)/BW  + Diet/BW ## ug/day/kg
  
  
  ############################## 3 ###################################
  ##################### Exposure scenario #Dermal ####################
  kpg = 11
  SA =  2.1 
  Fr  = 0.24 
  SE = (Cair * kpg * SA * Fr * ED)/BW  ## ug/day/kg
  
  exposure = IE + OE + SE
  
  return(c( IE, OE, SE))
  
}

############# DPHP
Pred.Product.DPHP <- function (pars){
  
  volume = pars[1]
  perc   = pars[2]
  Cp_ori = pars[3]
  
  sn = 1000000
  
  ##### percentage in plastic
  life = 10     # plastic service life
  dens = 2   #kg/m3
  dens_DPHP   <- 0.965 * 1E3 #kg/m3
  
  ######## household size ########
  household = 2.65
  
  A_p = volume * life / perc / dens / 328 /1000000
  A   = volume * life / perc / dens / 328 /1000000 * household 
  y0  = min(0.667, 0.667 * perc / dens_DPHP / (perc/dens_DPHP + (1-perc)/1.5E3) * 3.7) ## ug/m3,
  
  
  ## convective mass transfer coefficient log-normal
  h    =  3.6 
  ## indoor ventilation rate Q normal
  Q    = 221 
  ## inhalation rate normal
  InhR = 15 
  ## boday weight normal
  BW   = 80.8         ## kg
  ## tsp log-normal
  TSP  = 20             
  ## exposure duration percentage assume 90% indoor
  ED   = 0.9  ## -
  ## particle-air partition coefficient normal
  Kp   = 0.28
  ## Calculate DEHP conc in air 
  Cair = h * y0 * A/(h * A + (1 + Kp * TSP) * Q)    ## ug/m3
  ## inhalation exposure
  IE   = (Cair * InhR * ED + Cair * Kp * TSP * InhR * ED)/BW   ## ug/day/kg
  
  
  ######################################## 2 ####################################
  ################################ Exposure scenario #Oral ######################
  #diet oral
  ##kpf of DEHP
  kpf = 793917
  
  ## concentration in food contacting material Cp (%) log-normal
  #Cp_ori = 0.0005 
  Cp    = Cp_ori * 1E6   ## ug/g
  ## meat consumption (g/kg-day) normal
  meat  = 1.8 
  ## cheese consmption normal
  dairy = 16 
  ## fat consumption normal
  fat   = 1.2   
  Diet  = Cp/kpf * (meat*BW + dairy + fat*BW)      ##ug/day
  
  #dust oral
  IngR  =  50  ## mg/d
  Kdust <- Kp/8.32  
  OE = (Cair * Kdust * IngR)/BW  + Diet/BW ## ug/day/kg
  
  
  ############################## 3 ###################################
  ##################### Exposure scenario #Dermal ####################
  kpg = 28.9
  SA =  2.1 
  Fr  = 0.24 
  SE = (Cair * kpg * SA * Fr * ED)/BW  ## ug/day/kg
  
  exposure = IE + OE + SE
  
  return(c( IE, OE, SE))
  
}

######################## END OF EXPOSURE MODEL #######################

pars.DEP   = DEP.MCMC$bestpar
pars.DEHA  = DEHA.MCMC$bestpar
pars.DEHP  = DEHP.MCMC$bestpar
pars.DINCH = DINCH.MCMC$bestpar
pars.DPHP  = DPHP.MCMC$bestpar

product.DEP <- c(
  #PCT     =   
   # 0.1,                # percentage in coc
  #C_pcp   =   
    0.01* 1E6                  ## concentration in pcpp 
  )

product.DEHA <- c(
  #Vol    =   
    0.453592 * 1E6,    # Volume  
  #PCT     =   
    0.3,                # percentage in plastic
  #C_fcm   =   
    0.2                  ## concentration in food contacting material Cp 
  )

product.DEHP <- c(
  #Vol     =   
    0.453592 * 152694720,    # Volume  
 #PCT     =   
    0.3,                # percentage in plastic
 # C_fcm   =   
    0.0005                  ## concentration in food contacting material Cp 
  )

product.DINCH <- c(
  #Vol    =   
  0.453592 * 6.1 * 2E6 ,    # Volume  
  #PCT     =   
  0.3,                # percentage in plastic
  #C_fcm   =   
  0.1                  ## concentration in food contacting material Cp 
  )

product.DPHP <- c(
  #Vol     =   
    0.453592 * 2.12e8 ,    # Volume  
  #PCT     =   
    0.3,                # percentage in plastic
  #C_fcm   =   
    0.03                  ## concentration in food contacting material Cp 
  )

exposure.DEP    <- Pred.Product.DEP  (product.DEP)
exposure.DEHA   <- Pred.Product.DEHA (product.DEHA)
exposure.DEHP   <- Pred.Product.DEHP (product.DEHP)
exposure.DINCH  <- Pred.Product.DINCH(product.DINCH)
exposure.DPHP   <- Pred.Product.DPHP (product.DPHP)

R.DEP   = Pred.DEP  (exposure.DEP)
R.DEHA  = Pred.DEHA (exposure.DEHA)
R.DEHP  = Pred.DEHP (exposure.DEHP)
R.DINCH = Pred.DINCH(exposure.DINCH)
R.DPHP  = Pred.DPHP (exposure.DPHP)

################## Create the matrix for normalized sensitivity coefficient data ##################
NSC_DEP    = matrix(nrow=1,ncol=2)
NSC_DEHA   = matrix(nrow=3,ncol=2)
NSC_DEHP   = matrix(nrow=3,ncol=2)
NSC_DINCH  = matrix(nrow=3,ncol=2)
NSC_DPHP   = matrix(nrow=3,ncol=2)



for (i in 1:3) {
  
  #product.DEHA.new       <- c(product.DEHA[i] *1.01,product.DEHA[-i])
  #product.DEHP.new       <- c(product.DEHP[i] *1.01,product.DEHP[-i])
  #product.DINCH.new      <- c(product.DINCH[i]*1.01,product.DINCH[-i])
  #product.DPHP.new       <- c(product.DPHP[i] *1.01,product.DPHP[-i])
  
  product.DEHA.new    <- product.DEHA
  new.DEHA            <- product.DEHA[i] *1.01
  product.DEHA.new[i] <- new.DEHA
  
  product.DEHP.new    <- product.DEHP
  new.DEHP            <- product.DEHP[i] *1.01
  product.DEHP.new[i] <- new.DEHP
  
  product.DINCH.new    <- product.DINCH
  new.DINCH            <- product.DINCH[i] *1.01
  product.DINCH.new[i] <- new.DINCH
  
  product.DPHP.new    <- product.DPHP
  new.DPHP            <- product.DPHP[i] *1.01
  product.DPHP.new[i] <- new.DPHP
  # Each cycle, generate a new value of parameter i, and delete parameter i, so that you can proceed to the next parameter i+1

  exposure.DEHA.new   <- Pred.Product.DEHA (product.DEHA.new)
  exposure.DEHP.new   <- Pred.Product.DEHP (product.DEHP.new)
  exposure.DINCH.new  <- Pred.Product.DINCH(product.DINCH.new)
  exposure.DPHP.new   <- Pred.Product.DPHP (product.DPHP.new)
  
  Rnew.DEHA                <- Pred.DEHA (exposure.DEHA.new)
  Rnew.DEHP                <- Pred.DEHP (exposure.DEHP.new)
  Rnew.DINCH               <- Pred.DINCH(exposure.DINCH.new)
  Rnew.DPHP                <- Pred.DPHP (exposure.DPHP.new)
  
  delta.P.DEHA        <- 100 #exp(pars.DEP[i])/(exp(pars.DEP[i])) * 1#00
  delta.P.DEHP        <- 100 #exp(pars.DEHP[i])/(exp(pars.DEHP[i])) * 1#00
  delta.P.DINCH       <- 100 #exp(pars.DEP[i])/(exp(pars.DEP[i])) * 1#00
  delta.P.DPHP        <- 100 
  
  
  ## Estimated the AUC
  ## DEHA
  DEHA.AUC.CA.new      =  Rnew.DEHA$outdf.DEHA %>% filter (Time == 14) %>% select (AUC_CA)
  DEHA.AUC.CA.ori      =  R.DEHA$outdf.DEHA    %>% filter (Time == 14) %>% select (AUC_CA)
  DEHA.AUC.CAM.new     =  Rnew.DEHA$outdf.DEHA %>% filter (Time == 14) %>% select (AUC_CAM)
  DEHA.AUC.CAM.ori     =  R.DEHA$outdf.DEHA    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DEHA    =  DEHA.AUC.CA.new  - DEHA.AUC.CA.ori
  delta.AUC.CAM.DEHA   =  DEHA.AUC.CAM.new - DEHA.AUC.CAM.ori

  ## DEHP
  DEHP.AUC.CA.new      =  Rnew.DEHP$outdf.DEHP %>% filter (Time == 14) %>% select (AUC_CA)
  DEHP.AUC.CA.ori      =  R.DEHP$outdf.DEHP    %>% filter (Time == 14) %>% select (AUC_CA)
  DEHP.AUC.CAM.new     =  Rnew.DEHP$outdf.DEHP %>% filter (Time == 14) %>% select (AUC_CAM)
  DEHP.AUC.CAM.ori     =  R.DEHP$outdf.DEHP    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DEHP    =  DEHP.AUC.CA.new  - DEHP.AUC.CA.ori
  delta.AUC.CAM.DEHP   =  DEHP.AUC.CAM.new - DEHP.AUC.CAM.ori

  ## DINCH
  DINCH.AUC.CA.new      =  Rnew.DINCH$outdf.DINCH %>% filter (Time == 14) %>% select (AUC_CA)
  DINCH.AUC.CA.ori      =  R.DINCH$outdf.DINCH    %>% filter (Time == 14) %>% select (AUC_CA)
  DINCH.AUC.CAM.new     =  Rnew.DINCH$outdf.DINCH %>% filter (Time == 14) %>% select (AUC_CAM)
  DINCH.AUC.CAM.ori     =  R.DINCH$outdf.DINCH    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DINCH    =  DINCH.AUC.CA.new  - DINCH.AUC.CA.ori
  delta.AUC.CAM.DINCH   =  DINCH.AUC.CAM.new - DINCH.AUC.CAM.ori
  
  ## DPHP
  DPHP.AUC.CA.new      =  Rnew.DPHP$outdf.DPHP %>% filter (Time == 14) %>% select (AUC_CA)
  DPHP.AUC.CA.ori      =  R.DPHP$outdf.DPHP    %>% filter (Time == 14) %>% select (AUC_CA)
  DPHP.AUC.CAM.new     =  Rnew.DPHP$outdf.DPHP %>% filter (Time == 14) %>% select (AUC_CAM)
  DPHP.AUC.CAM.ori     =  R.DPHP$outdf.DPHP    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DPHP    =  DPHP.AUC.CA.new  - DPHP.AUC.CA.ori
  delta.AUC.CAM.DPHP   =  DPHP.AUC.CAM.new - DPHP.AUC.CAM.ori

  #####
  NSC_DEHA   [i, 1]   <- as.numeric((delta.AUC.CA.DEHA/DEHA.AUC.CA.ori)   * delta.P.DEHA)
  NSC_DEHA   [i, 2]   <- as.numeric((delta.AUC.CAM.DEHA/DEHA.AUC.CAM.ori) * delta.P.DEHA)

  NSC_DEHP   [i, 1]   <- as.numeric((delta.AUC.CA.DEHP/DEHP.AUC.CA.ori)   * delta.P.DEHP)
  NSC_DEHP   [i, 2]   <- as.numeric((delta.AUC.CAM.DEHP/DEHP.AUC.CAM.ori) * delta.P.DEHP)
  
  NSC_DINCH   [i, 1]   <- as.numeric((delta.AUC.CA.DINCH/DINCH.AUC.CA.ori)   * delta.P.DINCH)
  NSC_DINCH   [i, 2]   <- as.numeric((delta.AUC.CAM.DINCH/DINCH.AUC.CAM.ori) * delta.P.DINCH)
  
  NSC_DPHP   [i, 1]   <- as.numeric((delta.AUC.CA.DPHP/DPHP.AUC.CA.ori)   * delta.P.DPHP)
  NSC_DPHP   [i, 2]   <- as.numeric((delta.AUC.CAM.DPHP/DPHP.AUC.CAM.ori) * delta.P.DPHP)
  
  
  cat("iteration = ", i , "\n") 

}

#print(pars.DEHA.new)


## DEP
#for (i in 1:2) {

  print(product.DEP)
  product.DEP.new    <- product.DEP *1.01
  
  #new.DEP            <- product.DEP *1.01
  #product.DEP.new[2] <- new.DEP
  #product.DEP.new    <- c(product.DEP[i] *1.01,product.DEP[-i])
  
  print(product.DEP.new)
  exposure.DEP.new   <- Pred.Product.DEP(product.DEP.new)
  Rnew.DEP           <- Pred.DEP (exposure.DEP.new)
 
  delta.P.DEP        <- 100 # Basically, the ratio is 100 for each parameter.

  ## Estimated the AUC
  
  DEP.AUC.CA.new      =  Rnew.DEP$outdf.DEP %>% filter (Time == 14) %>% select (AUC_CA)
  DEP.AUC.CA.ori      =  R.DEP$outdf.DEP    %>% filter (Time == 14) %>% select (AUC_CA)
  DEP.AUC.CAM.new     =  Rnew.DEP$outdf.DEP %>% filter (Time == 14) %>% select (AUC_CAM)
  DEP.AUC.CAM.ori     =  R.DEP$outdf.DEP    %>% filter (Time == 14) %>% select (AUC_CAM)
  
  delta.AUC.CA.DEP    =  DEP.AUC.CA.new  - DEP.AUC.CA.ori
  delta.AUC.CAM.DEP   =  DEP.AUC.CAM.new - DEP.AUC.CAM.ori
 
  NSC_DEP    [1, 1]   <- as.numeric((delta.AUC.CA.DEP/DEP.AUC.CA.ori)     * delta.P.DEP)
  NSC_DEP    [1, 2]   <- as.numeric((delta.AUC.CAM.DEP/DEP.AUC.CAM.ori)   * delta.P.DEP)
  
 # pars.DEHA.new  [i]   <- DEHA.old 
  
 # cat("iteration = ", i , "\n") 
#}

 
colnames (NSC_DEP)   = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DEP)   = product.DEP.names[2]
colnames (NSC_DEHA)  = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DEHA)  = product.DEHA.names[1:3]
colnames (NSC_DEHP)  = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DEHP)  = product.DEHP.names[1:3]
colnames (NSC_DINCH) = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DINCH) = product.DINCH.names[1:3]
colnames (NSC_DPHP)  = c("NSC_AUC_CA","NSC_AUC_CAM") 
rownames (NSC_DPHP)  = product.DPHP.names[1:3]


##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################
#install.packages("RColorBrewer")
library("RColorBrewer")
#install.packages("ggsci")
library("ggsci")


###################### Fig. 8a; AUC of plasma_parent ####################
melt.DEP.CA        = melt(NSC_DEP[,1])
melt.DEP.CA$group  = c("DEP")
melt.DEP.CA$par    = rownames(NSC_DEP)

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


#p                = Circle.plot (melt.data.CA%>%filter(abs(value)>0.01))
#p

DEHP_pcp <- c(as.numeric(0),"DEHP","C_pcp")
names(DEHP_pcp) <-  c("value","group","par")

DPHP_pcp <- c(as.numeric(0),"DPHP","C_pcp")
names(DPHP_pcp) <-  c("value","group","par")

DINCH_pcp <- c(as.numeric(0),"DINCH","C_pcp")
names(DINCH_pcp) <-  c("value","group","par")

DEHA_pcp <- c(as.numeric(0),"DEHA","C_pcp")
names(DEHA_pcp) <-  c("value","group","par")

DEP_vol <- c(as.numeric(0),"DEP","Vol")
names(DEP_vol) <-  c("value","group","par")

DEP_fcm <- c(as.numeric(0),"DEP","C_fcm")
names(DEP_fcm) <-  c("value","group","par")

data.CA  <- rbind(melt.data.CA,DEHP_pcp,DPHP_pcp,DINCH_pcp,DEP_vol,DEHA_pcp , DEP_fcm )  
data.CA

library(tidyverse)
data.CA = data.CA %>%
  arrange(value) %>% 
  mutate(par=gsub("C_pcp", "italic(C[pcp])", par),
         par=gsub("C_fcm", "italic(C[fcm])", par),
         par=gsub("Vol", "italic(Vol)", par),
         par=gsub("PCT", "italic(PCT)", par))



#data.CA$par[data.CA$par=="C_pcp"] = c("C[pcp]")

data.CA$value <- as.numeric(as.character(data.CA$value))

data.CA$group <- as.character(data.CA$group)
data.CA$group <- factor(data.CA$group, levels=c("DEP","DEHA", "DEHP", "DINCH","DPHP"))
data.CA
#data.CA$par <- as.character(data.CA$par)
#data.CA$par <- factor(data.CA$par, levels=c("C[fcm]","C[pcp]", "Vol", "PCT"))


################ plot ########################
#install.packages("remotes")                    # Install remotes package
#remotes::install_github("coolbutuseless/ggpattern")
#library("ggpattern")
library("ggrepel")
windowsFonts(Times=windowsFont("Times New Roman")) 

pd = position_dodge(width = 0.7)

p1 = 
  ggplot (data.CA, aes(reorder(as.factor(par), -as.numeric(abs(value*100))), fill = group, 
          as.numeric(abs(value*100)))) + 
  geom_bar(stat = "identity",width=0.63,position=pd,size = 1) +
  scale_fill_jco() +
 # geom_text(data=data.CA, aes(x=par, y=abs(value*100), label=par), parse = TRUE) +

  #scale_x_discrete(labels=rev(parse(text=unique(data.CA$par))))
   scale_x_discrete(labels=parse(text = c("italic(C[fcm])","italic(C[pcp])","italic(PCT)","italic(Vol)")))
 # scale_x_discrete(labels=parse(levels(as.factor(par))))
 # labs(x = "par")
  #geom_text( parse = TRUE)
  #scale_fill_manual(values = c("#00AFBB", "#FC4E07"))
  #scale_color_manual(values = c("#0A0A0A","#0A0A0A","#0A0A0A")) + #+
  #scale_y_continuous(minor_breaks = seq(-3 , 7, 1), breaks = seq(-3, 7,2))
#scale_y_discrete(expand = c(0, 0), expression(paste("Log [Urinary Conc] (", mu,"g/L)"))) +

p1 <- p1+ scale_y_continuous(limits = c(0,110),expand = c(0,0))

p1

p1 <- p1 + 
  
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major.y      = element_blank(),#element_line(colour = "grey", size = 0.5),
    panel.grid.minor.y      = element_blank(),#element_line(colour = "grey", linetype = 4, size = 0.3),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 20, colour = "black"),    # tick labels along axes 
    axis.title.y            = element_text (size   = 20, colour = "black"),   # label of axes
    legend.title            = element_blank(),
    legend.justification    =  c("right", "top"),
    axis.title.x            = element_blank(),
    legend.position         = c(0.99, 0.99),
    legend.text             = element_text (size = 20),
    legend.text.align       = 0) + 
  labs (y = "Sensitivity (%)") 

p1





####### Save plot #######
ggsave("Sens_product_parent_barABc.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot",
       width = 15, height = 10, units = "cm",dpi=320)


#dev.off()


melt.DEP.CAM        = melt(NSC_DEP[,2])
melt.DEP.CAM$group  = c("MEP")
melt.DEP.CAM$par    = rownames(NSC_DEP)

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


#p2                = Circle.plot (melt.data.CAM%>%filter(abs(value)>0.01))
#p2

ggsave("Sens_product_meta.tiff",scale = 1,
       plot = p2,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot",
       width = 47, height = 47, units = "cm",dpi=320)

