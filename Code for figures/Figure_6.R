## loading R packages
## Load libraries
library(mrgsolve) # Needed to run the main PBPK code
library(magrittr) # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)  # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(dplyr)    # Needed for the pipe %>% operator
library(reshape)  # melt function to reshape the table, reshape2 is version 2 of reshape. reshape is more stable when using melt function. We decided to use reshape for all.
library(gridExtra)
mrgsolve::loadso
loadso(mod)
dev.off()

rm(list=ls())
## Input PK model
DEP.PK.code     <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/humanPK.RDS')
DEHA.PK.code    <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/humanPK.RDS')
DEHP.PK.code    <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/humanPK.RDS')
DINCH.PK.code   <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/humanPK.RDS')
DPHP.PK.code    <- readRDS (file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/humanPK.RDS')

# load pk model
mod.DEP         <- mcode ("humanPK",   DEP.PK.code)
mod.DEHA        <- mcode ("humanPK",  DEHA.PK.code)
mod.DEHP        <- mcode ("humanPK",  DEHP.PK.code)
mod.DINCH       <- mcode ("humanPK", DINCH.PK.code)
mod.DPHP        <- mcode ("humanPK",  DPHP.PK.code)


## Loading observed data
## DEP urine MEP ug; cumulative (ug)
DEP_urine  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEP/Dermal_eva_std.csv")
DOSE_DEP_eva   = 17.1
DEP_urine.obs  <- cbind.data.frame (Time = DEP_urine$time,
                                    MEP  = DEP_urine$DEP/DOSE_DEP_eva*100,
                                    SD   = DEP_urine$STD/DOSE_DEP_eva*100)

## DEHA urine 5cx ug = Cal #3; cumulative (ug)
DOSE_DEHA_oral   = 10 * 1000
DEHA_urine  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHA/Eva1.csv")
DEHA_cx_urine.obs  <- cbind.data.frame (Time = DEHA_urine$time_cx,
                                    cx  = DEHA_urine$average_cx/DOSE_DEHA_oral*100,
                                    SD  = DEHA_urine$std_cx/DOSE_DEHA_oral*100)
DEHA_OH_urine.obs  <- cbind.data.frame (Time = DEHA_urine$time_oh,
                                   OH   = DEHA_urine$average_oh/DOSE_DEHA_oral*100,
                                   SD   = DEHA_urine$std_oh/DOSE_DEHA_oral*100)
DEHA_oxo_urine.obs <- cbind.data.frame (Time = DEHA_urine$time_oxo,
                                  oxo   = DEHA_urine$average_oxo/DOSE_DEHA_oral*100,
                                  SD    = DEHA_urine$std_oxo/DOSE_DEHA_oral*100)
DEHA_cx_urine.obs  <- na.omit(DEHA_cx_urine.obs)
DEHA_OH_urine.obs  <- na.omit(DEHA_OH_urine.obs)
DEHA_oxo_urine.obs <- na.omit(DEHA_oxo_urine.obs)


## DEHP plasma conc 
DEHP_plasma.obs <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2012_Plasma_DEHP.csv")
names(DEHP_plasma.obs)=c("Time", "DEHP")
PlasmaMEHP_2012 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2012_Plasma_MEHP.csv")
names(PlasmaMEHP_2012)=c("Time", "CA")
UrineMEHP_2012  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2012_Urine_MEHP.csv")
names(UrineMEHP_2012)=c("Time", "Urine")

MEHP_plasma <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/STD_Plasma_MEHP.csv")
MEHP_plasma.obs  <- cbind.data.frame (        Time  = MEHP_plasma$Time,
                                            MEHP    = MEHP_plasma$DEHP,
                                            SD      = MEHP_plasma$STD)

## DINCH MHNCH urine ug  = #cAL2 57 mg oral dose; 57 kg BW; cumulative (ug/kg bw)
DINCH_urine   <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DINCH/Eva1.csv")
BW_DINCH          = 75.7
DOSEoral_DINCH   = BW_DINCH * 1000   # ug Oral dose; 
DINCH_MHNCH_urine.obs  <- cbind.data.frame (Time  = DINCH_urine$time,
                                    MHNCH   = DINCH_urine$ave_MHNCH * BW_DINCH /DOSEoral_DINCH * 100,
                                    SD      = DINCH_urine$std_MHNCH* BW_DINCH /DOSEoral_DINCH * 100)
DINCH_MINCH_urine.obs  <- cbind.data.frame (Time = DINCH_urine$time,
                                   MINCH   = DINCH_urine$ave_MINCH * BW_DINCH /DOSEoral_DINCH *100,
                                   SD      = DINCH_urine$std_MINCH* BW_DINCH /DOSEoral_DINCH * 100)
DINCH_oxo_urine.obs    <- cbind.data.frame  (Time = DINCH_urine$time,
                                      oxo   = DINCH_urine$ave_oxo_MINCH * BW_DINCH /DOSEoral_DINCH *100,
                                      SD    = DINCH_urine$std_oxo_MINCH * BW_DINCH /DOSEoral_DINCH * 100)
DINCH_cx_urine.obs     <- cbind.data.frame   (Time = DINCH_urine$time,
                                        cx   = DINCH_urine$ave_cx_MINCH * BW_DINCH /DOSEoral_DINCH *100,
                                        SD   = DINCH_urine$std_cx_MINCH * BW_DINCH /DOSEoral_DINCH * 100)
DINCH_CHDA_urine.obs   <- cbind.data.frame  (Time = DINCH_urine$time,
                                     CHDA   = DINCH_urine$ave_CHDA * BW_DINCH /DOSEoral_DINCH *100,
                                     SD     = DINCH_urine$std_CHDA * BW_DINCH /DOSEoral_DINCH * 100)


## DPHP plasma conc average + std
DPHP_plasma <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DPHP/Eva_DPHP_plasma.csv")
MPHP_plasma <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DPHP/Eva_MPHP_plasma.csv")
OH_plasma <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DPHP/Eva_OH_plasma.csv")
oxo_plasma <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DPHP/Eva_oxo_plasma.csv")

MPHP_urine  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DPHP/Eva_MPHP_urine.csv")
OH_urine    <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DPHP/Eva_OH_urine.csv")
oxo_urine   <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DPHP/Eva_oxo_urine.csv")

Eva_plasma_DPHP.obs <- cbind.data.frame (Time  = round(DPHP_plasma$time,2),
                                      DPHP = DPHP_plasma$DPHP,
                                        SD = DPHP_plasma$std)

Eva_plasma_MPHP.obs <- cbind.data.frame (Time  = round(MPHP_plasma$time,2),
                                     MPHP  = MPHP_plasma$MPHP,
                                        SD = MPHP_plasma$std)

Eva_plasma_OH.obs <- cbind.data.frame (Time  = round(OH_plasma$time,2),
                                         OH  = OH_plasma$OH,
                                          SD = OH_plasma$std)

Eva_plasma_oxo.obs <- cbind.data.frame (Time  = round(oxo_plasma$time,2),
                                       oxo    = oxo_plasma$oxo,
                                       SD     = oxo_plasma$std)

Eva_urine_MPHP.obs <- cbind.data.frame (Time  = round(MPHP_urine$time,2),
                                         MPHP = MPHP_urine$MPHP,
                                          SD  = MPHP_urine$std)

Eva_urine_OH.obs <- cbind.data.frame (Time  = round(OH_urine$time,2),
                                        OH  = OH_urine$OH,
                                        SD  = OH_urine$std)

Eva_urine_oxo.obs <- cbind.data.frame (Time  = round(oxo_urine$time,2),
                                        oxo = oxo_urine$oxo,
                                        SD  =  oxo_urine$std)


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

## Loading MCMC data
DEP.MCMC         <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/Human.MCMC.rds")
DEHA.MCMC        <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/Human.MCMC.rds")
DEHP.MCMC        <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/Human.MCMC.rds")
DINCH.MCMC       <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/Human.MCMC.rds")
DPHP.MCMC        <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/Human.MCMC.rds")

## Theme
windowsFonts(Times=windowsFont("Times New Roman"))

Theme.Fig <-theme(
  legend.position         = "none",
  legend.title            = element_blank(),
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border            = element_rect(colour = "black", fill=NA, size=3),
  panel.background        = element_rect (fill="white"),
  panel.grid              = element_blank(),
  axis.text               = element_text (size   = 25, colour = "black", face = "bold"),    # tick labels along axes 
  axis.line               = element_line(size    = 0.8, colour = "black"),
  axis.title              = element_text (size   = 25, colour = "black", face = "bold"),
  plot.title              = element_text(color="black", size=25,face="bold"))   # label of axes


Theme.Figb <-theme(
  legend.position         = "none",
  legend.title            = element_blank(),
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border            = element_rect(colour = "black", fill=NA, size=3),
  panel.background        = element_rect (fill="white"),
  panel.grid              = element_blank(),
  axis.text               = element_text (size   = 25, colour = "black", face = "bold"),    # tick labels along axes 
  axis.line               = element_line(size    = 0.8, colour = "black"),
  axis.title              = element_blank(),
  plot.title              = element_text(color="black", size=25,face="bold"))   # label of axes


Theme.Fig2 <-theme(
  legend.position         = "right",
  legend.title            = element_blank(),
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border            = element_rect(colour = "black", fill=NA, size=3),
  panel.background        = element_rect (fill="white"),
  panel.grid              = element_blank(),
  axis.text               = element_text (size   = 25, colour = "black", face = "bold"),    # tick labels along axes 
  axis.line               = element_line(size    = 0.8, colour = "black"),
  axis.title              = element_text (size   = 25, colour = "black", face = "bold"))   # label of axes



############################# prediction function ####################################
############################# DEP data ###############################################
######################################################################################

pred.DEP <- function(pars.DEP,pred=FALSE) {
  
  
  ## Get out of log domain
  pars.DEP %<>% lapply(exp)
  names(pars.DEP) <- names(pars.DEP)
  pars.DEP        <- pars.DEP[-which_sig.DEP]
  
  ## Exposure scenario for dermal exposue 
  tinterval             = 24                            ## hr, Time interval
  TDoses                = 1                             ## Dose times
  ## Set up the exposure time
  tsamp.DEP = tgrid(0,48,0.1)  # 48 hours
  ### DEP
  DOSE_DEP_eva   = 17.1
  Der_DEP    <- ev(ID=1, amt= DOSE_DEP_eva, ii=tinterval, addl=TDoses-1, cmt="AC", tinf = 3, replicate = FALSE)
  ex_DEP  <- Der_DEP 
  

  out.DEP <- 
    mod.DEP %>% 
    param(pars.DEP, BW = 85.7) %>%
    Req(Urine)%>%
    update(atol = 1E-8,maxsteps = 1000000) %>%
    mrgsim_d(data = ex_DEP, tgrid=tsamp.DEP)
  
  out.DEP <- cbind.data.frame (Time =out.DEP$time,
                              AU    =(out.DEP$Urine)/DOSE_DEP_eva*100)
  

  return (list( "out.DEP"  = out.DEP))
  
  
}

pars.DEP          = DEP.MCMC$bestpar[1:6]
fit.A             = pred.DEP (pars.DEP)$out.DEP

Newtime.r         = pred.DEP(pars.DEP)$out.DEP$Time
nrwo.r            = length (Newtime.r)

# Create the matrix 
MC.DEP.Curine    = matrix(nrow = nrwo.r, ncol = 2000)


## Input paramters
for(i in 1:2000){
  
  j = i * 2  
  pars                 = DEP.MCMC$pars    [j,]     # sample parameter set once every ten sets, so you will have 1000 sets from 1000 total sets
  
  MCdata               = pred.DEP (pars[1:6])
  
  MC.DEP.Curine    [,i]  = MCdata$out.DEP$AU
  
  
  cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
}


############# plot #############

MC.DEP.Curine.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DEP.Curine, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)
################################################
#### reset graphics
#dev.off()

################## plot ########################

p2.DEP.Curine <- 
  ggplot() + 
  geom_ribbon(data = MC.DEP.Curine.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DEP.Curine.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DEP.Curine.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "gold4") +
  geom_line(data= MC.DEP.Curine.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = DEP_urine.obs , aes(x=Time, y= MEP), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DEP_urine.obs , aes(x=Time,ymin= MEP-SD, ymax = MEP+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black") +
  ggtitle("MEP (DEP)")

p2.DEP.Curine

p2.DEP.Curine  <- p2.DEP.Curine +
  scale_y_continuous(limits = c(0,90),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,24),expand = c(0,0))

## Modified the theme
## plot

p2.DEP.Curine  = p2.DEP.Curine  + Theme.Figb + labs (x ="Time (hours)")


#p1.r   = p1.r + Theme.Fig + labs (x ="Dose (mg/kg-day)", y="")


p2.DEP.Curine


### Save figure #####

ggsave("Fig.DEP.tiff",scale = 1,
       plot = p2.DEP.Curine,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEP/PK/Code for plot",
       width = 20, height = 15, units = "cm",dpi=320)






############################# prediction function ####################################
############################# DEHA data ##############################################
######################################################################################


## Preiction function
pred.DEHA <- function(pars.DEHA,pred=FALSE) {
  
  
  ## Set up the exposure time
  tsamp.DEHA = tgrid(0,48,0.01)  # 48 hours
  tinterval             = 24                            ## hr, Time interval
  TDoses                = 1 
  
  ## Get out of log domain
  pars.DEHA %<>% lapply(exp)
  names(pars.DEHA) <- names(pars.DEHA)
  pars.DEHA        <- pars.DEHA[-which_sig.DEHA]
  
  ### DEHA
  DOSE_DEHA_oral   = 10 * 1000   # ug Oral dose; 
  Eva_DEHA    <- ev(ID=1, amt= DOSE_DEHA_oral, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_DEHA  <- Eva_DEHA
  
  out.DEHA <- 
    mod.DEHA %>% 
    param(pars.DEHA) %>%
    Req(Urine_OH, Urine_cx, Urine_oxo)%>%
    update(atol = 1E-8,rtol= 1e-8,maxsteps = 50000) %>%
    mrgsim_d(data = ex_DEHA, tgrid=tsamp.DEHA)
  
  out.DEHA <- cbind.data.frame (Time = out.DEHA$time,
                                OH   = (out.DEHA$Urine_OH)/DOSE_DEHA_oral*100,
                                cx   = (out.DEHA$Urine_cx)/DOSE_DEHA_oral*100,
                               oxo   = (out.DEHA$Urine_oxo)/DOSE_DEHA_oral*100)
  
  
  return (list( "out.DEHA"  = out.DEHA))
  
  
}

pars.DEHA         = DEHA.MCMC$bestpar[1:13]
fit.B             = pred.DEHA (pars.DEHA)$out.DEHA

Newtime.r         = pred.DEHA(pars.DEHA)$out.DEHA$Time
nrwo.r            = length (Newtime.r)

# Create the matrix 

MC.DEHA.cx        = matrix(nrow = nrwo.r, ncol = 2000)
MC.DEHA.OH        = matrix(nrow = nrwo.r, ncol = 2000)
MC.DEHA.oxo       = matrix(nrow = nrwo.r, ncol = 2000)


## Input paramters
#quan.Human  = summary(as.mcmc(DEHA.MCMC$pars))$quantiles  
#exp(quan.Human)

for(i in 1:2000){
  
  j = i *1  
  pars                 = DEHA.MCMC$pars    [j,]     # sample parameter set once every ten sets, so you will have 1000 sets from 1000 total sets
  
  MCdata               = pred.DEHA (pars[1:13])
  
  MC.DEHA.cx       [,i]  = MCdata $out.DEHA$cx
  MC.DEHA.OH       [,i]  = MCdata $out.DEHA$OH    # this dataset contains 5000 columns of CA. out.rat.a is one of the output objects of the function.
  MC.DEHA.oxo      [,i]  = MCdata $out.DEHA$oxo
  
  
  cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
}


############# plot #############


MC.DEHA.cx.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DEHA.cx, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.DEHA.OH.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DEHA.OH, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.DEHA.oxo.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DEHA.oxo, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)
################################################


MC.DEHA.cx  <- 
  ggplot() + 
  geom_ribbon(data = MC.DEHA.cx.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DEHA.cx.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DEHA.cx.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.DEHA.cx.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data= DEHA_cx_urine.obs, aes(x=Time, y= cx), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DEHA_cx_urine.obs, aes(x=Time,ymin= cx-SD, ymax = cx+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("5cx-MEPA (DEHA)")

MC.DEHA.cx

MC.DEHA.OH <- 
  ggplot() + 
  geom_ribbon(data = MC.DEHA.OH.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DEHA.OH.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DEHA.OH.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "black") +
  geom_line(data= MC.DEHA.OH.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = DEHA_OH_urine.obs, aes(x=Time, y= OH), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DEHA_OH_urine.obs, aes(x=Time,ymin= OH-SD, ymax = OH+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("5OH-MEHA (DEHA)")
  
MC.DEHA.OH

MC.DEHA.oxo <- 
  ggplot() + 
  geom_ribbon(data = MC.DEHA.oxo.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DEHA.oxo.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DEHA.oxo.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "gold4") +
  geom_line(data= MC.DEHA.oxo.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = DEHA_oxo_urine.obs, aes(x=Time, y= oxo), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DEHA_oxo_urine.obs, aes(x=Time,ymin= oxo-SD, ymax = oxo+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("5oxo-MEHA (DEHA)")
MC.DEHA.oxo

MC.DEHA.cx  <- MC.DEHA.cx +
  scale_y_continuous(limits = c(0,0.3),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,48),expand = c(0,0))

MC.DEHA.OH  <- MC.DEHA.OH +
  scale_y_continuous(limits = c(0,0.115),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,29),expand = c(0,0))

MC.DEHA.oxo  <- MC.DEHA.oxo +
  scale_y_continuous(limits = c(0,0.3),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,48),expand = c(0,0))

## Modified the theme
## plot

MC.DEHA.cx    = MC.DEHA.cx    + Theme.Fig + labs (x ="Time (hours)", y="Cumulative excretion [%]")
MC.DEHA.OH    = MC.DEHA.OH    + Theme.Figb + labs (x ="Time (hours)")
MC.DEHA.oxo   = MC.DEHA.oxo   + Theme.Fig + labs (x ="Time (hours)",  y="Cumulative excretion [%]")


#p1.r   = p1.r + Theme.Fig + labs (x ="Dose (mg/kg-day)", y="")


MC.DEHA.cx
MC.DEHA.OH
MC.DEHA.oxo 


### Save figure #####

ggsave("MC.DEHA.cx.tiff",scale = 1,
       plot = MC.DEHA.cx,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("MC.DEHA.OH.tiff",scale = 1,
       plot = MC.DEHA.OH,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("MC.DEHA.oxo.tiff",scale = 1,
       plot = MC.DEHA.oxo,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result",
       width = 20, height = 15, units = "cm",dpi=320)



############################# prediction function ####################################
############################# DEHP data ##############################################
######################################################################################


## Preiction function
pred.DEHP <- function(pars.DEHP,pred=FALSE) {
  
  tinterval             = 24                            ## hr, Time interval
  TDoses                = 1                             ## Dose times
  ## Exposure scenario for oral exposue to 2 mg/kg
  
  DOSEoral2012         = 50650                  ## Amount of oral dose average 2004, 2012 (average)
  ex2012               <- ev (ID=1, amt= DOSEoral2012, 
                              ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  
  ## Exposure scenario for oral exposue to 15 mg/kg
  
  
  
  ## Set up the exposure time
  tsamp.DEHP = tgrid(0,48,0.1)  # 48 hours
  
  ## Get out of log domain
  pars.DEHP %<>% lapply(exp)
  names(pars.DEHP) <- names(pars.DEHP)
  pars.DEHP        <- pars.DEHP[-which_sig.DEHP]
  
  
  ## Exposure scenario for oral exposue 
  out.2012 <- 
    mod.DEHP %>% 
    param(pars.DEHP, BW =82) %>%
    Req (Plasma, Plasma_MEHP, Urine)%>%
    update(atol = 1E-8,rtol= 1e-8,maxsteps = 50000) %>%
    mrgsim_d (data = ex2012, tgrid = tsamp.DEHP)
  
  out.2012 <- cbind.data.frame(Time        = out.2012$time, 
                                 CA_DEHP     = out.2012$Plasma,
                                 CA_MEHP     = out.2012$Plasma_MEHP,
                                 Urine_MEHP  = (out.2012$Urine/DOSEoral2012)*100)        ## %Dose, PFOS in urine
  
  
  return (list( "out.2012"  = out.2012))
  
}

pars.DEHP         = DEHP.MCMC$bestpar[1:6]
fit.B             = pred.DEHP (pars.DEHP)$out.2012

Newtime.r         = pred.DEHP(pars.DEHP)$out.2012$Time
nrwo.r            = length (Newtime.r)

# Create the matrix 

MC.b.CA        = matrix(nrow = nrwo.r, ncol = 2000)
MC.b.CA_DEHP   = matrix(nrow = nrwo.r, ncol = 2000)
MC.b.Curine    = matrix(nrow = nrwo.r, ncol = 2000)


## Input paramters

for(i in 1:2000){
  
  DOSEoral2012         = 50650    
  j = i *1  
  pars                 = DEHP.MCMC$pars    [j,]     # sample parameter set once every ten sets, so you will have 1000 sets from 1000 total sets
  
  MCdata               = pred.DEHP (pars[1:6])
 
  MC.b.CA        [,i]  = MCdata $out.2012$CA_MEHP
  MC.b.CA_DEHP   [,i]  = MCdata $out.2012$CA_DEHP    # this dataset contains 5000 columns of CA. out.rat.a is one of the output objects of the function.
  MC.b.Curine    [,i]  = MCdata $out.2012$Urine_MEHP
  
  
  cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
}


############# plot #############


MC.b.CA.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.b.CA, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.b.CA_DEHP.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.b.CA_DEHP, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.b.Curine.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.b.Curine, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)
################################################


p2.b.CA <- 
  ggplot() + 
  geom_ribbon(data = MC.b.CA.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="steelblue1", alpha=0.3) +
  geom_ribbon(data = MC.b.CA.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="steelblue3", alpha = 0.5) +
  geom_line(data= MC.b.CA.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.b.CA.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data= MEHP_plasma.obs, aes(x=Time, y= MEHP), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = MEHP_plasma.obs, aes(x=Time,ymin= MEHP-SD, ymax = MEHP+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle ("MEHP (DEHP)")
p2.b.CA 

p2.b.CA_DEHP <- 
  ggplot() + 
  geom_ribbon(data = MC.b.CA_DEHP.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="steelblue1", alpha=0.3) +
  geom_ribbon(data = MC.b.CA_DEHP.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="steelblue3", alpha = 0.5) +
  geom_line(data= MC.b.CA_DEHP.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.b.CA_DEHP.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data=DEHP_plasma.obs, aes(x=Time, y= DEHP), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  ggtitle ("DEHP")
p2.b.CA_DEHP

p2.b.Curine <- 
  ggplot() + 
  geom_ribbon(data = MC.b.Curine.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.b.Curine.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.b.Curine.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "gold4") +
  geom_line(data= MC.b.Curine.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data=UrineMEHP_2012, aes(x=Time, y= Urine), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  ggtitle ("MEHP (DEHP)")


p2.b.Curine


p2.b.CA  <- p2.b.CA +
  scale_y_continuous(limits = c(-10,5300),expand = c(0,0))+ 
  scale_x_continuous(limits = c(-0.5,24),expand = c(0,0))
p2.b.CA

p2.b.CA_DEHP  <- p2.b.CA_DEHP +
  scale_y_continuous(limits = c(-10,450),expand = c(0,0))+ 
  scale_x_continuous(limits = c(-0.5,24),expand = c(0,0))

p2.b.Curine  <- p2.b.Curine +
  scale_y_continuous(limits = c(0,10),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,24),expand = c(0,0))

## Modified the theme
## plot

p2.b.CA  = p2.b.CA  + Theme.Figb + labs (x ="Time (hours)")
p2.b.CA_DEHP  = p2.b.CA_DEHP  + Theme.Fig + labs (x ="Time (hours)",y=expression(paste("Concentration in plasma (", mu,"g/L)")))
p2.b.Curine  = p2.b.Curine  + Theme.Fig + labs (x ="Time (hours)",  y="Cumulative excretion [%]")


#p1.r   = p1.r + Theme.Fig + labs (x ="Dose (mg/kg-day)", y="")


p2.b.CA
p2.b.CA_DEHP 
p2.b.Curine


### Save figure #####

ggsave("p2.b.CA.tiff",scale = 1,
       plot = p2.b.CA,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("Fp2.b.CA_DEHP.tiff",scale = 1,
       plot = p2.b.CA_DEHP,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("p2.b.Curine.tiff",scale = 1,
       plot = p2.b.Curine,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result",
       width = 20, height = 15, units = "cm",dpi=320)


grid.arrange(p2.DEP.Curine, MC.DEHA.OH, p2.b.CA,nrow = 1)
ggsave("Fig.TIME_COURSE1.tiff",scale = 1,
       plot = grid.arrange(p2.DEP.Curine, MC.DEHA.OH, p2.b.CA,nrow = 1),
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot",
       width = 36, height = 12, units = "cm",dpi=320)

############################# prediction function ####################################
############################# DINCH data ##############################################
######################################################################################

## Preiction function
pred.DINCH <- function(pars.DINCH,pred=FALSE) {
  
  tinterval             = 24                            ## hr, Time interval
  TDoses                = 1  
  ## Set up the exposure time
  tsamp.DINCH = tgrid(0,72,0.01)  # 48 hours
  
  ## Get out of log domain
  pars.DINCH %<>% lapply(exp)
  names(pars.DINCH) <- names(pars.DINCH)
  pars.DINCH        <- pars.DINCH[-which_sig.DINCH]
  
  #### DINCH
  BW_DINCH          = 75.7
  DOSEoral_DINCH   = BW_DINCH * 1000   # ug Oral dose; 
  
  Eva_DINCH    <- ev(ID=1, amt= DOSEoral_DINCH, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_DINCH  <- Eva_DINCH
  
  
  out.DINCH <-
    
    mod.DINCH %>% 
    param(pars.DINCH,BW = 57) %>%  
    Req(Urine_MHNCH, Urine_MINCH, Urine_cx, Urine_oxo, Urine_CHDA)%>%
    update(atol = 1E-8,rtol= 1e-8,maxsteps = 50000) %>%
    mrgsim_d(data = ex_DINCH, tgrid=tsamp.DINCH)
  
  out.DINCH <-cbind.data.frame (Time = out.DINCH$time,
                             MHNCH   = (out.DINCH$Urine_MHNCH)/DOSEoral_DINCH * 100,
                             MINCH   = (out.DINCH$Urine_MINCH)/DOSEoral_DINCH* 100,
                             cx      = (out.DINCH$Urine_cx)/DOSEoral_DINCH * 100,
                             oxo     = (out.DINCH$Urine_oxo)/DOSEoral_DINCH * 100,
                             CHDA    = (out.DINCH$Urine_CHDA)/DOSEoral_DINCH* 100)
  
  
  return (list( "out.DINCH"  = out.DINCH))
  
  
}

pars.DINCH         = DINCH.MCMC$bestpar[1:18]
fit.B              = pred.DINCH (pars.DINCH)$out.DINCH

Newtime.r          = pred.DINCH(pars.DINCH)$out.DINCH$Time
nrwo.r             = length (Newtime.r)

# Create the matrix 

MC.DINCH.MHNCH      = matrix(nrow = nrwo.r, ncol = 2000)
MC.DINCH.MINCH      = matrix(nrow = nrwo.r, ncol = 2000)
MC.DINCH.cx         = matrix(nrow = nrwo.r, ncol = 2000)
MC.DINCH.oxo        = matrix(nrow = nrwo.r, ncol = 2000)
MC.DINCH.CHDA       = matrix(nrow = nrwo.r, ncol = 2000)


## Input paramters

for(i in 1:2000){
  
  BW_DINCH          = 75.7
  DOSEoral_DINCH    = BW_DINCH * 1000   # ug Oral dose; 
  
  j = i * 2  
  pars                 = DINCH.MCMC$pars    [j,]     # sample parameter set once every ten sets, so you will have 1000 sets from 1000 total sets
  
  MCdata               = pred.DINCH (pars[1:18])
  
  MC.DINCH.MHNCH     [,i]  = MCdata $out.DINCH$MHNCH
  MC.DINCH.MINCH     [,i]  = MCdata $out.DINCH$MINCH   
  MC.DINCH.cx        [,i]  = MCdata $out.DINCH$cx
  MC.DINCH.oxo       [,i]  = MCdata $out.DINCH$oxo   # this dataset contains 5000 columns of CA. out.rat.a is one of the output objects of the function.
  MC.DINCH.CHDA      [,i]  = MCdata $out.DINCH$CHDA
  
  
  cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
}


############# plot #############

MC.DINCH.MHNCH.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DINCH.MHNCH, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.DINCH.MINCH.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DINCH.MINCH, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.DINCH.cx.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DINCH.cx, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.DINCH.oxo.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DINCH.oxo, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.DINCH.CHDA.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.DINCH.CHDA, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)
################################################

MC.DINCH.MHNCH  <- 
  ggplot() + 
  geom_ribbon(data = MC.DINCH.MHNCH.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DINCH.MHNCH.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DINCH.MHNCH.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "black") +
  geom_line(data= MC.DINCH.MHNCH.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data= DINCH_MHNCH_urine.obs, aes(x=Time, y= MHNCH), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DINCH_MHNCH_urine.obs, aes(x=Time,ymin= MHNCH-SD, ymax = MHNCH+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("MHNCH (DINCH)")
MC.DINCH.MHNCH

MC.DINCH.MINCH  <- 
  ggplot() + 
  geom_ribbon(data = MC.DINCH.MINCH.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DINCH.MINCH.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DINCH.MINCH.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.DINCH.MINCH.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data= DINCH_MINCH_urine.obs, aes(x=Time, y= MINCH), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DINCH_MINCH_urine.obs, aes(x=Time,ymin= MINCH-SD, ymax = MINCH+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("MINCH (DINCH)")
MC.DINCH.MINCH


MC.DINCH.cx  <- 
  ggplot() + 
  geom_ribbon(data = MC.DINCH.cx.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DINCH.cx.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DINCH.cx.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.DINCH.cx.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data= DINCH_cx_urine.obs, aes(x=Time, y= cx), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DINCH_cx_urine.obs, aes(x=Time,ymin= cx-SD, ymax = cx+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("cx-MINCH (DINCH)")
MC.DINCH.cx

MC.DINCH.oxo <- 
  ggplot() + 
  geom_ribbon(data = MC.DINCH.oxo.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DINCH.oxo.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DINCH.oxo.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.DINCH.oxo.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = DINCH_oxo_urine.obs, aes(x=Time, y= oxo), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DINCH_oxo_urine.obs, aes(x=Time,ymin= oxo-SD, ymax = oxo+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black") +
  ggtitle("oxo-MINCH (DINCH)")
MC.DINCH.oxo

MC.DINCH.CHDA <- 
  ggplot() + 
  geom_ribbon(data = MC.DINCH.CHDA.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.DINCH.CHDA.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.DINCH.CHDA.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "gold4") +
  geom_line(data= MC.DINCH.CHDA.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = DINCH_CHDA_urine.obs, aes(x=Time, y= CHDA), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = DINCH_CHDA_urine.obs, aes(x=Time,ymin= CHDA-SD, ymax = CHDA+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("CHDA (DINCH)")
MC.DINCH.CHDA 

MC.DINCH.MHNCH  <- MC.DINCH.MHNCH +
  scale_y_continuous(limits = c(0,14),expand = c(0,0))+ 
  scale_x_continuous(limits = c(-0,48),expand = c(0,0))
MC.DINCH.MHNCH 

MC.DINCH.MINCH  <- MC.DINCH.MINCH +
  scale_y_continuous(limits = c(0,10),expand = c(0,0))+ 
  scale_x_continuous(limits = c(-0,48),expand = c(0,0))
MC.DINCH.MINCH

MC.DINCH.cx  <- MC.DINCH.cx +
  scale_y_continuous(limits = c(0,7),expand = c(0,0))+ 
  scale_x_continuous(limits = c(-0,48),expand = c(0,0))
MC.DINCH.cx

MC.DINCH.oxo  <- MC.DINCH.oxo +
  scale_y_continuous(limits = c(0,10),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,48),expand = c(0,0))
MC.DINCH.oxo

MC.DINCH.CHDA  <- MC.DINCH.CHDA +
  scale_y_continuous(limits = c(0,100),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,48),expand = c(0,0))
MC.DINCH.CHDA
## Modified the theme
## plot
MC.DINCH.MHNCH    = MC.DINCH.MHNCH + Theme.Figb + labs (x ="Time (hours)")
MC.DINCH.MINCH    = MC.DINCH.MINCH + Theme.Fig + labs (x ="Time (hours)", y="Cumulative excretion [%]")
MC.DINCH.cx       = MC.DINCH.cx + Theme.Fig + labs (x ="Time (hours)", y="Cumulative excretion [%]")
MC.DINCH.oxo      = MC.DINCH.oxo   + Theme.Fig + labs (x ="Time (hours)",  y="Cumulative excretion [%]")
MC.DINCH.CHDA     = MC.DINCH.CHDA  + Theme.Fig + labs (x ="Time (hours)", y="Cumulative excretion [%]")


#p1.r   = p1.r + Theme.Fig + labs (x ="Dose (mg/kg-day)", y="")

MC.DINCH.MHNCH 
MC.DINCH.MINCH 
MC.DINCH.cx
MC.DINCH.CHDA
MC.DINCH.oxo 


### Save figure #####

ggsave("MC.DINCH.MHNCH.tiff",scale = 1,
       plot = MC.DINCH.MHNCH,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result",
       width = 20, height = 15, units = "cm",dpi=320)

ggsave("MC.DINCH.MINCH.tiff",scale = 1,
       plot = MC.DINCH.MINCH,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("MC.DINCH.cx.tiff",scale = 1,
       plot = MC.DINCH.cx,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("MC.DINCH.oxo.tiff",scale = 1,
       plot = MC.DINCH.oxo,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("MC.DINCH.CHDA.tiff",scale = 1,
       plot = MC.DINCH.CHDA,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result",
       width = 20, height = 15, units = "cm",dpi=320)




############################# prediction function ####################################
############################# DPHP data ##############################################
######################################################################################


## Preiction function
pred.DPHP <- function(pars.DPHP,pred=FALSE) {
  
  
  ## Exposure scenario for oral exposue to 2 mg/kg
  
  tinterval             = 24                            ## hr, Time interval
  TDoses                = 1                             ## Dose times
  
  ## Set up the exposure time
  tsamp.DPHP = tgrid(0,48,0.01)  # 48 hours
  
  ## Get out of log domain
  pars.DPHP %<>% lapply(exp)
  names(pars.DPHP) <- names(pars.DPHP)
  pars.DPHP        <- pars.DPHP[-which_sig.DPHP]
  
  
  ##### DPHP average 738 ug/kg bw
  BW_DPHP          = 84    # AVERAGE WEIGHT
  DOSEoral_DPHP    = 738 * BW_DPHP  # ug Oral dose; 
  
  Eva_DPHP    <- ev(ID=1, amt= DOSEoral_DPHP, ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  ex_DPHP     <- Eva_DPHP
  
  out.DPHP <- 
    mod.DPHP %>% 
    param(pars.DPHP, BW = 84) %>%
    Req(Plasma_DPHP, Plasma_MPHP, Plasma_OH, Plasma_oxo, Urine_MPHP, Urine_OH, Urine_oxo)%>%
    update(atol = 1E-8,rtol= 1e-8,maxsteps = 50000) %>%
    mrgsim_d(data = ex_DPHP, tgrid=tsamp.DPHP)
  
  out.DPHP <- cbind.data.frame (Time  =out.DPHP$time,
                        Plasma_DPHP   = out.DPHP$Plasma_DPHP,
                        Plasma_MPHP   = out.DPHP$Plasma_MPHP,
                        Plasma_OH     = out.DPHP$Plasma_OH,
                        Plasma_oxo    = out.DPHP$Plasma_oxo,
                        Urine_MPHP    = (out.DPHP$Urine_MPHP)/DOSEoral_DPHP*100,
                        Urine_OH      = (out.DPHP$Urine_OH)/DOSEoral_DPHP*100,
                        Urine_oxo     = (out.DPHP$Urine_oxo)/DOSEoral_DPHP*100)

  
  return (list( "out.DPHP"  = out.DPHP))
  
  
}

pars.DPHP         = DPHP.MCMC$bestpar[1:14]
fit.B             = pred.DPHP (pars.DPHP)$out.DPHP

Newtime.r         = pred.DPHP(pars.DPHP)$out.DPHP$Time
nrwo.r            = length (Newtime.r)

# Create the matrix 

MC.plasma.DPHP       = matrix(nrow = nrwo.r, ncol = 2000)
MC.plasma.MPHP       = matrix(nrow = nrwo.r, ncol = 2000)
MC.plasma.OH         = matrix(nrow = nrwo.r, ncol = 2000)
MC.plasma.oxo        = matrix(nrow = nrwo.r, ncol = 2000)
 
MC.urine.MPHP        = matrix(nrow = nrwo.r, ncol = 2000)
MC.urine.OH          = matrix(nrow = nrwo.r, ncol = 2000)
MC.urine.oxo         = matrix(nrow = nrwo.r, ncol = 2000)

## Input paramters

for(i in 1:2000){
  
  j = i *2  
  pars                 = DPHP.MCMC$pars    [j,]     # sample parameter set once every ten sets, so you will have 1000 sets from 1000 total sets
  
  MCdata               = pred.DPHP (pars[1:14])
  
  MC.plasma.DPHP   [,i]  = MCdata $out.DPHP$Plasma_DPHP
  MC.plasma.MPHP   [,i]  = MCdata $out.DPHP$Plasma_MPHP    # this dataset contains 5000 columns of CA. out.rat.a is one of the output objects of the function.
  MC.plasma.OH     [,i]  = MCdata $out.DPHP$Plasma_OH
  MC.plasma.oxo    [,i]  = MCdata $out.DPHP$Plasma_oxo
  
  MC.urine.MPHP    [,i]  = MCdata $out.DPHP$Urine_MPHP
  MC.urine.OH      [,i]  = MCdata $out.DPHP$Urine_OH
  MC.urine.oxo     [,i]  = MCdata $out.DPHP$Urine_oxo
  
  
  cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
}


############# plot #############


MC.plasma.DPHP.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.plasma.DPHP, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.plasma.MPHP.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.plasma.MPHP, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.plasma.OH.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.plasma.OH, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.plasma.oxo.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.plasma.oxo, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)



MC.urine.MPHP.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.urine.MPHP, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.urine.OH.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.urine.OH, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

MC.urine.oxo.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.urine.oxo, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

################################################


MC.plasma.DPHP<- 
  ggplot() + 
  geom_ribbon(data = MC.plasma.DPHP.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="steelblue1", alpha=0.3) +
  geom_ribbon(data = MC.plasma.DPHP.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="steelblue3", alpha = 0.5) +
  geom_line(data= MC.plasma.DPHP.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.plasma.DPHP.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = Eva_plasma_DPHP.obs, aes(x=Time, y= DPHP), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Eva_plasma_DPHP.obs, aes(x=Time,ymin= DPHP-SD, ymax = DPHP+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("DPHP")
MC.plasma.DPHP

MC.plasma.MPHP<- 
  ggplot() + 
  geom_ribbon(data = MC.plasma.MPHP.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="steelblue1", alpha=0.3) +
  geom_ribbon(data = MC.plasma.MPHP.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="steelblue3", alpha = 0.5) +
  geom_line(data= MC.plasma.MPHP.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.plasma.MPHP.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = Eva_plasma_MPHP.obs, aes(x=Time, y= MPHP), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Eva_plasma_MPHP.obs, aes(x=Time,ymin= MPHP-SD, ymax = MPHP+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black") +
  ggtitle("MPHP (DPHP)")
MC.plasma.MPHP


MC.plasma.OH<- 
  ggplot() + 
  geom_ribbon(data = MC.plasma.OH.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="steelblue1", alpha=0.3) +
  geom_ribbon(data = MC.plasma.OH.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="steelblue3", alpha = 0.5) +
  geom_line(data= MC.plasma.OH.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.plasma.OH.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = Eva_plasma_OH.obs, aes(x=Time, y= OH), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Eva_plasma_OH.obs, aes(x=Time,ymin= OH-SD, ymax = OH+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black") +
  ggtitle("OH-MPHP (DPHP)")
MC.plasma.OH

MC.plasma.oxo<- 
  ggplot() + 
  geom_ribbon(data = MC.plasma.oxo.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="steelblue1", alpha=0.3) +
  geom_ribbon(data = MC.plasma.oxo.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="steelblue3", alpha = 0.5) +
  geom_line(data= MC.plasma.oxo.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "steelblue4") +
  geom_line(data= MC.plasma.oxo.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = Eva_plasma_oxo.obs, aes(x=Time, y= oxo), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Eva_plasma_oxo.obs, aes(x=Time,ymin= oxo-SD, ymax = oxo+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black") +
  ggtitle("oxo-MPHP (DPHP)")
MC.plasma.oxo


MC.urine.MPHP  <- 
  ggplot() + 
  geom_ribbon(data = MC.urine.MPHP.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.urine.MPHP.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.urine.MPHP.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "gold4") +
  geom_line(data= MC.urine.MPHP.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = Eva_urine_MPHP.obs, aes(x=Time, y= MPHP), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Eva_urine_MPHP.obs, aes(x=Time,ymin= MPHP-SD, ymax = MPHP+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("MPHP (DPHP)")
MC.urine.MPHP 


MC.urine.OH  <- 
  ggplot() + 
  geom_ribbon(data = MC.urine.OH.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.urine.OH.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.urine.OH.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "gold4") +
  geom_line(data= MC.urine.OH.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = Eva_urine_OH.obs, aes(x=Time, y= OH), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Eva_urine_OH.obs, aes(x=Time,ymin= OH-SD, ymax = OH+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black")+
  ggtitle("OH-MPHP (DPHP)")
MC.urine.OH 


MC.urine.oxo  <- 
  ggplot() + 
  geom_ribbon(data = MC.urine.oxo.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="gold", alpha=0.5) +
  geom_ribbon(data = MC.urine.oxo.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="gold3", alpha = 0.5) +
  geom_line(data= MC.urine.oxo.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "gold4") +
  geom_line(data= MC.urine.oxo.plot, aes(x = Time, y = median_est), 
            colour = "black") +
  geom_point(data = Eva_urine_oxo.obs, aes(x=Time, y= oxo), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  geom_errorbar(data = Eva_urine_oxo.obs, aes(x=Time,ymin= oxo-SD, ymax = oxo+SD), size = 0.8,width=0.5,
                position=position_dodge(0.05),colour="black") +
  ggtitle("oxo-MPHP (DPHP)")
MC.urine.oxo 

MC.plasma.DPHP  <- MC.plasma.DPHP+
  scale_y_continuous(limits = c(-50,320),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,48),expand = c(0,0))


MC.plasma.MPHP  <- MC.plasma.MPHP +
  scale_y_continuous(limits = c(-10,145),expand = c(0,0))+ 
  scale_x_continuous(limits = c(-0.5,29),expand = c(0,0))

MC.plasma.OH  <- MC.plasma.OH +
  scale_y_continuous(limits = c(-5,70),expand = c(0,0))+ 
  scale_x_continuous(limits = c(-0.5,29),expand = c(0,0))

MC.plasma.oxo  <- MC.plasma.oxo +
  scale_y_continuous(limits = c(-6,100),expand = c(0,0))+ 
  scale_x_continuous(limits = c(-0.2,29),expand = c(0,0))



MC.urine.MPHP    <- MC.urine.MPHP +
  scale_y_continuous(limits = c(0,0.15),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,48),expand = c(0,0))

MC.urine.OH    <- MC.urine.OH +
  scale_y_continuous(limits = c(0,5),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,48),expand = c(0,0))

MC.urine.oxo    <- MC.urine.oxo +
  scale_y_continuous(limits = c(0,5),expand = c(0,0))+ 
  scale_x_continuous(limits = c(0,48),expand = c(0,0))

## Modified the theme
## plot

MC.plasma.DPHP  = MC.plasma.DPHP + Theme.Fig + labs (x ="Time (hours)", y=expression(paste("Concentration in plasma (", mu,"g/L)")))
MC.plasma.MPHP  = MC.plasma.MPHP + Theme.Fig + labs (x ="Time (hours)", y=expression(paste("MPHP Concentration in plasma (", mu,"g/L)")))
MC.plasma.OH    = MC.plasma.OH + Theme.Fig + labs (x ="Time (hours)", y=expression(paste("Concentration in plasma (", mu,"g/L)")))
MC.plasma.oxo   = MC.plasma.oxo + Theme.Figb + labs (x ="Time (hours)", y=expression(paste("oxo Concentration in plasma (", mu,"g/L)")))

MC.urine.MPHP    = MC.urine.MPHP   + Theme.Fig + labs (x ="Time (hours)",  y="Cumulative excretion [%]")
MC.urine.OH      = MC.urine.OH     + Theme.Fig + labs (x ="Time (hours)",  y="Cumulative excretion [%]")
MC.urine.oxo     = MC.urine.oxo    + Theme.Fig + labs (x ="Time (hours)",  y="Cumulative excretion [%]")


#p1.r   = p1.r + Theme.Fig + labs (x ="Dose (mg/kg-day)", y="")


MC.plasma.DPHP
MC.plasma.MPHP
MC.plasma.OH
MC.plasma.oxo

MC.urine.MPHP 
MC.urine.OH 
MC.urine.oxo 


### Save figure #####

ggsave("MC.plasma.DPHP.tiff",scale = 1,
       plot = MC.plasma.DPHP,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("MC.plasma.MPHP.tiff",scale = 1,
       plot = MC.plasma.MPHP,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result",
       width = 20, height = 15, units = "cm",dpi=320)


ggsave("MC.plasma.OH.tiff",scale = 1,
       plot = MC.plasma.OH,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result",
       width = 20, height = 15, units = "cm",dpi=320)

ggsave("MC.plasma.oxo.tiff",scale = 1,
       plot = MC.plasma.oxo,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result",
       width = 20, height = 15, units = "cm",dpi=320)



ggsave("MC.urine.MPHP.tiff",scale = 1,
       plot = MC.urine.MPHP,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result",
       width = 20, height = 15, units = "cm",dpi=320)

ggsave("MC.urine.OH.tiff",scale = 1,
       plot = MC.urine.OH,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result",
       width = 20, height = 15, units = "cm",dpi=320)

ggsave("MC.urine.oxo.tiff",scale = 1,
       plot = MC.urine.oxo,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result",
       width = 20, height = 15, units = "cm",dpi=320)


grid.arrange(MC.DINCH.MHNCH, MC.plasma.oxo,nrow = 1)
ggsave("Fig.TIME_COURSE2A.tiff",scale = 1,
       plot = grid.arrange(MC.DINCH.MHNCH, MC.plasma.oxo,nrow = 1),
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot",
       width = 24, height = 12, units = "cm",dpi=320)

