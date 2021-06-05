## loading R packages
## Load libraries
library(mrgsolve) # Needed to run the main PBPK code
library(magrittr) # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)  # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(dplyr)    # Needed for the pipe %>% operator
library(reshape)  # melt function to reshape the table, reshape2 is version 2 of reshape. reshape is more stable when using melt function. We decided to use reshape for all.

rm(list=ls())
## Input PBPK model
humanPK.code    <- readRDS (file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/humanPK.RDS")

## Loading pbpk model 
mod <- mcode ("pk", humanPK.code)


## Loading human observed data
Inhal_2018  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2018_Inhalation.csv")
names(Inhal_2018)=c("Time", "Urine")
PlasmaDEHP_2012 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2012_Plasma_DEHP.csv")
names(PlasmaDEHP_2012)=c("Time", "CA_DEHP")
PlasmaMEHP_2012 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2012_Plasma_MEHP.csv")
names(PlasmaMEHP_2012)=c("Time", "CA")
UrineMEHP_2012  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2012_Urine_MEHP.csv")
names(UrineMEHP_2012)=c("Time", "Urine")
PlasmaMEHP_2004 <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2004_Plasma.csv")
names(PlasmaMEHP_2004)=c("Time", "CA")
UrineMEHP_2004  <- read.csv(file="C:/Users/Punkostrich/Dropbox/R/DEHP/2004_Urine.csv")
names(UrineMEHP_2004)=c("Time", "Urine")


## Loading MCMC data
Human.MCMC        <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/Human.MCMC.rds")

## loading the theta names
theta.names       <- readRDS(file = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/theta.names.rds")
which_sig         <- grep("sig", theta.names)

############################# prediction function ############################
############################################# Human data ################################################################


## Preiction function
pred <- function(pars,pred=FALSE) {
  
  ## Exposure scenario for oral exposue 2018
  
  tinterval              = 24                                        # Time interval
  TDoses                 = 1                                         # Dose times
  DOSEinhal2018          = 97.96                                     # Amount of oral dose
  ex2018               <- ev (ID = 1, amt= DOSEinhal2018, 
                                ii = tinterval, addl=TDoses-1, cmt="AC",tinf = 3, replicate = FALSE)
  
  
  ## Exposure scenario for oral exposue to 2 mg/kg
  
  DOSEoral2012         = 52800                   ## Amount of oral dose
  ex2012               <- ev (ID=1, amt= DOSEoral2012, 
                                ii=tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  
  ## Exposure scenario for oral exposue to 15 mg/kg
                          
  DOSEoral2004         = 48500                     ## mg, amount of oral dose                   
  ex2004               <- ev(ID=1,amt= DOSEoral2004, 
                               ii=tinterval, addl=TDoses-1,cmt="AST",replicate = FALSE)
  
  
 
  ## Set up the exposure time
  tsamp              = tgrid(0,tinterval*(TDoses-1)+24*2,0.01)   ## Siumuation 24*2 hours (2 days)
  tsamp2              = tgrid(0,tinterval*(TDoses-1)+24*2,0.1)   ## Siumuation 24*2 hours (2 days)
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  names(pars) <- names(pars)
  
  ## Exposure scenario for oral exposue to 4.2 mg/kg
  out.2018 <- 
    mod %>% 
    param(pars, BW = 60) %>%
    Req (Urine)%>%
    update(atol = 1E-10,maxsteps = 50000) %>%
    mrgsim_d (data = ex2018, tgrid = tsamp2)
  
  outdf.2018 <- cbind.data.frame(Time = out.2018$time, 
                                 Curine = (out.2018$Urine/DOSEinhal2018)*100)        ## %Dose, PFOS in urine
  
  ## Exposure scenario for oral exposue to 2 mg/kg
  out.2012 <- 
    mod %>% 
    param(pars, BW =82) %>%
    Req (Plasma, Plasma_MEHP, Urine)%>%
    update(atol = 1E-10,maxsteps = 50000) %>%
    mrgsim_d (data = ex2012, tgrid = tsamp)
  
  outdf.2012 <- cbind.data.frame(Time   = out.2012$time, 
                               CA_DEHP  = out.2012$Plasma,
                                 CA     = out.2012$Plasma_MEHP,
                                 Curine =(out.2012$Urine/DOSEoral2012)*100)        ## %Dose, PFOS in urine
  
  ## Exposure scenario for oral exposue to 15 mg/kg
  out.2004 <- 
    mod %>% 
    param(pars, BW = 75) %>%
    Req (Plasma_MEHP, Urine)%>%
    update(atol = 1E-10,maxsteps = 50000) %>%
    mrgsim_d (data = ex2004, tgrid = tsamp)
  
  outdf.2004 <- cbind.data.frame(Time   = out.2004$time, 
                                 CA     = out.2004$Plasma_MEHP,                            ## ug/ml,PFOS levels in plasma  
                                 Curine = (out.2004$Urine))        ## %Dose, PFOS in urine
    
  
  
  if (pred) return (list( "outdf.2018"  = outdf.2018, # if shows an error like this "　", just delete it.
                          "outdf.2012"  = outdf.2012,
                          "outdf.2004"  = outdf.2004))
    
  outdf.A1 = outdf.2018[which(outdf.2018$Time %in% Inhal_2018$Time), ]
  outdf.B1 = outdf.2012[which(outdf.2012$Time %in% PlasmaDEHP_2012$Time), ]
  outdf.B2 = outdf.2012[which(outdf.2012$Time %in% PlasmaMEHP_2012$Time), ]
  outdf.B3 = outdf.2012[which(outdf.2012$Time %in% UrineMEHP_2012$Time), ]
  outdf.C1 = outdf.2004[which(outdf.2004$Time %in% PlasmaMEHP_2004$Time), ]
  outdf.C2 = outdf.2004[which(outdf.2004$Time %in% UrineMEHP_2004$Time), ]
 
  
  return (list( "outdf.A1"  = outdf.A1, 
                "outdf.B1"  = outdf.B1,
                "outdf.B2"  = outdf.B2,
                "outdf.B3"  = outdf.B3,
                "outdf.C1"  = outdf.C1, 
                "outdf.C2"  = outdf.C2,
                "outdf.2018" = outdf.2018,
                "outdf.2004" = outdf.2004))

    
}


   

##################### Fig.4a Global fitting analysis ########################
## Human
pars.human        = Human.MCMC$bestpar[1:6]
predf             = pred (pars.human) # prediction results based on best-estimated parameters from MCMC
output.2018       <- predf$outdf.2018
output.2004       <- predf$outdf.2004

#output.2018 [ output.2018$Time == 3.8, ] 

# patch A1
time3.8 <- c(3.8, 2.3632360505) 
predf$outdf.A1 <- rbind(predf$outdf.A1[1:1,], time3.8, predf$outdf.A1[2:nrow(predf$outdf.A1),])
time6.6 <- c(6.6, 4.1743234571)
predf$outdf.A1 <- rbind(predf$outdf.A1[1:3,], time6.6, predf$outdf.A1[4:nrow(predf$outdf.A1),])
time12.6 <- c(12.6, 4.934801)
predf$outdf.A1 <- rbind(predf$outdf.A1[1:7,], time12.6, predf$outdf.A1[8:nrow(predf$outdf.A1),])
time14.2 <- c(14.2, 4.965789)
predf$outdf.A1 <- rbind(predf$outdf.A1[1:8,], time14.2, predf$outdf.A1[9:nrow(predf$outdf.A1),])
time17.9 <- c(17.9, 4.990687)
predf$outdf.A1 <- rbind(predf$outdf.A1[1:9,], time17.9, predf$outdf.A1[10:nrow(predf$outdf.A1),])
time21.9 <- c(21.9, 4.995905)
predf$outdf.A1 <- rbind(predf$outdf.A1[1:10,], time21.9, predf$outdf.A1[11:nrow(predf$outdf.A1),])

# patch C2
time10.2 <- c(10.2, 85.29180,	2375.854)
predf$outdf.C2 <- rbind(predf$outdf.C2[1:6,], time10.2, predf$outdf.C2[7:nrow(predf$outdf.C2),])


predf.A1          = cbind.data.frame(
  pre.value       = predf$outdf.A1$Curine,  
  obs.value       = Inhal_2018$Urine)

predf.B1          = cbind.data.frame(
  pre.value       = predf$outdf.B1$CA_DEHP[4:14],
  obs.value       = PlasmaDEHP_2012$CA[4:14]) # Exculde zero values

predf.B2          = cbind.data.frame(
  pre.value       = predf$outdf.B2$CA[1:14],
  obs.value       = PlasmaMEHP_2012$CA[1:14]) # Set the time points to match the prediction time points

predf.B3          = cbind.data.frame(
  pre.value       = predf$outdf.B3$Curine, 
  obs.value       = UrineMEHP_2012$Urine)

predf.C1          = cbind.data.frame(
  pre.value       = predf$outdf.C1$CA,
  obs.value       = PlasmaMEHP_2004$CA)

predf.C2          = cbind.data.frame(
  pre.value       = predf$outdf.C2$Curine/48500*100, 
  obs.value       = UrineMEHP_2004$Urine/48500*100)


ep_urine          =   rbind(predf.A1,predf.B3,predf.C2) # combination of many sets of predicted vs observed data
ep_plasma         =   rbind(predf.B1,predf.B2,predf.C1) # combination of many sets of predicted vs observed data


ep_urine$tissue  = rep("Urine",length(ep_urine $obs.value)) # Create a Species column, the number of rows is the same as the number of obs.value 
ep_plasma$tissue  = rep("Plasma",length(ep_plasma $obs.value))

plot(x=log10(ep_urine$pre.value),y = log10(ep_urine$obs.value))
abline(lm(log10(ep_urine$obs.value) ~ log10(ep_urine$pre.value)))


plot(x=log10(ep_plasma$pre.value),y = log10(ep_plasma$obs.value))
abline(lm(log10(ep_plasma$obs.value) ~ log10(ep_plasma$pre.value)))

############## Data of global evaluation of model fit for Figure 4a ################
## All data
ep.all =rbind(ep_urine,ep_plasma)
ep.all$log.obs = log(ep.all$obs.value,10)
ep.all$log.pre = log(ep.all$pre.value,10)
fit <- lm(log.obs ~ log.pre, data=ep.all)                   # Fit the model
ep.all$residuals = residuals(fit)                           # Save the residual values
ep.all$predicted = predict(fit)                             # Save the predicted values
ep.all$OPratio = ep.all$pre.value/ep.all$obs.value          # Estimated the predicted-to-observed ratio



############## plot #####################
p <- 
  ggplot(ep.all, aes(log.obs, log.pre, fill = tissue)) + 
  geom_segment (aes(xend    = log.obs,                   # connect the actual data points with their corresponding predicted value    
                    yend    = predicted),
                alpha   = .2) +
  geom_point   (aes(shape   = tissue, colour = tissue), size = 6)  +
  scale_colour_manual(labels = c(expression(paste("Plasma (", mu,"g/L)")), "Urine (%)"),values=c('#2E9FDF', '#E7B800')) +
  scale_shape_manual(labels = c(expression(paste("Plasma (", mu,"g/L)")), "Urine (%)"), values = c(19, 17)) +
  scale_fill_discrete(labels = c(expression(paste("Plasma (", mu,"g/L)")),"Urine (%)")) +
  #scale_fill_manual(values=c('blue', 'green')) + # legend name
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="dodgerblue4", size = 1) +
  geom_abline (intercept = log10(10), linetype = "dashed",
               slope     = 1,
               color     ="dodgerblue4", size = 1) +
  geom_abline (intercept = log10(0.1), linetype = "dashed", 
               slope     = 1,
               color     ="dodgerblue4", size = 1) +
  geom_abline (intercept = log10(3), linetype = "dotted",
               slope     = 1,
               color     ="dodgerblue4", size = 1) +
  geom_abline (intercept = log10(0.333), linetype = "dotted", 
               slope     = 1,
               color     ="dodgerblue4", size = 1) 



p<-
  p +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,4), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,4),labels = scales::math_format(10^.x))

windowsFonts(Times=windowsFont("Times New Roman"))

p1 <- p + 
    
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 20, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 20, colour = "black", face = "bold"),   # label of axes
    legend.title            = element_blank(),
    legend.justification    =  c("right", "bottom"),
    legend.position         = c(0.98, .05),
    legend.text             = element_text (size = 20),
    legend.text.align       = 0) + 

    labs (x = "Observed value", y = "Predicted value")
  #labs (x = expression(paste("Observed value (", mu,"g/L)")), y = expression(paste("Predicted value (", mu,"g/L)")))

p1


ggsave("Fig.5a.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result",
       width = 15, height = 15, units = "cm",dpi=320)


##### count #########


plasma = ep.all[which(ep.all$tissue == 'Plasma'), ]
plasma_ratio = plasma$OPratio
i_plasma = length(plasma_ratio)
j_plasma_3 = length(which(plasma_ratio >= 0.3 & plasma_ratio <= 3))
plasma_3 = j_plasma_3/i_plasma *100
j_plasma_10 = length(which(plasma_ratio >= 0.1 & plasma_ratio <= 10))
plasma_10 = j_plasma_10/i_plasma *100
p = c(plasma_3,plasma_10)

urine = ep.all[which(ep.all$tissue == 'Urine'), ]
urine_ratio = urine$OPratio
i_urine = length(urine_ratio)
j_urine_3 = length(which(urine_ratio >= 0.3 & urine_ratio <= 3))
urine_3 = j_urine_3/i_urine *100
j_urine_10 = length(which(urine_ratio >= 0.1 & urine_ratio <= 10))
urine_10 = j_urine_10/i_urine *100
u = c(urine_3,urine_10)


############# RMSE calculation ##############
#############################################
library("Metrics")

############# calculate RMSE ############
re <- rmse(ep.all$obs.value, ep.all$pre.value)
re
log10(re)

#### TS calculation####################
 residue <- ep.all$pre.value - ep.all$obs.value         # prediction - actual value
 N <- length(residue)
 abs_residue <- abs(ep.all$pre.value - ep.all$obs.value)

 s_res <- sum(residue)
 s_absres <- sum(abs_residue)
 
 TS <- s_res/s_absres * N
 TS
 
