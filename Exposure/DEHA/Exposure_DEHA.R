
library(truncnorm) 
library(ggridges)     # Used to create Figure
library(ggplot2)   # ggplot is the basic package for creating plots.
#install.packages("tkrplot")
library(tkrplot)
#install.packages("rriskDistributions")
library(rriskDistributions)
library(EnvStats)

##### DEHA
 ############################################## 1  ##############################################
 ################################ Exposure scenario #Inhalation #################################
 ## parameters
 sn <- 10000
 volume <- 0.453592 * 1E6   # kg
 dens_DEHA   <- 0.922 * 1E3 #kg/m3
 
 population <- httkpop_generate(method = 'virtual individuals',nsamp = N,
                                agelim_years=c(20,65))
 
 SA <- as.numeric((population$weight * population$height)^0.5/60)   #skin surface area # unshared (Us)
 BW <- as.numeric(population$weight)  # Body weight  # unshared (Us)
 
 ##### percentage in plastic, same as DEHP
 SD_perc <- 0.3*0.5
 perc <-rtruncnorm (n=10000, a = 0.1, 
                   b = qnorm(0.95, mean = 0.3, sd =  SD_perc), 
                   mean = 0.3, sd =  SD_perc) 
 quantile(perc, c(0.05,.5, 0.75, 0.9, 0.95)) 
 
 life <- runif(n=10000, min = 5, max = 30)     # plastic service life
 dens <- runif(n=10000, min = 1.1, max = 3)   #kg/m3
 
 ######## household size ########
 SD_household <- 2.65*1
 household<-rtruncnorm (n=10000, a = 1, 
                    b = qnorm(0.95, mean = 2.65, sd =  SD_household), 
                    mean = 2.65, sd =  SD_household) 
 quantile(household, c(0.05,.5, 0.75, 0.9, 0.95)) 
 
 A_p <- volume * life / perc / dens / 328 /10000  #328 million people in US
 A <- volume * life / perc / dens / 328 /10000 * household
 quantile(A, c(0.05,.5, 0.75, 0.9, 0.95))       #household emission area
 quantile(A_p, c(0.05,.5, 0.75, 0.9, 0.95)) 
 
 y0  = min(10.4,10.4 * perc / dens_DEHA / (perc/dens_DEHA + (1-perc)/1.5E3) * 3.7) ## ug/m3, vapor pressure 


 ## convective mass transfer coefficient log-normal
 dh <- c(3.6, 7.2)
 Percdh <- c(0.5, 0.95)
 b <- get.lnorm.par(p= Percdh,q=dh)
 h <- rlnorm(n=sn, meanlog = 1.2809338   , sdlog = 0.4214028)
 quantile(h, c(.5, 0.75, 0.9, 0.95))  # check ## m/h
 
 ## indoor ventilation rate Q normal
 SD_Q <- (88 - 221) / qnorm(0.1)
 Q <-rtruncnorm (n=10000, a = qnorm(0.05, mean = 221, sd = SD_Q), 
                          b = qnorm(0.95, mean = 221, sd = SD_Q), 
                          mean = 221, sd = SD_Q)              ## m3/h
 
 ## inhalation rate normal
 SD_InhR <- (21 - 15) / qnorm(0.95)
 InhR <-rtruncnorm (n=10000, a = qnorm(0.05, mean = 15, sd = SD_InhR ), 
                          b = qnorm(0.95, mean = 15, sd = SD_InhR ), 
                          mean = 15, sd = SD_InhR )              ## m3/day

 ## boday weight normal
# SD_BW <- (118 - 80.8) / qnorm(0.95)
# BW <-rtruncnorm (n=10000, a = qnorm(0.05, mean = 80.8, sd = SD_BW ), 
#                          b = qnorm(0.95, mean = 80.8, sd = SD_BW ), 
#                          mean = 80.8, sd = SD_BW )              ## kg
 ## tsp log-normal
 SD_TSP <- 20*5
 TSP <-rtruncnorm (n=10000, a = 0.1, 
                          b = qnorm(0.95, mean = 20, sd = SD_TSP), 
                          mean = 20, sd = SD_TSP)             
 dTSP <- c(20, 100)
 PercdTSP <- c(0.5, 0.95)
 c <- get.lnorm.par(p= PercdTSP,q=dTSP)
 TSP <-  rlnormTrunc(n=sn, meanlog = 2.9957515 , sdlog = 0.9785187, max = 1000) ## ug/m3
 quantile(TSP, c(.5, 0.75, 0.9, 0.95))  

 ## exposure duration percentage assume 80% indoor
 ED <- 0.9  ## -

 ## particle-air partition coefficient normal
 SD_kp <- 0.02*0.35
 Kp <-rtruncnorm (n=10000, a = qnorm(0.05, mean = 0.02, sd = SD_kp), 
                 b = qnorm(0.95, mean = 0.02, sd = SD_kp), 
                 mean = 0.02, sd = SD_kp)              ## m3/ug
 quantile(Kp, c(0.05,.5, 0.75, 0.9, 0.95)) 

 ## Calculate DEHP conc in air 
 Cair <- h * y0 * A/(h * A + (1 + Kp * TSP) * Q)    ## ug/m3
 quantile(Cair*1000, c(0.05,.5, 0.75, 0.9, 0.95)) 
 ## inhalation exposure
 IE <- (Cair * InhR * ED + Cair * Kp * TSP * InhR * ED)/BW   ## ug/day/kg

 
 
 ######################################## 2 ####################################
 ################################ Exposure scenario #Oral ######################
 #diet oral
 ##kpf of DEHA
 kpf <- runif(n=10000, min = 1410, max = 5475) 
   
 ## concentration in food contacting material Cp (%) log-normal
 dCp <- c(0.2, 0.5)    #3%
 PercdCp <- c(0.5, 0.95)   #50%
 c <- get.lnorm.par(p= PercdCp,q=dCp)
 Cp_ori <- rlnormTrunc(n=sn, meanlog = -1.6093735, sdlog = 0.5571095, max = 0.7) ## g/g
 max(Cp_ori)
 quantile(Cp_ori, c(.05, .1,.25,.5, 0.75, 0.9, 0.95,1))  
 Cp <- Cp_ori * 1E6   ## ug/g
 quantile(Cp, c(.05, .1,.25,.5, 0.75, 0.9, 0.95,1))  
 
 ## meat consumption (g/kg-day) normal
 SD_meat <- (4.1 - 1.8) / qnorm(0.95)
 meat <-rtruncnorm (n=10000, a = qnorm(0.05, mean = 1.8, sd = SD_meat ), 
                  b = qnorm(0.95, mean = 4.1, sd = SD_meat ), 
                  mean = 1.8, sd = SD_meat )              ## g/kg-day
 
 ## cheese consmption normal
 SD_dairy <- 16 * 1
 dairy <- rtruncnorm (n=10000, a = 1, 
                    b = qnorm(0.95, mean = 16 , sd = SD_dairy ), 
                    mean = 16 , sd = SD_dairy )              ## g/d
 
 ## fat consumption normal
 SD_fat <- (2.3 - 1.2) / qnorm(0.95)
 fat <- rtruncnorm (n=10000, a = 1.2, 
                      b = qnorm(0.95, mean = 1.2 , sd = SD_fat ), 
                      mean = 1.2 , sd = SD_fat )              ## g/kg-d
 
 Diet <- Cp/kpf * (meat*BW + dairy + fat*BW)      ##ug/day
 quantile(Diet/BW, c(.5, 0.75, 0.9, 0.95, 0.975))  # check
 
 
 #dust oral
 IngR  <- runif(n=10000, min = 30, max = 70)  ## mg/d

 #SD_Kdust <- 2.5E-3 * 0.35
 #Kdust <-rtruncnorm (n=10000, a = qnorm(0.05, mean = 2.5E-3, sd = SD_Kdust), 
 #                         b = qnorm(0.95, mean = 2.5E-3, sd = SD_Kdust ), 
 #                         mean = 2.5E-3, sd = SD_Kdust)              ## m3/mg
 
 Kdust <- Kp/8.32            ## m3/ug 


 OE <- (Cair * Kdust * IngR)/BW  + Diet/BW ## ug/day/kg
 quantile((Cair * Kdust * IngR)/BW, c(.5, 0.75, 0.9, 0.95)) 
 quantile(OE, c(.5, 0.75, 0.9, 0.95))
 

 ############################## 3 ###################################
 ##################### Exposure scenario #Dermal ####################

 SD_kpg <- 1.77 * 5
 kpg <- rtruncnorm (n=10000, a = 0.1e-03, 
                          b = qnorm(0.95, mean = 1.91, sd = SD_kpg ), 
                          mean = 1.91 , sd = SD_kpg )              ## m/d


 #SD_SA <- (2.5 - 2.1) / qnorm(0.95)
# SA <-rtruncnorm (n=10000, a = qnorm(0.05, mean = 2.1, sd = SD_SA), 
#                          b = qnorm(0.95, mean = 2.1, sd = SD_SA), 
#                          mean = 2.1, sd = SD_SA)              ## m2

 SD_Fr <- 0.24 * 1
 Fr  <- rtruncnorm (n=10000, a = 0.11, 
                    b = qnorm(0.95, mean =0.24, sd = SD_Fr ), 
                    mean = 0.24 , sd = SD_Fr ) ## - fraction exposed to air
 quantile(Fr, c(0.05,.5, 0.75, 0.9, 0.95))

 SE <- (Cair * kpg * SA * Fr * ED)/BW  ## ug/day/kg


 ## Summary results
 ## Inhalation
 df.inhal<-as.data.frame(IE)
 df.inhal$DEHA <- c("Inhalation") 
 names(df.inhal)=c("Value","DEHA")

 ## Oral
 df.oral<-as.data.frame(OE)
 df.oral$DEHA <- c("Oral") 
 names(df.oral)=c("Value","DEHA")

 ## Dermal
 df.dermal<-as.data.frame(SE)
 df.dermal$DEHA <- c("Dermal") 
 names(df.dermal)=c("Value","DEHA")

 ################## summary #####################
 df <- rbind.data.frame (df.inhal,df.oral,df.dermal)
 df$log.value <- log10(df$Value)
 df$ggjoy = c("A")
 ######### end of exposure estimates ############

##################################################################
########## Figure: Exposure estimate distributions  ##############
##################################################################

windowsFonts(Times=windowsFont("Times New Roman")) # Abbreviate the font Times New Roman as Times

df$ggjoy = c("A") # Create a new column ggjoy in df table. Without this, the plot p1 will have four layers for four species overlapped in each panel, which does not look good.
# With this, y axis is A for all species in one layer in one panel.

p1 = 
  ggplot (df, aes(x = as.numeric(log.value), y = as.factor(ggjoy),fill = DEHA)) + # fill, fill one color for each species
  geom_density_ridges (scale = 28, size = 1, rel_min_height = 0.001, alpha = 0.6) + # over size of the entire ggridge plot
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0), name = expression(paste("Log [Daily Uptake] (", mu,"g/kg/d)"))) +
  scale_fill_manual(values = c("grey", "#2E9FDF","#7CAE00")) + 
  scale_colour_manual(values = c("grey", "#2E9FDF","#7CAE00"))

p1=
  p1 + theme_ridges()+ # change the theme of the plot p1 (incorporate p1 into the theme)
  theme (
    plot.background         = element_blank(),
    text                    = element_text (family = "Times",face="bold",size = 24),
    #panel.border            = element_rect (colour = "black", fill=NA, size=2),
    #panel.background        = element_blank(),
    panel.background        = element_rect (fill   = "#f7f7f7"),
    #panel.grid.major.y      = element_line (size   = 0.5, colour = "grey"),
    #panel.grid.minor.y      = element_blank(),
    panel.grid.major.x      = element_blank(), 
    #panel.grid.minor.x      = element_blank(), 
    axis.text               = element_text (size   = 24, colour = "black",face="bold"),  
    axis.title              = element_text ( size   = 24, colour = "black", face = "bold"),   # label of axes
    axis.ticks.x            = element_line (size   = 1, colour = "black"), 
    axis.ticks.y            = element_blank(), 
    axis.text.y             = element_blank(),
    strip.text              = element_text (size   = 24),
    axis.title.x            = element_text (hjust = 0.5, face="bold", size = 24),
    legend.title            = element_text (size   = 24, face="bold"),
    legend.justification    =  c("left", "top"),,
    legend.position         = c(.05, .95),
    legend.text             = element_text (size = 24,face="bold")) + 
    labs (x = "", y = "") 

p1

ggsave("Exposure_updatedA.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result",
       width = 25, height = 20, units = "cm",dpi=320)



###### get mean & percentile ##########
## Inhalation ##
mean_inhal = mean(IE)
quantile(IE, c(.05,0.5, .95))

## Oral ###
mean_oral = mean(OE)
quantile(OE, c(.05,0.5, .95))

## Dermal ##
mean_dermal = mean(SE)
quantile(SE, c(.05,0.5, .95))

## Overall ##
all = IE + OE + SE
quantile(all, c(.05,0.5, .95))

IE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/IE.rds')
OE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/OE.rds')
SE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/SE.rds')

all = IE + OE + SE
quantile(all, c(.05,0.5, .95))

saveRDS(IE,file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/IE.rds')
saveRDS(OE,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/OE.rds')
saveRDS(SE,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/SE.rds')
saveRDS(BW,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/BW.rds')
saveRDS(df,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/df.rds')
saveRDS(all,file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/all.rds')


IE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/IE.rds')
OE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/OE.rds')
SE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/SE.rds')
all = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/all.rds')

quantile(all, c(.05,0.5, .95))
