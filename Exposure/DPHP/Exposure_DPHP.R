
library(truncnorm) 
library(ggridges)     # Used to create Figure
library(ggplot2)   # ggplot is the basic package for creating plots.
#install.packages("tkrplot")
library(tkrplot)
#install.packages("rriskDistributions")
library(rriskDistributions)
library(EnvStats)
library(httk)
library(mc2d)
library(fitdistrplus)

rm(list = ls())
M <- 100    #M sets 
N <- 10000   #N exposure

set.seed(123)
############################################# outer loop ##############################################
## production volume #(Sc, shared among members of the entire cohort)
SD_volume <- 0.3*  0.453592 * 2.12e8    # 30% CV
volume.sample <-rtruncnorm (n=M, a = qnorm(0.05, mean =  0.453592 * 2.12e8    , sd = SD_volume  ), 
                            b = qnorm(0.95, mean = 0.453592 * 2.12e8 , sd = SD_volume ), 
                            mean = 0.453592 * 2.12e8 , sd = SD_volume )             # kg 
quantile(volume.sample, c(0,0.05,.5, 0.75, 0.9, 0.95,1)) 

##### percentage in plastic (UcSc)
perc_mean.sample <-  runif(n=M, min = 0.3, max = 0.6)     # plastic percentage (ScUc)
perc_CV.sample   <-  runif(n=M, min = 0.3, max = 1)     # plastic percentage, (ScUc)

life_mean.sample <- runif(n=M, min = 5, max = 30)     # plastic service life, (ScUc)
life_CV.sample <-  runif(n=M, min = 0.3, max = 1)     # plastic service life, (ScUc)

dens.sample <- runif(n=M, min = 1.1, max = 3)   #kg/m2 (Sc, shared among members of the entire cohort)
dens_DPHP   <- 0.965 * 1E3 #kg/m3

## particle-air partition coefficient normal  #(Sc, shared among members of the entire cohort)
SD_kp <- 0.28 * 0.85
Kp.sample <-rtruncnorm (n=M, a = 0.01, 
                 b = qnorm(0.95, mean = 0.28, sd = SD_kp), 
                 mean = 0.28, sd = SD_kp)              ## m3/ug
quantile(Kp.sample, c(0.0001,0.05,.5, 0.75, 0.9, 0.95)) 



## vapor from product #(Sc, shared among members of the entire cohort)
#SD_y0  <- 0.3 * 0.667  # 30% CV
#y0.sample  <- rtruncnorm (n=M, a = 0.5, 
#                          b = qnorm(0.95, mean = 0.667 , sd = SD_y0), 
#                          mean = 0.667, sd = SD_y0)             ## ug/m3 

##kpf of DPHP  #(ScUc, highly uncertain)
kpf_mean.sample <- runif(n=M,  min = 793917, max = 5487500)   #40-50% etoh content 
kpf_CV.sample <-  runif(n=M, min = 0.3, max = 1)     #  (ScUc)


## concentration in food contacting material Cp (%) log-normal # # (ScUc)
Cp_mean.sample <-  runif(n=M,  min = 0.01, max = 0.2) * 1E6   #40-50% etoh content 
Cp_CV.sample   <-  runif(n=M, min = 0.3, max = 1)     #  (ScUc)

#dust oral # unshared (ScUc)
IngR_mean.sample  <- runif(n=N, min = 30, max = 70)  ## mg/d
IngR_CV.sample   <-  runif(n=M, min = 0.3, max = 1)     # plastic service life, (ScUc)


occurance.sample <- runif(n=M, min = 0.03, max = 0.3) 

Sc.sample <- cbind(volume.sample, perc_mean.sample,perc_CV.sample,life_mean.sample, life_CV.sample ,  occurance.sample ,
                   dens.sample, Kp.sample, kpf_mean.sample, kpf_CV.sample , 
                   perc_mean.sample, perc_CV.sample, Cp_mean.sample, Cp_CV.sample,
                   IngR_mean.sample, IngR_CV.sample)



############################################# inner loop ##############################################
df       = matrix(0, N, 1)
df.IE    = matrix(0, N, 1)
df.OE    = matrix(0, N, 1)
df.SE    = matrix(0, N, 1)

dim(df)
######################################## creat visual population ######################################
for (i in 1:M){
  
  fixed = Sc.sample[i,]
  #print(fixed)
  
  volume = fixed['volume.sample']
  perc   = fixed['perc.sample']
  life_mean = fixed['life_mean.sample'] 
  life_CV   = fixed['life_CV.sample'] 
  dens = fixed['dens.sample'] 
  Kp   = fixed['Kp.sample']
  #y0  = fixed['y0.sample'] 
  kpf_mean = fixed['kpf_mean.sample']
  kpf_CV = fixed['kpf_CV.sample']
  Cp  = fixed['Cp.sample']
  perc_mean = fixed['perc_mean.sample']
  perc_CV   = fixed['perc_CV.sample']
  Cp_mean   = fixed['Cp_mean.sample']
  Cp_CV     = fixed['Cp_CV.sample']
  occurance = fixed['occurance.sample']
  IngR_mean = fixed['IngR_mean.sample']
  IngR_CV = fixed['IngR_CV.sample']
  
  
  set.seed(1)
  # httkpop_generate(method='direct resampling', nsamp=100)
  #Generate a population using the virtual-individuals method,
  #including 80 females and 20 males,
  #including only ages 20-65,
  #including non-obese individuals
  population <- httkpop_generate(method = 'virtual individuals',nsamp = N,
                                 agelim_years=c(20,65))
  
  BSA <- as.numeric((population$weight * population$height)^0.5/60)   #skin surface area # unshared (Us)
  BW <- as.numeric(population$weight)  # Body weight  # unshared (Us)
  #hist(as.numeric(unlist(BSA)))
  #hist(as.numeric(unlist(population$weight)))
  #mean(as.numeric(unlist(population$weight)))
  #mean(as.numeric(unlist(BSA)))
  ############################################## 1  ##############################################
  ################################ Exposure scenario #Inhalation #################################
  ## parameters
  life <- rtruncnorm (n=N, a = 1, 
                      b = 70, 
                      mean = life_mean, sd =  life_CV * life_mean)  # ScUc
  
  kpf <- rtruncnorm (n=N, a = 1E2, 
                     b = qnorm(0.95, mean = kpf_mean, sd =  kpf_CV * kpf_mean), 
                     mean = kpf_mean, sd =  kpf_CV * kpf_mean)  # ScUc
  
  perc <- rtruncnorm (n=N, a = 0.01, 
                     b = 0.8, 
                     mean = perc_mean, sd =  perc_CV * perc_mean)  # ScUc
  
  Cp <- rtruncnorm (n=N, a = 0.01, 
                      b = 0.8*1E6, 
                      mean = Cp_mean, sd =  Cp_CV * Cp_mean)  # ScUc
  
  IngR <- rtruncnorm (n=N, a = 0.01, 
                      b = 1000, 
                      mean = IngR_mean, sd =  IngR_CV * IngR_mean)  # ScUc
  
  Cp_ori = Cp / 1E6
  y0  = min(0.667, 0.667 * perc / dens_DPHP / (perc/dens_DPHP + (1-perc)/1.5E3) * 3.7) ## ug/m3,
  
  ######## household size ######## # unshared (Us)
  SD_household <- 2.65*1 
  household<-rtruncnorm (n=N, a = 1, 
                         b = qnorm(0.95, mean = 2.65, sd =  SD_household), 
                         mean = 2.65, sd =  SD_household)  # unshared (Us)
  #quantile(household, c(0.05,.5, 0.75, 0.9, 0.95,1)) 
  
  
  ## convective mass transfer coefficient log-normal # unshared (Us)
  #dh <- c(3.6, 7.2)
  #Percdh <- c(0.5, 0.95)
  #b <- get.lnorm.par(p= Percdh,q=dh)
  h <- rlnorm(n=N, meanlog = 1.2809338   , sdlog = 0.4214028) # # unshared (Us)
  #quantile(h, c(.5, 0.75, 0.9, 0.95,1))  # check ## m/h
  
  ## indoor ventilation rate Q normal  # unshared (Us)
  SD_Q <- (88 - 221) / qnorm(0.1)
  Q <-rtruncnorm (n=N, a = qnorm(0.05, mean = 221, sd = SD_Q), 
                  b = qnorm(0.95, mean = 221, sd = SD_Q), 
                  mean = 221, sd = SD_Q)              ## m3/h  # unshared (Us)
  
  ## inhalation rate normal  # unshared (Us)
  SD_InhR <- (21 - 15) / qnorm(0.95)
  InhR <-rtruncnorm (n=N, a = qnorm(0.05, mean = 15, sd = SD_InhR ), 
                     b = qnorm(0.95, mean = 15, sd = SD_InhR ), 
                     mean = 15, sd = SD_InhR )              ## m3/day  # unshared (Us)
  
  ## tsp log-normal  # unshared (Us)
  #dTSP <- c(20, 100)
  #PercdTSP <- c(0.5, 0.95)
  #c <- get.lnorm.par(p= PercdTSP,q=dTSP)
  TSP <-  rlnormTrunc(n=N, meanlog = 2.9957515 , sdlog = 0.9785187, max = 1000) ## ug/m3
  #quantile(TSP, c(.5, 0.75, 0.9, 0.95,1))  
  
  ## exposure duration percentage assume 90% indoor  # unshared (Us)
  SD_ED <- 0.3 * 0.9
  ED <-rtruncnorm (n=N, a = 0.4, 
                   b = 1, 
                   mean = 0.9, sd = SD_ED )   # unshared (Us)
  #quantile(ED, c(.5, 0.75, 0.9, 0.95,1))
  #hist(ED)
  
  
  ## Calculate DPHP conc in air 
  
  A_p <- volume * life / perc / dens / 328 /1000000
  A <- volume * life / perc / dens / 328 /1000000 * household
  #quantile(A, c(0.05,.5, 0.75, 0.9, 0.95)) 
  #quantile(A_p, c(0.05,.5, 0.75, 0.9, 0.95)) 
  
  Cair <- h * y0 * A/(h * A + (1 + Kp * TSP) * Q)    ## ug/m3
  ## inhalation exposure
  IE <- (Cair * InhR * ED + Cair * Kp * TSP * InhR * ED)/BW   ## ug/day/kg
  # print(min(IE))
  
  
  
  ######################################## 2 ####################################
  ################################ Exposure scenario #Oral ######################
  #diet oral
  
  
  
  ## meat consumption (g/kg-day) normal # unshared (Us)
  SD_meat <- (4.1 - 1.8) / qnorm(0.95)
  meat <-rtruncnorm (n=N, a = qnorm(0.05, mean = 1.8, sd = SD_meat ), 
                     b = qnorm(0.95, mean = 4.1, sd = SD_meat ), 
                     mean = 1.8, sd = SD_meat )              ## g/kg-day
  
  ## cheese consmption normal # unshared (Us)
  SD_dairy <- 16 * 1
  dairy <- rtruncnorm (n=N, a = 1, 
                       b = qnorm(0.95, mean = 16 , sd = SD_dairy ), 
                       mean = 16 , sd = SD_dairy )              ## g/d
  
  ## fat consumption normal # unshared (Us)
  SD_fat <- (2.3 - 1.2) / qnorm(0.95)
  fat <- rtruncnorm (n=N, a = 1.2, 
                     b = qnorm(0.95, mean = 1.2 , sd = SD_fat ), 
                     mean = 1.2 , sd = SD_fat )              ## g/kg-d
  
  Diet <- Cp/kpf * occurance * (meat*BW + dairy + fat*BW)      ##ug/day
  #quantile(Diet/BW, c(.5, 0.75, 0.9, 0.95, 0.975))  # check
  
  
  #dust oral # unshared (Us)
  IngR  <- runif(n=N, min = 30, max = 70)  ## mg/d
  
  Kdust <- Kp/8.32            ## m3/ug #(Sc, shared among members of the entire cohort)
  #quantile(Kdust, c(0.05,.5, 0.75, 0.9, 0.95, 0.975))  # check
  
  OE <- (Cair * Kdust * IngR)/BW  + Diet/BW ## ug/day/kg
  #quantile((Cair * Kdust * IngR)/BW, c(.5, 0.75, 0.9, 0.95)) 
  #quantile(OE, c(.5, 0.75, 0.9, 0.95))
  
  
  ############################## 3 ###################################
  ##################### Exposure scenario #Dermal ####################
  ## dermal absorption coefficient # unshared (Us)
  SD_kpg <- 7.55 * 5
  kpg <- rtruncnorm (n=N, a = 0.1e-03, 
                     b = qnorm(0.95, mean = 7.55, sd = SD_kpg ), 
                     mean = 7.55, sd = SD_kpg )              ## m/d
  
  
  ## fraction of skin exposed  # unshared (Us)
  SD_Fr <- 0.24 * 1
  Fr  <- rtruncnorm (n=N, a = 0.11, 
                     b = qnorm(0.95, mean =0.24, sd = SD_Fr ), 
                     mean = 0.24 , sd = SD_Fr ) ## - # unshared (Us)
  #quantile(Fr, c(0.05,.5, 0.75, 0.9, 0.95))
  
  SE <- (Cair * kpg * BSA * Fr * ED)/BW  ## ug/day/kg
  
  
  
  all = IE + OE + SE
  
  df       = cbind (df, all)
  df.IE    = cbind (df.IE, IE)
  df.OE    = cbind (df.OE, OE)
  df.SE    = cbind (df.SE, SE)
  
  
  cat("loop= ", i , "\n")
  
}


## Summary results
## Inhalation
df.inhal<-as.data.frame(df.IE)
df.inhal = subset(df.inhal, select = -c(V1))
df.inhal$DPHP <- c("Inhalation") 
#names(df.inhal)=c("Value","DPHP")

## Oral
df.oral<-as.data.frame(df.OE)
df.oral = subset(df.oral, select = -c(V1) )
df.oral$DPHP <- c("Oral") 
#names(df.oral)=c("Value","DPHP")

## Dermal
df.dermal<-as.data.frame(df.SE)
df.dermal = subset(df.dermal, select = -c(V1) )
df.dermal$DPHP <- c("Dermal") 
#names(df.dermal)=c("Value","DPHP")

## all
df.all <-as.data.frame(df)
df.all = as.data.frame(subset(df.all, select = -c(V1) ))

# get measured values from literature
x.median <- c(0.04)
x.95 <- c(0.24)
survey.median <- quantile(x.median, c(.95))
survey.95 <- quantile(x.95, c(.95))

SD_survey <- (survey.95 - survey.median) / qnorm(0.95)
survey <- rlnormTrunc(n=N, meanlog = survey.median, sdlog =  SD_survey, min = 1e-10, max = 5000)

lst.ks <- lapply(1:ncol(df.all), function(i)
 # z.test(log(df.all[, i]), sigma.x=1, log(survey), sigma.y=1, conf.level=0.90)$statistic)
  #round(bhattacharyya.dist(log(df.all[, i]), log(survey) ,diag(N),diag(N)),digits=2))
  JSD(rbind(df.all[, i], survey), unit = 'log2'))
lst.ks <- as.data.frame(lst.ks)
#A = z.test(log(df.all[, 2]), sigma.x=1, log(survey), sigma.y=1, conf.level=0.90)
# get residues
df.2.3 <- as.data.frame(abs(lst.ks))

# order by residues, for quality check only
y = df.2.3[,order(df.2.3[nrow(df.2.3),], decreasing = FALSE)]
y

df.2.3$DPHP <- c("quantile")
df.all$DPHP <- c("all") 

library(data.table) #data.table_1.9.5
df.2 <- rbindlist(list(df.inhal, df.oral, df.dermal, df.all, df.2.3),use.names=FALSE )
df.2R <- as.data.frame(df.2)
df.2A = as.data.frame(subset(df.2R, select = -c(DPHP) ))

df.ks <- df.2A[,order(df.2A[(nrow(df.2A)),], decreasing = FALSE)]       # order column based on p value of kstest
# select first 5 columns
df.f <- df.ks[, 1]
df.f2 <- cbind.data.frame(df.f, df.2R$DPHP)
names(df.f2)[length(names(df.f2))]<-"DPHP" 

# extrace subset data frame for inhal, dermal,  oral, all sepearately
df.inhal <- df.f2[df.f2$DPHP == 'Inhalation', ]
df.dermal <- df.f2[df.f2$DPHP == 'Dermal', ]
df.oral <- df.f2[df.f2$DPHP == 'Oral', ]
df.all <- df.f2[df.f2$DPHP == 'all', ]


# mix all data within each category
df.inhal.A = as.data.frame(subset(df.inhal, select = -c(DPHP) ))
df.inhal <- data.frame(x=unlist(df.inhal.A))
df.inhal$DPHP <- c("Inhalation") 
names(df.inhal)=c("Value","DPHP")

df.dermal.A = as.data.frame(subset(df.dermal, select = -c(DPHP) ))
df.dermal <- data.frame(x=unlist(df.dermal.A))
df.dermal$DPHP <- c("Dermal") 
names(df.dermal)=c("Value","DPHP")


df.oral.A = as.data.frame(subset(df.oral, select = -c(DPHP) ))
df.oral <- data.frame(x=unlist(df.oral.A))
df.oral$DPHP <- c("Oral") 
names(df.oral)=c("Value","DPHP")


df.all.A = as.data.frame(subset(df.all, select = -c(DPHP) ))
df.all <- data.frame(x=unlist(df.all.A))
df.all$DPHP <- c("all") 
names(df.all)=c("Value","DPHP")



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
  ggplot (df, aes(x = as.numeric(log.value), y = as.factor(ggjoy),fill = DPHP)) + # fill, fill one color for each species
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

ggsave("Exposure_updateda1.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result",
       width = 25, height = 20, units = "cm",dpi=320)



###### get mean & percentile ##########
## Inhalation ##
IE = df.inhal$Value
mean_inhal = mean(df.inhal$Value)
quantile(df.inhal$Value, c(.05,0.5, .95))

## Oral ###
OE = df.oral$Value
mean_oral = mean(df.oral$Value)
quantile(df.oral$Value, c(.05,0.5, .95))

## Dermal ##
SE = df.dermal$Value
mean_dermal = mean(df.dermal$Value)
quantile(df.dermal$Value, c(.05,0.5, .95))

## Overall ##
df.all <- df.all$Value
quantile(df.all, c(0.01,.05,0.25, 0.5, 0.75,.95,0.99,1))

saveRDS(IE,file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/IE.rds')
saveRDS(OE,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/OE.rds')
saveRDS(SE,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/SE.rds')
saveRDS(BW,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/BW.rds')
saveRDS(df,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/df.rds')
saveRDS(df.all,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/df.all.rds')

df.all = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/df.all.rds')
OE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/OE.rds')

quantile(df.all, c(0.01,.05,0.25, 0.5, 0.75,.95,0.99,1))
quantile(OE, c(0.01,.05,0.25, 0.5, 0.75,.95,0.99,1))
