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
M <- 500    #M sets 
N <- 1000   #N exposure


set.seed(123)
############################################# outer loop ##############################################
## production volume #(Sc, shared among members of the entire cohort)
SD_volume <- 0.3* 0.453592 * 20 * 1E6   # kg # 30% CV
volume.sample <-rtruncnorm (n=M, a = qnorm(0.05, mean = 0.453592 *20 * 1E6 , sd = SD_volume  ), 
                            b = qnorm(0.95, mean = 0.453592 * 20 * 1E6 , sd = SD_volume ), 
                            mean = 0.453592 * 20 * 1E6 , sd = SD_volume )             # kg 
quantile(volume.sample, c(0,0.05,.5, 0.75, 0.9, 0.95,1)) 


## particle-air partition coefficient normal  #(Sc, shared among members of the entire cohort)
SD_kp <- (9.1E-5 - 4.6E-5) / qnorm(0.75)
Kp.sample <-rtruncnorm (n=M, a = 1E-7, 
                 b = qnorm(0.95, mean = 4.6E-5, sd = SD_kp), 
                 mean = 4.6E-5, sd = SD_kp)              ## m3/ug
quantile(Kp.sample, c(0.05,.5, 0.75, 0.9, 0.95)) 


## dermal absorption percentage 
abs_mean.sample  <-  runif(n=M, min = 0.007, max = 0.01)   ## %
abs_CV.sample    <-  runif(n=M, min = 0.3, max = 1)     #  (ScUc)

######## household size ######## # unshared (ScUs)
##### emission surface are per household
A_mean.sample  <-  runif(n=M, min = 0.001, max = 5)   ## m2
A_CV.sample    <-  runif(n=M, min = 0.3, max = 1)     #  (ScUc)

occurance.sample <- runif(n=M, min = 0.3, max = 0.8) 

dens_DEP <- 1.12 * 1E3 #kg/m3

Sc.sample <- cbind(volume.sample,abs_mean.sample, abs_CV.sample , occurance.sample ,
                   Kp.sample, A_mean.sample, A_CV.sample)



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
  
  volume   = fixed['volume.sample']
 # C_pcp    = fixed['C_pcp.sample']
  abs_mean = fixed['abs_mean.sample'] 
  abs_CV   = fixed['abs_CV.sample'] 
  A_mean = fixed['A_mean.sample'] 
  A_CV   = fixed['A_CV.sample'] 
  Kp     = fixed['Kp.sample']
  occurance = fixed['occurance.sample']
 # y0     = fixed['y0.sample'] 
 
  
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
  abs <- rtruncnorm (n=N, a = 0.0007, 
                      b = 1, 
                      mean = abs_mean, sd =  abs_CV * abs_mean)  # ScUc
  
  A <- rtruncnorm (n=N, a = 0.0001, 
                     b = qnorm(0.95, mean = A_mean, sd =  A_CV * A_mean), 
                     mean = A_mean, sd =  A_CV * abs_mean)  # ScUc
  
  
  ##### percentage in personal care products log-normal (Uc, shared among members of the entire cohort)
  #dCpcp <- c(0.002, 0.03)    #1% 28.6%
  #PercdCpcp <- c(0.5, 0.95)   #50th, 97.5th percentile
  #c <- get.lnorm.par(p= PercdCpcp,q=dCpcp)
  Cpcp_ori <- rlnormTrunc(n=M, meanlog = -6.214137, sdlog = 1.646638,max = 0.2)## g/g
  #quantile(Cpcp_ori*100, c(.05, .1,.25,.5, 0.75, 0.9, 0.95,0.975,0.999,1)) 
  C_pcp <- Cpcp_ori * 1E6                            # % to ug/g (ppm)
  #quantile(C_pcp.sample, c(0.05,.5, 0.75, 0.9, 0.95)) 

  ## vapor from product #(Sc, shared among members of the entire cohort)
  y0  = min(19113, 19113 * Cpcp_ori / dens_DEP / (Cpcp_ori/dens_DEP + (1-Cpcp_ori)/1E3) * 3.7) ## ug/m3, 19113 is the vapor pressure of DEP
  #quantile(y0.sample, c(0.05,.5, 0.75, 0.9, 0.95,1)) 
  
  
  ## convective mass transfer coefficient log-normal # unshared (Us)
  #dh <- c(3.6, 7.2)
  #Percdh <- c(0.5, 0.95)
  #b <- get.lnorm.par(p= Percdh,q=dh)
  h <- rlnorm(n=N, meanlog = 1.2809338 , sdlog = 0.4214028)
  #quantile(h, c(.5, 0.75, 0.9, 0.95))  # check ## m/h
  
  
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
  
  ## Calculate DEP conc in air 
  Cair <- h * y0 * A/(h * A + (1 + Kp * TSP) * Q)    ## ug/m3
  #quantile(Cair , c(0.05,.5, 0.75, 0.9, 0.95))
  ## inhalation exposure
  IE <- (Cair * InhR * ED + Cair * Kp * TSP * InhR * ED)/BW   ## ug/day/kg
  
  
  
  
  ######################################## 2 ####################################
  ################################ Exposure scenario #Oral ######################
  
  
  #dust oral # unshared (Us)
  IngR  <- runif(n=N, min = 30, max = 70)  ## mg/d
  
  Kdust <- Kp/8.32            ## m3/ug #(Sc, shared among members of the entire cohort)
  #quantile(Kdust, c(0.05,.5, 0.75, 0.9, 0.95, 0.975))  # check
  
  OE <- (Cair * Kdust * IngR)/BW   ## ug/day/kg
  #quantile((Cair * Kdust * IngR)/BW, c(.5, 0.75, 0.9, 0.95)) 
  #quantile(OE, c(.5, 0.75, 0.9, 0.95))

  
  ############################## 3 ###################################
  ##################### Exposure scenario #Dermal ####################
  ## dermal absorption coefficient # unshared (Us)
  SD_kpg <- 1.56 * 5
  kpg <- rtruncnorm (n=N, a = 0.1e-03, 
                     b = qnorm(0.95, mean =1.56, sd = SD_kpg ), 
                     mean = 1.56 , sd = SD_kpg )              ## m/d
  
  ## fraction of skin exposed  # unshared (Us)
  SD_Fr <- 0.24 * 1
  Fr  <- rtruncnorm (n=N, a = 0.11, 
                     b = 1, 
                     mean = 0.24 , sd = SD_Fr ) ## - # unshared (Us)
  #quantile(Fr, c(0.05,.5, 0.75, 0.9, 0.95))
  
  ### g/day cosmatics log-normal (unshared)
  
  #dcos <- c(5.6, 18.7)
  #Percdcos <- c(0.5, 0.90)
  #c <- get.lnorm.par(p= Percdcos,q=dcos)
  cosmatic <- rlnormTrunc(n=N, meanlog = 1.7227691, sdlog = 0.9408541, max = 25) ## g/d
  #quantile(cosmatic, c(0.05,.5, 0.75, 0.9, 0.95,1)) 
  
  SE <- (Cair * kpg * BSA * Fr * ED)/BW + (C_pcp * cosmatic * occurance * abs)/BW ## ug/day/kg 5% absorption
  quantile(SE, c(0.05,.5, 0.75, 0.9, 0.95))
  
  
  
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
df.inhal = subset(df.inhal, select = -c(V1) )
df.inhal$DEP <- c("Inhalation") 
#names(df.inhal)=c("Value","DEP")

## Oral
df.oral<-as.data.frame(df.OE)
df.oral = subset(df.oral, select = -c(V1) )
df.oral$DEP <- c("Oral") 
#names(df.oral)=c("Value","DEP")

## Dermal
df.dermal<-as.data.frame(df.SE)
df.dermal = subset(df.dermal, select = -c(V1) )
df.dermal$DEP <- c("Dermal") 
#names(df.dermal)=c("Value","DEP")

## all
df.all <-as.data.frame(df)
df.all = as.data.frame(subset(df.all, select = -c(V1) ))

quants <- c(0.50,0.95)
df.stat = apply( df.all, 2 , quantile , probs = quants ) # get median and 95% quantile
df.stat <- as.data.frame(df.stat)

# get measured values from literature
x.median <- c(2)
x.95 <- c(20)
survey.median <- quantile(x.median, c(.95))
survey.95 <- quantile(x.95, c(.95))

# get residues
df.2.1 <- abs(df.stat[c("50%"),] - survey.median)
df.2.2 <- abs(df.stat[c("95%"),] - survey.95)
df.2.3 <- as.data.frame(rbind(df.2.1, df.2.2))

# order by residues, for quality check
y = df.2.3[,order(df.2.3[nrow(df.2.3),])]

df.2.3$DEP <- c("quantile")
df.all$DEP <- c("all") 

library(data.table) #data.table_1.9.5
df.2 <- rbindlist(list(df.inhal, df.oral, df.dermal, df.all, df.2.3),use.names=FALSE )
df.2R <- as.data.frame(df.2)
df.2A = as.data.frame(subset(df.2R, select = -c(DEP) ))

df.95 <- df.2A[,order(df.2A[(nrow(df.2A)-1),])]       # order column based on 50% quantile
# select first 100 columns
df.95_50c <- df.95[, 1:100]

df.50 <- df.95_50c[,order(df.95_50c[(nrow(df.95_50c)),])] # order column based on 95% quantile
# select first 50 columns
df.f <- df.50[, 1:50]
df.f2 <- cbind.data.frame(df.f, df.2R$DEP)
names(df.f2)[length(names(df.f2))]<-"DEP" 

# extrace subset data frame for inhal, dermal,  oral, all sepearately
df.inhal <- df.f2[df.f2$DEP == 'Inhalation', ]
df.dermal <- df.f2[df.f2$DEP == 'Dermal', ]
df.oral <- df.f2[df.f2$DEP == 'Oral', ]
df.all <- df.f2[df.f2$DEP == 'all', ]

# mix all data within each category
df.inhal.A = as.data.frame(subset(df.inhal, select = -c(DEP) ))
df.inhal <- data.frame(x=unlist(df.inhal.A))
df.inhal$DEP <- c("Inhalation") 
names(df.inhal)=c("Value","DEP")

df.dermal.A = as.data.frame(subset(df.dermal, select = -c(DEP) ))
df.dermal <- data.frame(x=unlist(df.dermal.A))
df.dermal$DEP <- c("Dermal") 
names(df.dermal)=c("Value","DEP")


df.oral.A = as.data.frame(subset(df.oral, select = -c(DEP) ))
df.oral <- data.frame(x=unlist(df.oral.A))
df.oral$DEP <- c("Oral") 
names(df.oral)=c("Value","DEP")


df.all.A = as.data.frame(subset(df.all, select = -c(DEP) ))
df.all <- data.frame(x=unlist(df.all.A))
df.all$DEP <- c("all") 
names(df.all)=c("Value","DEP")



#df.all$DEP <- c("All") 
#names(df.all)=c("Value","DEP")




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
  ggplot (df, aes(x = as.numeric(log.value), y = as.factor(ggjoy),fill = DEP)) + # fill, fill one color for each species
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

ggsave("Exposure_updated3.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP",
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
quantile(df.all$Value, c(0.01,.05,0.25, 0.5, 0.75,.95,0.99,1))

saveRDS(IE,file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/IE.rds')
saveRDS(OE,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/OE.rds')
saveRDS(SE,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/SE.rds')
saveRDS(BW,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/BW.rds')
saveRDS(df,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/df.rds')
saveRDS(df.all,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/df.all.rds')


IE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/IE.rds')
OE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/OE.rds')
SE = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/SE.rds')
df = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/df.rds')
df.all = readRDS(file ='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/df.all.rds')

quantile(df.all, c(0.01,.05,0.25, 0.5, 0.75,.95,0.99,1))

