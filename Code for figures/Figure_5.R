
## Load libraries
library(mrgsolve) # Needed to run the main PBPK code
library(magrittr) # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2) # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(dplyr) # Needed for the pipe %>% operator
library(reshape) # melt function to reshape the table
library(scales) # for plotting the figure
library(gridExtra) # for plotting the figure
library(grid) # for plotting the figure
library(lattice) # for plotting the figure
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
#install.packages("FME")
library(FME)         # Package for MCMC simulation and model fitting
#install.packages("minpack.lm")
library(minpack.lm)  # Package for model fitting
#install.packages("reshape")
library(reshape)     # Package for melt function to reshape the table
#install.packages("truncnorm")
library(truncnorm)   # Package for the truncated normal distribution function 
#install.packages("EnvStats")
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
#install.packages("invgamma")
library(invgamma)    # Package for inverse gamma distribution function
#install.packages("foreach")
library(foreach)     # Package for parallel computing
#install.packages("doParalle")
library(doParallel)  # Package for parallel computing
#install.packages("bayesplot")
library(bayesplot)   # Package for MCMC traceplot
library(dplyr)

rm(list = ls())


########### load exposure data
df_DEP <- readRDS(file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/df.rds')
df_DEHA <- readRDS(file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/df.rds')
df_DEHP <- readRDS(file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/df.rds')
df_DINCH <- readRDS(file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/df.rds')
df_DPHP <- readRDS(file = 'C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/df.rds')



df_DEP$Chemical <- c("DEP")
df_DEHA$Chemical <- c("DEHA")
df_DEHP$Chemical <- c("DEHP")
df_DINCH$Chemical <- c("DINCH")
df_DPHP$Chemical <- c("DPHP")

colnames(df_DEP)[2] <- "Route"
colnames(df_DEHA)[2] <- "Route"
colnames(df_DEHP)[2] <- "Route"
colnames(df_DINCH)[2] <- "Route"
colnames(df_DPHP)[2] <- "Route"

head(df)

df          <- rbind.data.frame (df_DEP, df_DEHA, df_DEHP, df_DINCH, df_DPHP)
df$Chemical <- as.character(df$Chemical)

df$Chemical <- factor(df$Chemical, levels=c("DPHP", "DINCH", "DEHP","DEHA","DEP"))

######################## plot ###############################

library(ggplot2)
library(ggridges)



windowsFonts(Times=windowsFont("Times New Roman")) # Abbreviate the font Times New Roman as Times

#df$ggjoy = c("A") # Create a new column ggjoy in df table. Without this, the plot p1 will have four layers for four species overlapped in each panel, which does not look good.
# With this, y axis is A for all species in one layer in one panel.

p1 = 
  ggplot (df, aes(x = as.numeric(log.value), y = as.factor(Chemical), color = Route, fill = Route)) + 
  geom_density_ridges (
    scale = 0.95, size = 1, rel_min_height = 0.01, alpha = 0.6
    ) + # over size of the entire ggridge plot
  scale_y_discrete(expand = c(0, 0), name = "") +
  scale_x_continuous(
    breaks = c(-9:4), limits = c(-8, 4),
    expand = c(0.01, 0), name = expression(paste("Log [Daily Intake] (", mu,"g/kg bw/d)"))) +
  scale_fill_manual(values = c("grey", "#2E9FDF","#7CAE00")) +  
  scale_color_manual(values = c("#0A0A0A","#0A0A0A","#0A0A0A")) +
  guides(fill = guide_legend(
                override.aes = list(
                        fill = c("grey", "#2E9FDF","#7CAE00"),
                        color = NA, point_color = NA)
                       )
  )
# labs (x = expression(paste("Log [Concentration] (", mu,"g/L)")))



p1=
  p1 + theme_ridges()+ 
  theme (
    plot.background         = element_blank(),
    text                    = element_text (family = "Times",size = 20),
    panel.background        = element_rect (fill   = "white"),
    panel.grid.major.x      = element_blank(),
    axis.text               = element_text (size   = 20, colour = "black"),  
    axis.title.x            = element_text (hjust = 0.5, face="bold", size = 20),
    axis.ticks              = element_line(color="#5A5A5A", size = 1) , 
    #axis.line               = element_line(color = "#5A5A5A", size = 1),
    axis.ticks.length       = unit(.2, "cm"),
    legend.title            = element_blank(),
    #legend.justification    =  c("right", "top"),
    legend.position         = "top",
    legend.text             = element_text (size = 20,face="bold")) 

p1


ggsave("Exposure_summaryABCDEFGHIJ.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot",
       width = 20, height = 20, units = "cm",dpi=800)

