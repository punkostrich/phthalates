dev.off()
library(ggplot2)
library(dplyr)

rm(list=ls())

df.DPHP2   <- readRDS (file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/result/df24.rds')
df.DPHP   <- na.omit(df.DPHP2)
df.DPHP.model   <-  as.data.frame(df.DPHP [df.DPHP $type == '#1 Model Prediction',])
df.DPHP.survey  <-  as.data.frame(df.DPHP [df.DPHP $type == 'ESB 2012',])
df.DPHP.model   <-  as.data.frame(sample_n(df.DPHP.model, 10000))
df.DPHP         <-  rbind.data.frame (df.DPHP.model, df.DPHP.survey)

quantile(df.DPHP.model$Value, c(.05,0.5, .95))
quantile(df.DPHP.survey$Value, c(.05,0.5, .95))


df.DEHP2   <- readRDS (file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/result/df24.rds')
df.DEHP   <- na.omit(df.DEHP2)
df.DEHP.model   <-  as.data.frame(df.DEHP [df.DEHP $type == '#1 Model Prediction',])
df.DEHP.survey  <-  as.data.frame(df.DEHP [df.DEHP $type == 'NHANES 2015-2016',])
df.DEHP.model   <-  sample_n(df.DEHP.model, 10000)
df.DEHP         <-  rbind.data.frame (df.DEHP.model, df.DEHP.survey)

quantile(df.DEHP.model$Value, c(.05,0.5, .95))
quantile(df.DEHP.survey$Value, c(.05,0.5, .95))

df.DEHA2  <- readRDS (file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/result/df24.rds')
df.DEHA   <- na.omit(df.DEHA2)
df.DEHA.model   <-  as.data.frame(df.DEHA[df.DEHA $type == '#1 Model Prediction',])
df.DEHA.model   <-  sample_n(df.DEHA.model, 10000)
df.DEHA         <-  rbind.data.frame (df.DEHA.model)


DEHA_Survey <- c(as.numeric(0),"NHANES 2015-2016","Survey 2015-2016",as.numeric(0), "DEHA")
names(DEHA_Survey ) <-  c("Value","type","year", "log.value","compound")
df.DEHA   <- rbind(df.DEHA, DEHA_Survey )  

quantile(df.DEHA.model$Value, c(.05,0.5, .95))


df.DINCH2  <- readRDS (file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/result/df24.rds')
df.DINCH   <- na.omit(df.DINCH2)
df.DINCH.model   <-  as.data.frame(df.DINCH[df.DINCH $type == '#1 Model Prediction',])
df.DINCH.survey  <-  as.data.frame(df.DINCH[df.DINCH $type == 'NHANES 2015-2016',])
df.DINCH.model   <-  sample_n(df.DINCH.model, 10000)
df.DINCH         <-  rbind.data.frame (df.DINCH.model, df.DINCH.survey)

quantile(df.DINCH.model$Value, c(.05,0.5, .95))
quantile(df.DINCH.survey$Value, c(.05,0.5, .95))


df.DEP2    <- readRDS (file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/result/df24.rds')
df.DEP   <- na.omit(df.DEP2)
df.DEP.model   <-  as.data.frame(df.DEP[df.DEP $type == '#1 Model Prediction',])
df.DEP.survey  <-  as.data.frame(df.DEP[df.DEP $type == 'NHANES 2015-2016',])
df.DEP.model   <-  sample_n(df.DEP.model, 10000)
df.DEP         <-  rbind.data.frame (df.DEP.model, df.DEP.survey)

quantile(df.DEP.model$Value, c(.05,0.5, .95))
quantile(df.DEP.survey$Value, c(.05,0.5, .95))


df.DEP$compound   <- c("MEP")
df.DEHP$compound  <- c("MEHP")
df.DINCH$compound <- c("MHNCH")
df.DPHP$compound  <- c("oxo-MPHP")
df.DEHA$compound  <- c("5OH-MEHA")


df <- rbind.data.frame (df.DEP, df.DEHA, df.DEHP, df.DINCH, df.DPHP)
df$compound <- as.character(df$compound)
df$compound <- factor(df$compound, levels=c("MEP","5OH-MEHA", "MEHP", "MHNCH","oxo-MPHP"))

df$type[df$type == '#1 Model Prediction'] <- c("Model Prediction")
df$type[df$type == 'ESB 2012'] <- c("Biomarker Survey")
df$type[df$type == 'NHANES 2015-2016'] <- c("Biomarker Survey")
df$type <- as.character(df$type)
df$type <- factor(df$type, levels=c("Model Prediction", "Biomarker Survey"))




############# PLOT #################

# define the summary function
f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

library(ggplot2)
library(CRSSIO)

windowsFonts(Times=windowsFont("Times New Roman")) 

#install.packages("remotes")
#remotes::install_github("BoulderCodeHub/CRSSIO")

pd = position_dodge(width = 0.4)

p1 = 
  ggplot (df, aes(x = as.factor(compound), y = as.numeric(log.value), color = type, fill = type)) + 
  stat_boxplot_custom(geom="errorbar", position=pd, width=0.3,size = 1) +
  stat_boxplot_custom(geom="boxplot", position=pd, width=0.4,size = 1,
                      qs = c(0.05, 0.25, 0.5, 0.75, 0.95),outlier.shape = NA) +
 # geom_boxplot(width=0.4,position=pd,size = 0.7,outlier.shape = NA) +
  scale_fill_manual(values = c("#00AFBB", "#FC4E07"))+
  scale_color_manual(values = c("#0A0A0A","#0A0A0A","#0A0A0A")) + 
  scale_y_continuous(minor_breaks = seq(-3 , 7, 1), breaks = seq(-3, 7,2))+
  coord_cartesian(ylim = c(-2, 3.2))
  #scale_y_discrete(expand = c(0, 0), expression(paste("Log [Urinary Conc] (", mu,"g/L)"))) +

#p1 <- p1+ scale_y_continuous(limits = c(-3,5.2),expand = c(0,0))


p1 <- p1 + 
  
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major.y      = element_line(colour = "grey", size = 0.5),
    panel.grid.minor.y      = element_line(colour = "grey", linetype = 4, size = 0.3),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black"),    # tick labels along axes 
    axis.title.y            = element_text (size   = 20, colour = "black"),   # label of axes
    legend.title            = element_blank(),
    legend.justification    =  c("right", "top"),
    axis.title.x            = element_blank(),
    legend.position         = c(0.5, 0.2),
    legend.text             = element_text (size = 20),
    legend.text.align       = 0) + 
    labs (y = expression(paste("Log [Urine Conc] (", mu,"g/L)")))

#labs (x = expression(paste("Observed value (", mu,"g/L)")), y = expression(paste("Predicted value (", mu,"g/L)")))

p1

 #= MC.DPHP.plasma  + Theme.Fig + labs (x ="Time (hours)", y=expression(paste("Concentration in plasma (", mu,"g/L)")))
ggsave("BoxplotABCDEFGHIjk.tiff",scale = 1,
       plot = p1,
       path = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot",
       width = 15, height = 12, units = "cm",dpi=600)

