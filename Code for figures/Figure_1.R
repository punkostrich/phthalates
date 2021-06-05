
library(circlize)

rm(list=ls())
circos.clear()

#png(filename = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot/my_plot.png") 
tiff(filename = "C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/Code for plot/ccplotABCDEFG.jpg", units="in", width=5, height=5, res=500)

group = c("IN SILICO", "IN VIVO", "IN VITRO", "Human")

set.seed(999)
n = 2000
df = data.frame(sectors = sample(group[1:4], n, replace = TRUE),
                x = rnorm(n), y = runif(n))

circos.par("track.height" = 0.1)
circos.initialize(df$sectors, x = df$x)

circos.track(df$sectors, y = df$y)#,
 #            panel.fun = function(x, y) {
# #              circos.text(CELL_META$xcenter, 
#                           CELL_META$cell.ylim[2] + mm_y(5), 
#                           CELL_META$sector.index)
#            })

##########     blue, yellow,  red , green,
col = rep(c("#E7B800", "cornflowerblue","#FC4E07","#00AFBB")) ##2E9FDF"), 1)
#col = rep(c("#FC4E07", "#00AFBB","#FC4E07","#00AFBB")) ##2E9FDF"), 1)
circos.trackPoints(df$sectors, df$x, df$y, col = col, pch = 16, cex = 0.5)
#circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)

bgcol = rep(c("#EFEFEF", "#CCCCCC"), 2)
circos.track(df$sectors, y = df$y)     ############ 2ND LAYER
circos.update(sector.index = "IN VIVO", track.index = 2,bg.col = "#EFEFEF")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "IN SILICO", col = "black", facing = "bending.inside")
circos.update(sector.index = "IN SILICO", track.index = 2,bg.col = "#EFEFEF")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "IN VITRO", col = "black", facing = "bending.outside")
circos.update(sector.index = "IN VITRO", track.index = 2,bg.col = "#EFEFEF")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "IN VIVO", col = "black", facing = "bending.inside")
circos.update(sector.index = "Human", track.index = 2,bg.col = "#EFEFEF")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "IN SILICO", col = "black", facing = "bending.outside")


circos.track(df$sectors, y = df$y)   ###### 3RD LAYER
circos.update(sector.index = "IN VIVO", track.index = 3)#,bg.col = "#CCCCCC")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "Exposure Modeling", col = "black", facing = "bending.inside")
circos.update(sector.index = "Human", track.index = 3)#,bg.col = "#EFEFEF")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "PK Modeling", col = "black", facing = "bending.outside")
circos.update(sector.index = "IN SILICO", track.index = 3)#, bg.col = "#CCCCCC")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "ToxCast Assays", col = "black", facing = "bending.outside")
circos.update(sector.index = "IN VITRO", track.index = 3)#, bg.col = "#EFEFEF")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "Population Module", col = "black", facing = "bending.inside")


circos.track(df$sectors, y = df$y)   ###### 4th LAYER
circos.update(sector.index = "IN VIVO", track.index = 4,bg.col = "#00AFBB", bg.border = "black")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "Product content data", col = "white", facing = "bending.inside")

circos.update(sector.index = "Human", track.index = 4)
circos.text(CELL_META$xcenter, CELL_META$ycenter, "Pharmacokinetic data", col = "black", facing = "bending.outside")



circos.update(sector.index = "IN SILICO", track.index = 4, bg.col = "cornflowerblue", bg.border = "black")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "Endpoint data ", col = "white", facing = "bending.outside")

circos.update(sector.index = "IN VITRO", track.index = 4)
circos.text(CELL_META$xcenter, CELL_META$ycenter, "Biomonitoring data", col = "black", facing = "bending.inside")



#bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
#circos.trackHist(df$sectors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)


#circos.track(df$sectors, x = df$x, y = df$y,
 #            panel.fun = function(x, y) {
 #              ind = sample(length(x), 10)
 #              x2 = x[ind]
 #              y2 = y[ind]
 #              od = order(x2)
 #              circos.lines(x2[od], y2[od])
  #           })



#circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#  xlim = CELL_META$xlim
#  ylim = CELL_META$ylim
#  breaks = seq(xlim[1], xlim[2], by = 0.1)
#  n_breaks = length(breaks)
#  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
#              breaks[-1], rep(ylim[2], n_breaks - 1),
#              col = rand_color(n_breaks), border = NA)
#})



#circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
ylim = c(0, 1)
 xlim = CELL_META$xlim
 ylim = CELL_META$ylim
 breaks = seq(xlim[1], xlim[2], by = 0.1)
 n_breaks = length(breaks)
            col = rand_color(n_breaks)
#})

library(shape)            
      #      col = rep(c("#E7B800", "cornflowerblue","#FC4E07","#00AFBB")           
#circos.link("IN VIVO", c(0.5, 2.2), "Human", c(-1,-2.6), col = rand_color(8, transparency = 0.7), #"#00AFBB", 
       #     border = NA, h = 0.8)

circos.link("IN VIVO", c(-0.8, 1), "Human", c(-1,-2.6), col = "#00AFBB80",
            border = NA, h = 0.8)

circos.link("IN SILICO", c(-1, 1.4), "Human", c(1.2,2.9), col = "#E7B80080",
            border = NA, h = 0.8)

circos.link("IN VITRO", c(0,-2.4), "Human", c(-1,1.2), col = "#E7B80080",#col = "#E7B800",
            border = NA, h = 0.8)

#circos.link("IN VIVO", -1.2, "IN SILICO", 2, col = "#00AFBB80", h = 0.8, lwd = 4, lty = 2,
           # directional = 1, arr.length = ifelse( 0.02, 0.4),arr.width = 0.8, arr.type = "big.arrow")


circos.link("IN VIVO", c(-1.9,-0.8), "IN SILICO", c(1.6,2.8), col = "#00AFBB40", h = 0.8, lwd = 2, lty = 2,border = "cyan4",
            directional = 1, arr.length = ifelse(0.02, 0.1),arr.width = 0.8, arr.type = "big.arrow")

circos.link("IN VIVO", c(-1.35), "IN SILICO", c(2.2), col = "#00AFBB", h = 0.8, lwd = 1, lty = 0,border = NA,
            directional = 1, arr.length = ifelse(0.02, 0.5),arr.width = 1)

circos.link("IN VIVO", c(-1.9,-0.8), "IN VITRO", c(0.4,1.6), col = "#00AFBB40",h = 0.8, lwd = 2, lty = 2,border = "cyan4",
            directional = 1, arr.length = ifelse( 0.2, 0.1),arr.width = 0.8, arr.type = "big.arrow")

circos.link("IN VIVO", c(-1.35), "IN VITRO", c(1), col = "#00AFBB",h = 0.8, lwd = 2, lty = 0,border = NA,
            directional = 1, arr.length = ifelse( 0.2, 0.5),arr.width = 1)

#circos.link("a", 0, "b", 0, h = 0.4)
#circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red",
#            border = "blue", h = 0.2)
#circos.link("e", 0, "g", c(-1,1), col = "green", border = "black", lwd = 2, lty = 2)


dev.off()

