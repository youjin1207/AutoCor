library(geosphere)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(spdep)
library(gridExtra)
###
source("Code/MoranI.R")
source("Code/Phi.R")
###
load("Data/analysis_dat.RData")
#names(analysis_dat)
## locations
latitude = analysis_dat$Fac.Latitude
longitude = analysis_dat$Fac.Longitude
dist.matrix = distm(cbind(longitude, latitude), cbind(longitude, latitude), fun = distHaversine)
weights = max(dist.matrix) / dist.matrix
diag(weights) = 0
summary(as.numeric(weights))
weights = ifelse(weights > 10, 10, weights) # weight not exceeding 10

white.moran = make.permute.moran(weights, analysis_dat$PctWhite, 500)
hispanic.moran = make.permute.moran(weights,  analysis_dat$PctHisp, 500)
black.moran = make.permute.moran(weights,  analysis_dat$PctBlack, 500)

### Phi
## category 1 
summary(analysis_dat$PctBlack)
summary(analysis_dat$PctHisp)
summary(analysis_dat$PctWhite)
max.percentage = pmax(analysis_dat$PctBlack, analysis_dat$PctHisp, analysis_dat$PctWhite)
largest.race = ifelse(max.percentage == analysis_dat$PctBlack, 2, 1)
largest.race = ifelse(max.percentage == analysis_dat$PctHisp, 3, largest.race)
race.Phi = make.permute.Phi(weights, largest.race, 500) 
## category 2 
ethnic.indi = ifelse(analysis_dat$PctBlack > 0.10 & analysis_dat$PctHisp > 0.10, 1, 4)
ethnic.indi = ifelse(analysis_dat$PctBlack > 0.10 & analysis_dat$PctHisp <= 0.10, 2, ethnic.indi)
ethnic.indi = ifelse(analysis_dat$PctBlack <= 0.10 & analysis_dat$PctHisp > 0.10, 3, ethnic.indi)
ethnic.indi = ifelse(analysis_dat$PctBlack <= 0.10 & analysis_dat$PctHisp <= 0.10, 4, ethnic.indi)
ethnic.indi.Phi = make.permute.Phi(weights, ethnic.indi, 500)


### join count
nb_15nn = knn2nb(knearneigh( cbind(analysis_dat$Fac.Longitude, analysis_dat$Fac.Latitude), k = 15 ))
nblist = nb2listw(nb_15nn)
## category 1
mc.test = joincount.mc(largest.race, nblist, nsim = 500)
BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, mc.test[[3]]$statistic)
BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, mc.test[[3]]$p.value)
BB.stat = formatC(BB.stat, 2, format = "f")
BB.pval = formatC(BB.pval, 4, format = "f")
mat = rbind(as.integer(table(largest.race)), BB.stat, BB.pval)
rownames(mat) = c("n","Join-count statistic", "P-value (permutation)")
colnames(mat) = c("White", "Hispanic", "African-American")
print(xtable(mat))
## category 2
mc.test = joincount.mc(ethnic.indi, nblist, nsim = 500)
BB.stat = c(mc.test[[1]]$statistic, mc.test[[2]]$statistic, 
            mc.test[[3]]$statistic, mc.test[[4]]$statistic)
BB.pval = c(mc.test[[1]]$p.value, mc.test[[2]]$p.value, 
            mc.test[[3]]$p.value, mc.test[[4]]$p.value)
BB.stat = formatC(BB.stat, 2, format = "f")
BB.pval = formatC(BB.pval, 4, format = "f")
mat = rbind(as.integer(table(ethnic.indi)), BB.stat, BB.pval)
rownames(mat) = c("n", "Join-count statistic", "P-value (permutation)")
colnames(mat) = c(expression(paste(AA>10*"%", " & ", HP > 10*"%")),
                  expression(paste(AA>10*"%", " & ", HP <= 10*"%")),
                  expression(paste(AA<=10*"%", " & ", HP > 10*"%")),
                  expression(paste(AA<=10*"%", " & ", HP <= 10*"%")))
print(xtable(mat))

#### Figure ####
usa = map_data("usa")
states = map_data("state")
gg1 = ggplot(data = states) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "snow2", color = "white") + 
  coord_fixed(1.3) +
  guides(fill=FALSE)

labs = data.frame(
  long = analysis_dat$Fac.Longitude,
  lat = analysis_dat$Fac.Latitude,
  white = analysis_dat$PctWhite,
  hispanic = analysis_dat$PctHisp,
  black = analysis_dat$PctBlack, 
  small = analysis_dat$small_nunits, 
  medium = analysis_dat$med_nunits
)

summary(labs$white)  
ditch_the_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank()
)

breaks = seq(0, 1, 0.2)
# Create the textGrobs (Moran's I or Phi)
Moran1 = textGrob(expression(paste("Moran's I (", Phi, ") : 30.99")))
Moran2 = textGrob(expression(paste("Moran's I (", Phi, ") : 93.36")))
Moran3 = textGrob(expression(paste("Moran's I (", Phi, ") : 20.63")))

Pval1 = textGrob(paste("P-value (permutation) :", formatC(white.moran[3], 4, format = "f")))
Pval2 = textGrob(paste("P-value (permutation) :", formatC(white.moran[3], 4, format = "f")))  
Pval3 = textGrob(paste("P-value (permutation) :", formatC(white.moran[3], 4, format = "f")))


white = gg1 + theme_bw() +  ditch_the_axes + 
  ggtitle("Proportion of White") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  geom_point(data = labs, aes(x = long, y = lat, colour = white),size = 2) +  
  scale_colour_gradient(low = "powderblue", high = "black", breaks = breaks) +
  theme(legend.position="none") +  
  annotate("text", x = -100, y =22, size = 7,
           label = as.character(expression(paste("Moran's I:",  " 30.99"))), parse = TRUE) + 
  annotate("text", x = -100, y = 19, size = 7,
           label = paste("P-value (permutation) :",  formatC(white.moran[3], 4, format = "f")), parse = TRUE)

hispanic = gg1 + theme_bw() +  ditch_the_axes + 
  ggtitle("Proportion of Hispanic") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  geom_point(data = labs, aes(x = long, y = lat, colour = hispanic),size = 2) +
  scale_colour_gradient(low = "powderblue", high = "black", breaks = breaks) +
  theme(legend.position="none") + 
  annotate("text", x = -100, y =22, size = 7,
           label = as.character(expression(paste("Moran's I:",  " 93.36"))), parse = TRUE) + 
  annotate("text", x = -100, y = 19, size = 7,
           label = paste("P-value (permutation) :",  formatC(hispanic.moran[3], 4, format = "f")), parse = TRUE)

black = gg1 + theme_bw() +  ditch_the_axes + 
  ggtitle("Proportion of African American") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  geom_point(data = labs, aes(x = long, y = lat, colour = black),size = 2) +
  scale_colour_gradient(low = "powderblue", high = "black", breaks = breaks) +
  theme(legend.position="none") + 
  annotate("text", x = -100, y =22, size = 7,
           label = as.character(expression(paste("Moran's I:",  " 20.63"))), parse = TRUE) + 
  annotate("text", x = -100, y = 19, size = 7,
           label = paste("P-value (permutation) :",  formatC(hispanic.moran[3], 4, format = "f")), parse = TRUE)


grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position,
                                     legend.title=element_blank(), 
                                     legend.text=element_text(size=20),
                                     legend.key = element_rect(size = 10),
                                     legend.key.size = unit(2, 'lines')))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

pdf("Figure/threerace.pdf", height = 6, width = 25)
grid_arrange_shared_legend(white, hispanic, black, ncol = 3, nrow = 1,
                           position = "right") 
dev.off()  


## category 1
max.percentage = pmax(analysis_dat$PctBlack, analysis_dat$PctHisp, analysis_dat$PctWhite)
largest.race = ifelse(max.percentage == analysis_dat$PctBlack, 2, 1)
largest.race = ifelse(max.percentage == analysis_dat$PctHisp, 3, largest.race)
largest.race <- factor(largest.race, levels=c(1, 3, 2), labels=c("White", "Hispanic", "African-American"))
cate1 = gg1 + theme_bw() +  ditch_the_axes + 
  ggtitle("Dominant ethnic group") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  geom_point(data = labs, aes(x = long, y = lat, colour = as.factor(largest.race)),size = 2) +
  theme(legend.title=element_blank(),  
        legend.text=element_text(size=20),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(2, 'lines')) +
  guides(color = guide_legend(override.aes = list(size=5))) + 
  annotate("text", x = -100, y =22, size = 7,
           label = as.character(expression(paste(Phi,  ": 9.17"))), parse = TRUE) + 
  annotate("text", x = -100, y = 19, size = 7,
           label = paste("P-value (permutation) :",  formatC(race.Phi[3], 4, format = "f")), parse = TRUE)


## category 2
ethnic.indi = ifelse(analysis_dat$PctBlack > 0.10 & analysis_dat$PctHisp > 0.10, 1, 4)
ethnic.indi = ifelse(analysis_dat$PctBlack > 0.10 & analysis_dat$PctHisp <= 0.10, 2, ethnic.indi)
ethnic.indi = ifelse(analysis_dat$PctBlack <= 0.10 & analysis_dat$PctHisp > 0.10, 3, ethnic.indi)
ethnic.indi = ifelse(analysis_dat$PctBlack <= 0.10 & analysis_dat$PctHisp <= 0.10, 4, ethnic.indi)
ethnic.indi = factor(ethnic.indi, levels=c(1, 2, 3, 4), 
                     labels=c(paste(expression(AA>10),"% & ", expression(HP>10),"%"), 
                              paste(expression(AA>10),"% & ", expression(HP<=10),"%"), 
                              paste(expression(AA<=10),"% & ", expression(HP>10),"%"), 
                              paste(expression(AA<=10), "% & ", expression(HP<=10),"%")))
cate2 = gg1 + theme_bw() +  ditch_the_axes + 
  ggtitle("Classification by non-White group") + 
  theme(plot.title = element_text(hjust = 0.5, size = 25)) + 
  geom_point(data = labs, aes(x = long, y = lat, colour = as.factor(ethnic.indi)),size = 2) +
  theme(legend.title=element_blank(),  
        legend.text=element_text(size=20),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(2, 'lines')) +
  guides(color = guide_legend(override.aes = list(size=5))) + 
  annotate("text", x = -100, y =22, size = 7,
           label = as.character(expression(paste(Phi,  ": 22.72"))), parse = TRUE) + 
  annotate("text", x = -100, y = 19, size = 7,
           label = paste("P-value (permutation) :",  formatC(ethnic.indi.Phi[3], 4, format = "f")), parse = TRUE)


pdf("Figure/cateplot.pdf", height = 6, width = 25)
multiplot(cate1, cate2, cols=2)
dev.off()