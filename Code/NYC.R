library(geosphere)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(spdep)
###
source("Code/MoranI.R")
source("Code/Phi.R")
###
nyc = read.csv("Data/NYC_Transit_Subway_Entrance_And_Exit_Data.csv", sep = ",", header = TRUE)
names(nyc)

latitude = nyc$Entrance.Latitude
longitude = nyc$Entrance.Longitude
longitude[which(longitude > 0)] = -longitude[which(longitude > 0)] #correct typo
dist.matrix = distm(cbind(longitude, latitude), cbind(longitude, latitude), fun = distHaversine)
nyc.weight.mat = max(dist.matrix) / dist.matrix
diag(nyc.weight.mat) = 0
nyc.weight.mat = ifelse(nyc.weight.mat > 2000, 2000, nyc.weight.mat)

nyc.entrance.type.Phi = make.permute.Phi(nyc.weight.mat, as.integer(nyc$Entrance.Type), 500)
### map ###
ditch_the_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank()
)
mapnyc = data.frame(latitude, longitude,
                    entrance.type = factor(nyc$Entrance.Type, levels = c("Stair", "Ramp", "Walkway", "Door", "Easement", "Elevator", "Escalator")))
nyc_base = ggmap::get_map("New York City", zoom = 10, 
                          maptype = "terrain-background")
nyc.map = ggmap(nyc_base) + 
  scale_y_continuous(limits = c(min(latitude)-0.05, max(latitude)+0.01), expand = c(0, 0)) +
  scale_x_continuous(limits = c(min(longitude)-0.01, max(longitude)+0.01), expand = c(0, 0)) + 
  theme_bw() +  ditch_the_axes +  
  geom_point(data=mapnyc, aes(x=longitude, y=latitude, color = entrance.type),size= 3) + 
  ggtitle("NYC Transit Subway Entrance Type") + 
  theme(plot.title = element_text(hjust = 0.5, size = 40)) + 
  theme(legend.title= element_blank(),  
        legend.text=element_text(size= 25),
        legend.key = element_rect(size = 50),
        legend.key.size = unit(3, 'lines')) +
  guides(color = guide_legend(override.aes = list(size=5))) + 
  annotate("text", x = -73.90, y = 40.545, size = 10,
           label = as.character(expression(paste(Phi,  ": 1.82"))), parse = TRUE) + 
  annotate("text", x = -73.90, y = 40.53, size = 10,
           label = paste("P-value (permutation) :",  formatC(nyc.entrance.type.Phi[3], 4, format = "f")), parse = TRUE)


pdf("Figure/nycmap.pdf", height = 17, width = 18)
nyc.map
dev.off()  
