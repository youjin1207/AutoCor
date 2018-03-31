load("Data/analysis_dat.RData")
names(analysis_dat)
# small sized : power plants consisted by 1 or 2 EGUs
# medium sized : power plants consistec by 3,4, or 5 EGUs
# EGU : energy generating unit
size.indicator = rep(1, nrow(analysis_dat))
size.indicator = ifelse(analysis_dat$small_nunits == 1, 2, size.indicator)
size.indicator = ifelse(analysis_dat$med_nunits == 1, 3, size.indicator)

# locations
library(geosphere)
latitude = analysis_dat$Fac.Latitude
longitude = analysis_dat$Fac.Longitude
dist.matrix = distm(cbind(longitude, latitude), cbind(longitude, latitude), fun = distHaversine)
weights = max(dist.matrix) / dist.matrix
diag(weights) = 0
summary(as.numeric(weights))
weights = ifelse(weights > 10, 10, weights)