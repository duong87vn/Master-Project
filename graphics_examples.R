library(sp)
library(spdep)
library(autoimage)
library(data.table)
library(maptools)
library(fields)

# data.frame containging all grid locations IN TERMS OF POSITION
g = expand.grid(x = 1:6, y = 1:6)

# a function to create a polygon for each grid point
create_polys = function(x){
  xlocs = x[1] + c(-0.5, -0.5, 0.5, 0.5, -0.5)
  ylocs = x[2] + c(-0.5, 0.5, 0.5, -0.5, -0.5)
  data.frame(x = c(xlocs, NA), y = c(ylocs, NA))
}

# create polygon for grid locations
polys = apply(g, 1, create_polys)
polys_df = data.table::rbindlist(polys)
polys_list = list(x = polys_df$x, y = polys_df$y) # convert data frame to list format

# number of regions
nloc = length(polys)

# color certain polygons
sppoly <- map2SpatialPolygons(polys_list, IDs = seq_len(nloc))

# create SpatialPolygonsDataFrame

# replace rpois(36, 1000000) with the populations you generated previously
polydf = SpatialPolygonsDataFrame(sppoly, data = data.frame(pop = ex_pop_kept))

# plot polygons
plot(sppoly)

# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[c(20, 26, 32)] = "blue"
mycol[c(5, 6, 10, 11, 12, 17)] = "orange"

plot(sppoly, col = mycol)
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons

# plot of population sizes
textlab = list("sp.text", coordinates(sppoly), row.names(sppoly)) # to label regions
# to plot regions with colors corresponding to population size
spplot(polydf, main = "population size", col.regions = tim.colors(), sp.layout = list(textlab))

# determine neighborhood structure
kn = poly2nb(sppoly, queen = FALSE)
listw <- nb2listw(kn, style = "B")

# you would change this to something like:
listw <- mat2listw(ex_admatrix) #where x is the matrix specifying the neighbors of each region
plot(listw, coordinates(sppoly))
text(coordinates(sppoly), labels=row.names(sppoly))


## Demonstration
# Flexible
# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[c(11)] = "green"
mycol[c(5,10,12,17)] = "blue"
mycol[c(16)] = "orange"

plot(sppoly, col = mycol)
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons

#ULS
# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[c(21,26,31,33,19)] = "green"
mycol[c(5,10,17)] = "blue"
mycol[c(11)] = "orange"

plot(sppoly, col = mycol)
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons

#Circular
# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[c(4,6,16,18)] = "green"
mycol[c(5,10,17,12)] = "blue"
mycol[c(11)] = "orange"

plot(sppoly, col = mycol)
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons
