library(smerc)

## generate population
ex_pop <- rpois(36,1000000)

## generate observed cases for each regions
ex_cases_1 <- c()
region <- c(1:36)

clusters_1 <- c(11,12,17,6)
for ( i in clusters_1) { ex_cases_1[i] <- rpois(1,0.005*ex_pop_kept_1[i])}
cluster_2 <- c(32,20,26) ## change 21 to 32
for ( i in cluster_2) { ex_cases_1[i] <- rpois(1,0.005*ex_pop_kept_1[i])}

regions_0.2 <- c(4,5,10,18,21,24,23,27,15,14)
for ( i in regions_0.2) { ex_cases_1[i] <- rpois(1,0.002*ex_pop_kept_1[i])}

regions_0.1 <- c(3,2,1,7,8,9,13,19,16,25,31,33,34,35,36,30,28,29,22)
for ( j in regions_0.1) { ex_cases_1[j] <- rpois(1,0.001*ex_pop_kept_1[j])}

ex_matrix_1 <- data.frame(cbind(ex_cases_1,ex_pop_kept,region))


## tests for new example
## create distances between 36 regions using x-y coordinates ( look at nydf data
x <- rep((1:6),6) # x coordinate 
y <- c() # y coordinate
for ( i in 1:6) {y <- append(y,rep((i),6))}

cen <- data.frame(cbind(x,y,region))   # region with its x-y coordinates



library(autoimage) ## get images of the density by colors
autoimage(1:6, 1:6, matrix(ex_matrix_1$ex_cases_1, nrow = 6))


## Nearest neigbors 
library(FNN)
# distance matrix 

dis <- data.frame() ## distance matrix 

for ( i in 1:36)
{
  for ( j in 1:36)
    
  {
    dis[i,i] <- 0
    dis[i,j] <- sqrt((cen[i,]$y - cen[j,]$y)^2 + (cen[i,]$x - cen[j,]$x)^2)
  }
}

## Adjacency matrix 
## making an adjacency matrix 
library(igraph)
ex_admatrix <- matrix(rep(c(0),36*36),nrow = 36,ncol = 36)
set1 <- c(1,2:5,7,13,19,25,8:11,14:17,20:23,26:29) # set take i+1 and i+6
set2 <- c(31,32:35) # set take i+1
set3 <- c(6,12,18,24,30) # set take i+6

for ( i in set1)
{
  ex_admatrix[i,i+1] <- 1 
  ex_admatrix[i,i+6] <- 1
  
  for ( i in set2)
  {
    ex_admatrix[i,i+1] <- 1 
    
    
    for ( i in set3)
    {
      
      ex_admatrix[i,i+6] <- 1
    }
  }}
for ( i in 1:36){ ex_admatrix[,i] <- ex_admatrix[i,]}
for ( i in set4) # add 4 more neighbors to 26 
{
  ex_admatrix[26,i] <- 1
  ex_admatrix[i,26] <- 1
}


## Circular test
coor <- cbind(x,y) # x-y coordinate of each region, distance from 1 to 6. 

test1_ex <- scan.test(coor,floor(ex_matrix_1$ex_cases_1),
                      ex_matrix_1$ex_pop_kept_1, nsim = 99,
                      alpha = 0.01, lonlat = FALSE)


## ULS test

test2_ex = uls.test(coor,floor(ex_matrix_1$ex_cases_1),
                    ex_matrix_1$ex_pop_kept_1, w = ex_admatrix,
                    alpha = 0.01, lonlat = TRUE,
                    nsim = 99, ubpop = 0.2)

## Flexible test
test3_ex = flex.test(coor,floor(ex_matrix_1$ex_cases_1),
                     w = ex_admatrix, k = 4,
                     ex_matrix_1$ex_pop_kept_1, nsim = 99,
                     alpha = 0.01, lonlat = TRUE)
## making up 6x6 grid
## making a graph using igraph
g <- graph.adjacency(ex_admatrix,mode = "Undirected")
plot.igraph(g,vertex.color="yellow")

## clusters of test 1 

cluster_test1 <- c()
for ( i in 1:length(test1_ex$clusters))
{
  cluster_test1 <- append(cluster_test1,test1_ex$clusters[[i]]$locids)
}

cluster_test2 <- c()
for ( i in 1:length(test2_ex$clusters))
{
  cluster_test2 <-append(cluster_test2, test2_ex$clusters[[i]]$locids)
}

cluster_test3 <- c()
for ( i in 1:length(test3_ex$clusters))
{
  cluster_test3 <-append(cluster_test3, test3_ex$clusters[[i]]$locids)
}

plot.igraph(g,vertex.color="yellow")
## plot clusters 

## plot Circular 
## new plots from French

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
polydf = SpatialPolygonsDataFrame(sppoly, data = data.frame(pop = rpois(36, 1000000)))

# plot polygons
plot(sppoly)

# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[cluster_test1] = "blue"
#mycol[c(5, 6, 10, 11, 12, 17)] = "orange"

plot(sppoly, col = mycol, main="Clusters of The Circular Scanning Test")
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons




## Plot the ULS 
## new plots from French

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
polydf = SpatialPolygonsDataFrame(sppoly, data = data.frame(pop = rpois(36, 1000000)))

# plot polygons
plot(sppoly)

# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[cluster_test2] = "blue"
#mycol[c(5, 6, 10, 11, 12, 17)] = "orange"

plot(sppoly, col = mycol, main="Clusters of The ULS Scanning Test")
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons

## Plot for the Flexible 
## new plots from French

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
polydf = SpatialPolygonsDataFrame(sppoly, data = data.frame(pop = rpois(36, 1000000)))

# plot polygons
plot(sppoly)

# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[cluster_test3] = "blue"
#mycol[c(5, 6, 10, 11, 12, 17)] = "orange"

plot(sppoly, col = mycol, main="Clusters of Flexible Scanning Test")
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons
