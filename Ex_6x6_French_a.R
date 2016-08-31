library(smerc)


## Adjacency matrix 
## making an adjacency matrix 
library(igraph)
ex_admatrix <- matrix(rep(c(0),36*36),nrow = 36,ncol = 36)
set1 <- c(1,2:5,7,13,19,25,8:11,14:17,20:23,26:29) # set take i+1 and i+6
set2 <- c(31,32:35) # set take i+1
set3 <- c(6,12,18,24,30) # set take i+6
set4 <- c(19,21,31,33)

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

## generate population
ex_pop <- rpois(36,1000000)

## generate observed cases for each regions
ex_cases_1a <- c()
region <- c(1:36)

cluster_1a <- c(5,10,11,17)
for ( i in cluster_1a) { ex_cases_1a[i] <- rpois(1,0.005*ex_pop[i])}
cluster_2a <- c(21,19,26,31,33) ## change 21 to 32
for ( i in cluster_2a) { ex_cases_1a[i] <- rpois(1,0.005*ex_pop[i])}

regions_0.2a <- c(4,6,12,16,24,25,23,24,27,15,14)
for ( i in regions_0.2a) { ex_cases_1a[i] <- rpois(1,0.002*ex_pop[i])}

regions_0.1a <- c(3,2,1,7,8,9,13,18,20,34,35,36,30,32,28,29,22)
for ( j in regions_0.1a) { ex_cases_1a[j] <- rpois(1,0.001*ex_pop[j])}

ex_matrix_1 <- data.frame(cbind(ex_cases_1a,ex_pop,region))


## tests for new example
## create distances between 36 regions using x-y coordinates ( look at nydf data
x <- rep((1:6),6) # x coordinate 
y <- c() # y coordinate
for ( i in 1:6) {y <- append(y,rep((i),6))}

cen <- data.frame(cbind(x,y,region))   # region with its x-y coordinates



library(autoimage) ## get images of the density by colors
autoimage(1:6, 1:6, matrix(ex_matrix_1$ex_cases_1a, nrow = 6))


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



## Circular test
coor <- cbind(x,y) # x-y coordinate of each region, distance from 1 to 6. 

test1_exa <- scan.test(coor,floor(ex_matrix_1$ex_cases_1a),
                      ex_matrix_1$ex_pop, nsim = 99,
                      alpha = 0.01, lonlat = FALSE)


## ULS test

test2_exa = uls.test(coor,floor(ex_matrix_1$ex_cases_1a),
                    ex_matrix_1$ex_pop, w = ex_admatrix,
                    alpha = 0.01, lonlat =FALSE,
                    nsim = 99, ubpop = 0.3)

## Flexible test
test3_exa = flex.test(coor,floor(ex_matrix_1$ex_cases_1a),
                     w = ex_admatrix, k = 9,
                     ex_matrix_1$ex_pop, nsim = 99,
                     alpha = 0.01, lonlat = FALSE)
## making up 6x6 grid
## making a graph using igraph
g <- graph.adjacency(ex_admatrix,mode = "Undirected")
plot.igraph(g,vertex.color="yellow")

## clusters of test 1 

cluster_test1a <- c()
for ( i in 1:length(test1_exa$clusters))
{
  cluster_test1a <- append(cluster_test1a,test1_exa$clusters[[i]]$locids)
}

cluster_test2a <- c()
for ( i in 1:length(test2_exa$clusters))
{
  cluster_test2a <-append(cluster_test2a, test2_exa$clusters[[i]]$locids)
}

cluster_test3a <- c()
for ( i in 1:length(test3_exa$clusters))
{
  cluster_test3a <-append(cluster_test3a, test3_exa$clusters[[i]]$locids)
}

plot.igraph(g,vertex.color="yellow")
## plot clusters 
c1a <- cluster_test1a
for ( i in c1)
  
{V(g)[i]$color <- "red"} 


plot(g)

# plot not in igraph
## Circular Method
plot(y ~ x, type = "n", data = cen,main="A Demo of The Circular Method")
textcol_demo = rep("black",36)
# for each radius scanned, the window includes all the regions within that radius
# with radius of 1, the window includes 11,12,17,5
textcol_demo[11] = "green" # radius = 1 
textcol_demo[12] = "blue" 
textcol_demo[17] = "blue" 
textcol_demo[5] = "blue" 
# for radius of 1.41, the window gains 10,16,18
textcol_demo[10] = "blue"# radius = 1.41
textcol_demo[16] = "red" 
textcol_demo[18] = "red" 
textcol_demo[4] = "red"
textcol_demo[6] = "red"

text(y ~ x, lab = 1:36, data = cen, col = textcol_demo)


plot(y ~ x, type = "n", data = cen,main="Example 2: Clusters of The Circular Method")

textcol = rep("black", 36)
textcol[cluster_test1a] = "blue"

text(y ~ x, lab = 1:36, data = cen, col = textcol)

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
mycol[cluster_test1a] = "blue"
#mycol[c(5, 6, 10, 11, 12, 17)] = "orange"

plot(sppoly, col = mycol, main="Clusters of The Circular Method")
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons












## ULS Method

ratio_a <- c()
for ( i in 1:36)
{
  ratio_a[i] <- ex_cases_1a[i]/ex_pop[i]
}
ratio_a_list <- c()
ratio_a_list <- append(ratio_a_list, which(ratio_a >= sum(ex_cases_1a)/sum(ex_pop)))
plot(y ~ x, type = "n", data = cen,main="A Demo of The ULS Method")
textcol1_demo = rep("black", 36)
textcol1_demo[ratio_a_list] = "orange"
text(y ~ x, lab = 1:36, data = cen, col = textcol1_demo)


plot(y ~ x, type = "n", data = cen,main="Example 2: Clusters of The ULS Method")
textcol1 = rep("black", 36)
textcol1[cluster_test2a] = "red"

text(y ~ x, lab = 1:36, data = cen, col = textcol1)
# new plots from French
# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[cluster_test2a] = "yellow"


plot(sppoly, col = mycol, main="Clusters of The ULS Method")
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons





## Flexible Method 
## For each time, the window pick one region that is the nearest to the original centroid
plot(y ~ x, type = "n", data = cen,main="A Demo of The Flexible Method,k=7")
textcol2_demo = rep("black", 36)
textcol2_demo[11] = "green"  
textcol2_demo[12] = "blue" 
textcol2_demo[17] = "blue" 
textcol2_demo[5] = "blue" 
textcol2_demo[10] = "blue"
textcol2_demo[16] = "red" 
textcol2_demo[18] = "red" 
# in this scenario, k = 7  
text(y ~ x, lab = 1:36, data = cen, col = textcol2_demo)


## plot the clusters
plot(y ~ x, type = "n", data = cen,main="Example 2: Clusters of The Flexible Method,k=9")
textcol2 = rep("black", 36)
textcol2[cluster_test3a] = "green"

text(y ~ x, lab = 1:36, data = cen, col = textcol2)
## plot the true clusters
plot(y ~ x, type = "n", data = cen,main="Example 2: True Clusters")
textcol2 = rep("black", 36)
textcol2[clusters_1a] = "green"
textcol2[cluster_2a] = "green"

text(y ~ x, lab = 1:36, data = cen, col = textcol2)
# new plots from French
# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[cluster_test3a] = "red"


plot(sppoly, col = mycol, main="Clusters of The Flexible Method")
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons



# true cluster plot
# color certain polygons
mycol = rep("white", nloc)
# change color of true clusters to a different color
mycol[cluster_1a] = "green"
mycol[cluster_2a] = "green"

plot(sppoly, col = mycol, main="True Clusters")
text(coordinates(sppoly), labels=row.names(sppoly)) # name polygons
