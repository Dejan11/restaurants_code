library(sp)
library(class)
library(kohonen)
library(dummies)
library(ggplot2)
library(maptools)
library(reshape2)
library(rgeos)
library(rgdal)
library(ggmap)
library(leaflet)
library(RColorBrewer)
library(corrplot)
library(ggdendro)
library(vegan)

source("C:/GEOTECH/NOVA/GPS/2014-01-SOM-Example-code_release/coolBlueHotRed.R")
# source(file = "F:/Master thesis/Data/thesis_files.R")
lisbon_map = readShapePoly("F:/Master thesis/Data/New_Lisboa/Final_Lisboa12.shp")
results = data.frame(lisbon_map)
# lm = lisbon_map
spearman = cor(results[-1], method = "spearman")
pearson = cor(results[-1], method = "pearson")
corrplot(spearman, type = "upper", order="hclust", tl.col="black", tl.srt=60,title = "Spearman's correlation",mar = c(0,0,1,0), cl.cex = 1.5, tl.cex = 1 )
corrplot(pearson, type = "upper", order="hclust", tl.col="black", tl.srt=60,main = "Pearson's correlation",mar = c(0,0,1,0),cl.cex = 1.5, tl.cex = 1)

color_palette = colors()[c(26,36,254,552,176,394,261)]
som_entry = results[c(-1:-3,-5,-6, -9,-11,-15,-16, -21,-22) ]

# data_som_matrix = scale(data_som,center = FALSE, scale = apply(data_som,2,sd, na.rm=TRUE))
x = som_entry
data_som_matrix = apply(som_entry,2,function(x)(x-min(x))/(max(x)-min(x)))

names(data_som_matrix) <- names(som_entry)
# Execute SOM
set.seed(7)
som_model = som(data_som_matrix, 
                grid=somgrid( xdim = 16, ydim = 16, topo = "hexagonal" ), 
                rlen=100, 
                alpha=c(0.05,0.01),
                toroidal = TRUE,
                n.hood = "circular",
                keep.data = TRUE,
                init = data_som_matrix[seq(4,1024,4),]
)

par(mar = c(2, 2, 2, 2))
plot(som_model, main = "Platform")

# list of two matrices, containing codebook vectors for X and Y, respectively(shows the codebook vectors).
plot(som_model, type = "codes")
# shows the mean distance to the closest codebook vector during training.
plot(som_model, type = "changes")
#counts nodes shows the number of objects mapped to the individual units. Empty units are depicted in gray
plot(som_model, type = "counts", main="Node Counts - No of objects per unit", palette.name=coolBlueHotRed)
#map quality
# shows the sum of the distances to all immediate neighbours.  This kind of visualisation is also known as a U-matrix plot.  Units near a class boundary can be expected to
# have higher average distances to their neighbours.  Only available for the "som" and "super-
# som" maps, for the moment.
plot(som_model, type="dist.neighbours", main = "SOM neighbour distances", palette.name=grey.colors)
#shows the mean distance of objects mapped to a unit to the codebook vector of that unit.
# The smaller the distances, the better the objects are represented by the codebook vectors.
plot(som_model, type = "quality", main="Node Quality- Mean Distance", palette.name=coolBlueHotRed)

# properties of each unit can be calculated and shown in colour code.  It can be used to visualise the similarity of one particular object to all units in the map, to show the mean
# similarity of all units and the objects mapped to them,  etcetera.   The parameter property contains the numerical values
par(mfrow=c(4,4))
for (i in 1:14) 
  plot(som_model, type = "property", property = som_model$codes[,i], main=names(som_model$data)[i], palette.name=grey.colors)
par(mfrow=c(1,1))
#Clustering SOM - som_model$codes - values for each cell - as a result from SOM algorithm
# wss significates within cluster sum of squares, hence first for each variable (column) we are calculating variance, 
# then summary of them multiplied with the number of cells (12*20 - 1) 
mydata_clust <- som_model$codes
# within cluster sum of squares in case that is only one centroid. wcss is square summary of distances from 
# entity point to cluster centroid
wcss = (nrow(mydata_clust)-1)*sum(apply(mydata_clust,2,var))

# The data given by x are clustered by the k-means method, which aims to partition the points into k groups such that the sum of squares from points 
# to the assigned cluster centres is minimized. At the minimum, all cluster centres are at the mean of their Voronoi 
# sets (the set of data points which are nearest to the cluster centre).
# within cluster sum of squares for 2 up to 8 possible clusters, where tot.withinss is summary of sum of squars
for (i in 2:8) wcss[i] = kmeans(mydata_clust, centers=i)$tot.withinss
# adjusting plot location
par(mar=c(5.1,4.1,4.1,2.1))
par(mfcol=c(1,1), mar=c(3,4,1,0.5), oma=c(0,0,1,0))
# WCSS and number of clusters
plot(wcss, xlab="Number of Clusters", ylab="Within groups sum of squares", main="Within cluster sum of squares (WCSS)")
# Form clusters on grid
## use hierarchical clustering to cluster the codebook vectors

som_cluster <- cutree(hclust(dist(som_model$codes)), 7)
dendo = hclust(dist(som_model$codes))
dhc <- as.dendrogram(dendo)
# Rectangular lines
ddata <- dendro_data(dhc, type = "triangle")
p <- ggplot(segment(ddata)) + ggtitle("Dendogram branches - SOM codebook values")+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) 
p
ddata <- dendro_data(dhc, type = "rectangle")
p <- ggplot(segment(ddata)) + ggtitle("Dendogram branches - SOM codebook values")+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0))
p
# Show the map with different colours for every cluster
color_palette = colors()[c(176,254,26,552,36,585,394,261)]
plot(som_model, type="mapping", bgcol = color_palette[som_cluster], main = "Clusters")
add.cluster.boundaries(som_model, som_cluster)
#show the same plot with the codes instead of just colours
plot(som_model, type= "codes", bgcol = color_palette[som_cluster], main = "Clusters")
add.cluster.boundaries(som_model, som_cluster)
# Plot the map of Lisbon, coloured by the clusters the map to show locations
cluster_details = data.frame(id=results$seccao, cluster = som_cluster[som_model$unit.classif])
# Adding clusters to each city section
lisbon_map@data[26] = seq(1:1053)
cluster_details[3] = seq(1:1053)
rr2 = merge(lisbon_map@data,cluster_details, by.x = "V26", by.y = "V3")
lisbon_map@data = subset(rr2, select = c(-1,-27))
cluster_details[3]= NULL
# saving shapefile
#writePolyShape(lisbon_map,"F:/Master thesis/Data/New_Lisboa/Lisboa_clusters_entry")
# This part is related to GGPLOT, we have to make data.frame from spatial data frame
# it's possible to fortify by "cluster" and then in fill=factor(id) we would have results of cluster, but no borders(no need to merge)
lisbon = fortify(lisbon_map, region="seccao")
# merging clusters with lisbon data.frame by id's 
lisbon = merge(lisbon, cluster_details, by="id")
g = ggplot(lisbon) + aes(long, lat, group=group, fill=factor(cluster)) + 
	geom_polygon() + geom_path( color = "white") + 
	coord_equal() + scale_fill_manual(values = color_palette, breaks=c(1,2,3,4,5,6,7,8), labels=c("Random - residential","Tax-Houshold 1&2","Older", "Low tax - residential", "Low density","No residential - red", "Higher class", "NULL")) + 
	ggtitle("The map of Lisbon - Clusters with entry variables") + theme(legend.key = element_rect(colour = "black"), legend.title = element_text(face = "italic")) + 
	guides(fill = guide_legend(title = "CLUSTERS", title.position = "top", label.position = "bottom")) 
plot(g)
# spatial autocorrelation test for spatial clustering
l_m = readShapePoly("F:/Master thesis/Data/New_Lisboa/Lisbon_with_centroids.shp")
l_m = data.frame(l_m)
x = apply(som_model$codes,1,sum)
y = x[som_model$unit.classif]
l_m$value = y - som_model$distances
# mantel test for spatial autocorrelation
# SOM_Entry test gives 0.012 p-value
w=(dist(l_m[,c(l_m$x,l_m$y)]))
u=(dist(l_m$value))
mantel(w,u)
mantel.correlog(u, w)
