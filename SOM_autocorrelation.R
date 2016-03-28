library(sp)
library(class)
library(kohonen)
library(maptools)
library(reshape2)
library(rgeos)
library(rgdal)
library(vegan)

l_m_entry = readShapePoly("F:/Master thesis/Data/New_Lisboa/Lisboa_clusters_entry1.shp")
l_m_mixed = readShapePoly("F:/Master thesis/Data/New_Lisboa/Lisboa_clusters_mixed.shp")
l_m_population = readShapePoly("F:/Master thesis/Data/New_Lisboa/Lisboa_clusters_population.shp")
l_m_household = readShapePoly("F:/Master thesis/Data/New_Lisboa/Lisboa_clusters_houshold1.shp")

d_entry = l_m_entry@data
d_mixed = l_m_mixed@data
d_pop = l_m_population@data
d_houshold = l_m_household@data
#breaks= dataset$cluster 
#colors <- c("red", "blue", "lightpink", "skyblue2", "white","black","yellow")
#plot(lisbon_map, col = colors[breaks])#colors[np] manually set the color for each region
#mtext("Map", cex=1.5, side = 3, line = 1)

# SOM_Entry test gives 0.012 p-value
w=(dist(d_entry[,c(d_entry$x,d_entry$y)]))
u=(dist(d_entry$cluster))
mantel(w,u)
# SOM_Mixed test gives 0.011 p-value
w=(dist(d_entry[,c(d_entry$x,d_entry$y)]))
u=(dist(d_mixed$cluster))
mantel(w,u)
# SOM_Population test gives 0.23 p-value
w=(dist(d_entry[,c(d_entry$x,d_entry$y)]))
u=(dist(d_pop$cluster))
mantel(w,u)
# SOM_Household test gives 0.838 p-value
w=(dist(d_entry[,c(d_entry$x,d_entry$y)]))
u=(dist(d_houshold$cluster))
mantel(w,u)


