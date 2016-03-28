# making funcitons for cut off duplicates
x = read.csv(file = "F:/Master thesis/Data/Original_Tourist_here_maps.csv")
duplicat = function(y){
  b = x[x$feature_type==y,]
  subset(b,!duplicated(b[3:4]))
}
z = read.csv(file = "F:/Master thesis/Data/restaurants_POI_Rest_Coffe_Lx.csv")

duplicat1 = function(y){
  k = z[z$facility_type_desc==y,]
  subset(k,!duplicated(k[8:9]))
}
         
