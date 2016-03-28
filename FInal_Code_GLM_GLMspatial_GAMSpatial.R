library(sp)
library(class)
library(ggplot2)
library(maptools)
library(reshape2)
library(rgeos)
library(spdep)
library(spacemakeR) 
library(tripack)
library(mgcv)
library(MASS)
require(Rcmdr)
library(pscl)
library(cluster)
library(ggdendro)
library(eeptools)
library(stats)
library(rgdal)
library(GISTools)

l_m = readShapePoly("F:/Master thesis/Data/New_Lisboa/Lisbon_with_centroids.shp")
class(l_m)
names(l_m)
dataset=data.frame(l_m)

xy=cbind(l_m$x,l_m$y)
plot(l_m,border="gray")
points(xy,pch="+",col="blue")

#==========================#
###EXPLANORATORY ANALYSIS###
#==========================#
# Histogram
#===========================
#histomap(l_m,"Restaurant")####be careful how you use
#===========================
# BoxMap
#===========================

#boxplotmap(l_m,"z")####be careful how you use

#====================================================
# Physical contiguity criteria for irregular l_m
#====================================================
op=par(mfrow=c(1,2))
rooknb1=poly2nb(l_m,queen=FALSE)
plot(l_m,border="gray")
plot(rooknb1,xy,add=T,col="blue")
title(main="Rook")
queennb1=poly2nb(l_m,queen=TRUE)
plot(l_m,border="gray")
plot(queennb1,xy,add=T,col="blue")
title(main="Queen")
par(op)

#==============================
# Criterions based on graphycs
#==============================
op=par(mfrow=c(2,2))
trinb=tri2nb(xy)
plot(l_m,border="gray")
plot(trinb,xy,add=T,col="blue")
title(main="Triangulation Delaunay")
soinb=graph2nb(soi.graph(trinb,xy))
plot(l_m,border="gray")
plot(soinb,xy,add=T,col="blue")
title(main="Sphere of influence")
gabrielnb=graph2nb(gabrielneigh(xy),sym=TRUE)
plot(l_m,border="gray")
plot(gabrielnb,xy,add=T,col="blue")
title(main="Gabriel Graphics")
relativenb=graph2nb(relativeneigh(xy),sym=TRUE)
plot(l_m,border="gray")
plot(relativenb,xy,add=T,col="blue")
title(main="Relative neighboors")
par(op)
#=============================
# Criterions based on distance
#=============================
op=par(mfrow=c(2,2))
knn1_nb=knn2nb(knearneigh(xy, k = 1))
plot(l_m,border="gray")
plot(knn1_nb,xy,add=T,col="blue")
title(main="Nearest neighbor")
knn2_nb=knn2nb(knearneigh(xy, k = 2))
plot(l_m,border="gray")
plot(knn2_nb,xy,add=T,col="blue")
title(main="2 Nearest neighbors")
knn3_nb=knn2nb(knearneigh(xy, k = 3))
plot(l_m,border="gray")
plot(knn3_nb,xy,add=T,col="blue")
title(main="3 Nearest neighbors")
knn4_nb=knn2nb(knearneigh(xy, k = 4))
plot(l_m,border="gray")
plot(knn4_nb,xy,add=T,col="blue")
title(main="4 Nearest neighbors")
par(op)

#==========================================================================
### SELECt matrix based on neighborhood by (PrincipaCoordinatesNeighMatrix) 
#==========================================================================
summary(test.W(l_m$Restaurant,rooknb1))#testing which criterion is better by AIC- select best of them with lower AIC
summary(test.W(l_m$Restaurant,queennb1))#physical criterion with queen
summary(test.W(l_m$Restaurant,trinb))#delaneuy criterion
summary(test.W(l_m$Restaurant,soinb))#sphere of influence
summary(test.W(l_m$Restaurant,gabrielnb))#gabriel graphycs
summary(test.W(l_m$Restaurant,relativenb))#Relative neighbors
summary(test.W(l_m$Restaurant,knn1_nb))
summary(test.W(l_m$Restaurant,knn2_nb))
summary(test.W(l_m$Restaurant,knn3_nb))
summary(test.W(l_m$Restaurant,knn4_nb))

#============================#
###Spatial weights matrices###
#============================#
knn2W1=nb2listw(knn2_nb,style="W")
summary(knn2W1)
#================================================#
###Global indicators - spatial autocorrelation ###
#================================================#
moran.plot(as.vector(scale(l_m$Restaurant)),knn2W1,xlab="Restaurants",ylab="No of rest of neighbor's cell ",main="KNN 2")
moran.test(l_m$Restaurant,knn2W1,alternative="two.sided")
geary.test(l_m$Restaurant,knn2W1,alternative="two.sided")
#===============================================#
###Local indicators - spatial autocorrelation ###
#===============================================#
locm3=localmoran(l_m$Restaurant,knn2W1,alternative="two.sided")#local spatial autocorrelation
#========#
###MAPS###
#========#
l_m$sz <- scale(l_m$Restaurant)#sz = scale for no of restarurants - 0 is new mean and new 1 is sd 
l_m$lag_sz <- lag.listw(knn2W1,l_m$sz)#lag_sz is the same for neighbors of the number of restaurants
l_m$quad_sig <- NA
l_m@data[(l_m$sz >= 0 & l_m$lag_sz >= 0) & (locm3[, 5] <= 0.05), "quad_sig"] <- 1
l_m@data[(l_m$sz <= 0 & l_m$lag_sz <= 0) & (locm3[, 5] <= 0.05), "quad_sig"] <- 2
l_m@data[(l_m$sz >= 0 & l_m$lag_sz <= 0) & (locm3[, 5] <= 0.05), "quad_sig"] <- 3
l_m@data[(l_m$sz >= 0 & l_m$lag_sz <= 0) & (locm3[, 5] <= 0.05), "quad_sig"] <- 4
l_m@data[(l_m$sz <= 0 & l_m$lag_sz >= 0) & (locm3[, 5] <= 0.05), "quad_sig"] <- 5

breaks=seq(1, 5, 1)
labels=c("High-High", "Low-Low", "High-Low", "Low-High", "Not Signif.")
np <- findInterval(l_m$quad_sig, breaks)
colors <- c("red", "blue", "lightpink", "skyblue2", "white")
plot(l_m, col = colors[np])#colors[np] manually set the color for each region
mtext("Local neighborhoods of Moran test", cex=1.5, side = 3, line = 1)
legend("bottomright", legend = labels, fill = colors, bty = "n",cex=0.7, inset = c(0.05,0))
#=======================================================#
###Development of (spatial)GLM and (spatial)GAM models###
#=======================================================#
###Before the models have been created, the dataset is split it into 80% of training and 20% of validating
#20% at the end is used for testing and the best model is chosen to predict the rest of the values
names(dataset)
#==========================
#Full GLM model
#==========================
GLM.1 <- glm(Restaurant ~ BuildAft01 + Colle_deg + Employed + Exlus_resi + Hous1or2 + Hous3or4 + HousNoUnem + Men_20_24 + Men_25_64 + 
               Men_64_ + Pensioners + Res_No_act + Tax_index + Tourst_idx + Woman_64_ + Womn_20_24 + Womn_25_64 + Work_in_Te + Work_Sec, 
             family=poisson(log), data=dataset)
summary(GLM.1)
#=============================
# Exponentiated coefficients
#=============================
exp(coef(GLM.1))
#============================================================
#Stepwise model selection based on backward/forward selection
#============================================================
#AIC - finds the model that gives the best prediction - for model comparison
#BIC - assumes that one of the models is the true model and find the "true" model - variable selection
#============================================================
#Creation of the base GLM model
#============================================================
GLM.2=stepwise(GLM.1, direction='forward/backward', criterion='BIC')
summary(GLM.2, cor=FALSE)
cor(dataset$Restaurant,predict(GLM.2, type = "response"))
#============================================================
#evaluation models
#============================================================
AIC(GLM.2, GLM.1)
anova(GLM.2,GLM.1, test = "Chisq") #anova proves it that reduced model is good
#=================================================================================
#Developing spatial component with Moran eigenvector and spatial weights from KNN2
#=================================================================================
GLM.3 <- ME(Restaurant ~ Hous1or2+HousNoUnem+Indv20_24+Indv25_64+Post_sec_e+Colle_deg+Res_No_act+Res_20to64+Work_in_Te+Work_Sec+Tax_index+Tourst_idx+Womn_20_24+Men_20_24+Men_25_64+Womn_25_64+Woman_64_+Men_64_+Pensioners+Employed+BuildAft01+Hous3or4, data=dataset, family="poisson", listw=knn2W1, alpha=0.5)
GLM.4 <- glm(Restaurant ~ Hous1or2+HousNoUnem+Indv20_24+Indv25_64+Post_sec_e+Colle_deg+Res_No_act+Res_20to64+Work_in_Te+Work_Sec+Tax_index+Tourst_idx+Womn_20_24+Men_20_24+Men_25_64+Womn_25_64+Woman_64_+Men_64_+Pensioners+Employed+BuildAft01+Hous3or4+fitted(GLM.3), data=dataset,family=poisson(log))
summary(GLM.4)
cor(dataset$Restaurant,predict(GLM.4, type = "response"))
#=================================================================================
#Involving spatial component into base GLM model
#=================================================================================
GLM.6 <- glm(Restaurant ~ Tourst_idx+Tax_index+Hous3or4+Exlus_resi+Men_25_64+Woman_64_+Men_64_+Work_Sec+ fitted(GLM.3),data = dataset,family = poisson(log))

#GLM.5=stepwise(GLM.4, direction='forward/backward', criterion='AIC')
summary(GLM.6, cor=FALSE)
AIC(GLM.2,GLM.6)
anova(GLM.6, GLM.2, test = "Chisq")
anova(GLM.6)
cor(dataset$Restaurant,predict(GLM.6, type = "response"))
shapiro.test(residuals(GLM.6, type = "deviance"))
#============================================================
#Creation of base GAM model first step
#============================================================

lGAMp <- gam(Restaurant~s(Hous1or2, bs = "cr"),data=dataset,family="poisson")
lGAMp1<- gam(Restaurant~s(Tax_index, bs = "cr"),data=dataset,family="poisson")#+HousNoUnem+Indv20_24+Indv25_64+Post_sec_e+Colle_deg+Res_No_act+Res_20to64+Work_in_Te+Work_Sec+Tax_index+Tourst_idx+Womn_20_24+Men_20_24+Men_25_64+Womn_25_64+Woman_64_+Men_64_+Pensioners+Employed+BuildAft01+Hous3or4+s(x,y),data=dataset,family=poisson(log))
lGAMp2<- gam(Restaurant~s(I(Tourst_idx^0.5)),data=dataset,family="poisson")#+HousNoUnem+Indv20_24+Indv25_64+Post_sec_e+Colle_deg+Res_No_act+Res_20to64+Work_in_Te+Work_Sec+Tax_index+Tourst_idx+Womn_20_24+Men_20_24+Men_25_64+Womn_25_64+Woman_64_+Men_64_+Pensioners+Employed+BuildAft01+Hous3or4+s(x,y),data=dataset,family=poisson(log))
lGAMp3<- gam(Restaurant~s(Work_in_Te, bs = "cr"),data=dataset,family="poisson")
lGAMp4<- gam(Restaurant~s(Colle_deg),data=dataset,family="poisson")
lGAMp5<- gam(Restaurant~s((Employed), bs = "cr"),data=dataset,family="poisson")
lGAMp6<- gam(Restaurant~s(I(Womn_20_24^0.5)),data=dataset,family="poisson")
lGAMp7<- gam(Restaurant~s(Res_No_act),data=dataset,family="poisson")
lGAMp8<- gam(Restaurant~s((Woman_64_)),data=dataset,family="poisson")
lGAMp9<- gam(Restaurant~s((Indv25_64), bs = "cr"),data=dataset,family="poisson")

BICall = BIC(lGAMp,lGAMp1,lGAMp2,lGAMp3,lGAMp4,lGAMp5,lGAMp6,lGAMp7,lGAMp8,lGAMp9)
which.min(BICall[,2])
#============================================================
#Creation of base GAM model - 2nd step
#============================================================
lGAMp10<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Hous1or2, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp11<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp12<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Work_in_Te, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp13<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Colle_deg),data=dataset,family="poisson", method = "REML", select = T)
lGAMp14<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Employed, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp15<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(I(Womn_20_24^0.5)),data=dataset,family="poisson", method = "REML", select = T)
lGAMp16<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Res_No_act),data=dataset,family="poisson", method = "REML", select = T)
lGAMp17<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Indv25_64, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp18<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Woman_64_),data=dataset,family="poisson", method = "REML", select = T)

BICall1 = BIC(lGAMp10,lGAMp11,lGAMp12,lGAMp13,lGAMp14,lGAMp15,lGAMp16,lGAMp17,lGAMp18)
which.min(BICall1[,2])
#============================================================
#Creation of base GAM model - 3rd step
#============================================================
lGAMp19<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp21<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Colle_deg),data=dataset,family="poisson", method = "REML", select = T)
lGAMp22<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Employed, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp23<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(I(Womn_20_24^0.5)),data=dataset,family="poisson", method = "REML", select = T)
lGAMp24<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Res_No_act),data=dataset,family="poisson", method = "REML", select = T)
lGAMp25<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Indv25_64, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp26<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Woman_64_, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
BICall2 = BIC(lGAMp19,lGAMp21,lGAMp22,lGAMp23,lGAMp24,lGAMp25,lGAMp26)
which.min(BICall2[,2])
#============================================================
#Creation of base GAM model - 4th step
#============================================================
lGAMp27<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg),data=dataset,family="poisson", method = "REML", select = T)
lGAMp28<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+s(Employed, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp29<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+s(I(Womn_20_24^0.5)),data=dataset,family="poisson", method = "REML", select = T)
lGAMp30<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+s(Res_No_act),data=dataset,family="poisson", method = "REML", select = T)
lGAMp31<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+s(Indv25_64, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp32<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+s(Woman_64_, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
BICall3 = BIC(lGAMp27,lGAMp28,lGAMp29,lGAMp30,lGAMp31,lGAMp32)
which.min(BICall3[,2])
lGAMp33<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(Employed, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp34<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5)),data=dataset,family="poisson", method = "REML", select = T)
lGAMp35<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(Res_No_act),data=dataset,family="poisson", method = "REML", select = T)
lGAMp36<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(Indv25_64, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp37<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(Woman_64_, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
BICall4 = BIC(lGAMp33,lGAMp34,lGAMp35,lGAMp36,lGAMp37)
which.min(BICall4[,2])
lGAMp38<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp39<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Res_No_act),data=dataset,family="poisson", method = "REML", select = T)
lGAMp40<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Indv25_64, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp41<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Woman_64_, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
BICall5 = BIC(lGAMp38,lGAMp39,lGAMp40,lGAMp41)
which.min(BICall5[,2])
lGAMp42<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+s(Res_No_act),data=dataset,family="poisson", method = "REML", select = T)
lGAMp43<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+s(Indv25_64, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
lGAMp44<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+s(Woman_64_, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
BICall6 = BIC(lGAMp42,lGAMp43,lGAMp44)
which.min(BICall6[,2])
lGAMp45<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+s(Indv25_64, bs = "cr")+s(Res_No_act),data=dataset,family="poisson", method = "REML", select = T)
lGAMp46<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+s(Indv25_64, bs = "cr")+s(Woman_64_, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
BICall7 = BIC(lGAMp45,lGAMp46)
which.min(BICall7[,2])
lGAMp47<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+s(Indv25_64, bs = "cr")+s(Res_No_act)+s(Woman_64_, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)

BICall8 = BIC(lGAMp,lGAMp1,lGAMp2,lGAMp3,lGAMp4,lGAMp5,lGAMp6,lGAMp7,lGAMp8,lGAMp9,lGAMp10,lGAMp11,lGAMp12,lGAMp13,lGAMp14,lGAMp15,lGAMp16,lGAMp17,lGAMp18,
             lGAMp19,lGAMp21,lGAMp22,lGAMp23,lGAMp24,lGAMp25,lGAMp26,lGAMp27,lGAMp28,lGAMp29,lGAMp30,lGAMp31,lGAMp32,lGAMp33,lGAMp34,lGAMp35,lGAMp36,lGAMp37, 
              lGAMp38,lGAMp39,lGAMp40,lGAMp41,lGAMp42,lGAMp43,lGAMp44,lGAMp45,lGAMp46,lGAMp47)
which.min(BICall8[,2])
#============================================================
#Chosen base model for GAM
#============================================================
lGAMp38<- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr"),data=dataset,family="poisson", method = "REML", select = T)
anova(lGAMp47,lGAMp38,test = "Chisq")
cor(dataset$Restaurant,predict(lGAMp38, type = "response"))
summary(lGAMp38)
#============================================================
#Base GAM model with spatial component
#============================================================
lGAMspat_coor <- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+ s(x,y, bs = "tp"),data=dataset,family="poisson", method = "REML", select = T)
# lGAMspat_neigh <- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+ fitted(GLM.3),data=dataset,family="poisson", method = "REML", select = T)
# lGAMspat_neigh_coor <- gam(Restaurant~s(I(Tourst_idx^0.5))+ s(Tax_index, bs = "cr")+s(Hous1or2, bs = "cr")+ s(Colle_deg)+s(I(Womn_20_24^0.5))+s(Employed, bs = "cr")+ s(x,y, bs = "tp")+ fitted(GLM.3),data=dataset,family="poisson", method = "REML", select = T)
# BIC(lGAMspat_coor, lGAMspat_neigh, lGAMspat_neigh_coor)
#============================================================
#Comparison between GAM base and GAM with spatial component
#============================================================
AIC(lGAMspat_coor, lGAMp38)
cor(dataset$Restaurant,predict(lGAMspat_coor, type = "response"))
summary(lGAMspat_coor)

# dataset_red = dataset[,c(3:26)]
# gams = apply(dataset_red, 2, function (q) gam(Restaurant ~ s(q), data = dataset_red, family= "poisson", method = "REML"))
# for (i in c(1:24)){
# gamsBIC = BIC(gams[[i]]) 
# }
# which.min(gamsBIC)

summary(lGAMspat_coor)
gam.check(lGAMspat_coor)
concurvity(lGAMspat_coor) #https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/concurvity.html
# sampleCV =  dataset[sample(nrow(dataset), 200),]
shapiro.test(residuals(lGAMspat_coor,type='deviance'))
ad.test(residuals(lGAMspat_coor,type='deviance'))#residuals don't follow normality distribution

plot(lGAMspat_coor)
plot(lGAMp38)
ad.test(residuals(lGAMp38,type='deviance'))#residuals don't follow normality distribution
# s(I(Womn_20_24^0.5)) is almost penalized out
anova(lGAMp38,lGAMspat_coor, test = "Chisq") # to analyse the variance of the model
anova(lGAMp38) 
anova(lGAMspat_coor,lGAMp38, test = "Chisq") # to analyse components of the model

exp(coef(lGAMspat_coor))  # Exponentiated coefficients
#============================================================
# Cutt-offs
#============================================================
dataset$dev_GAM_spat <- with(dataset, residuals(lGAMspat_coor,type='deviance'))
dataset$rest_GAM_spat <- with(dataset, predict(lGAMspat_coor,type='response'))
#============================================================
# defining number of clusters
#============================================================
ggdendrogram(diana(dataset$dev_GAM_spat, metric = "euclidian")) + coord_flip() 
# showData(dataset, placement='-20+200', font=getRcmdr('logFont'), maxwidth=80, maxheight=30)
# with(dataset, Hist(dev_GAM_spat, scale="percent", breaks="Sturges", col="darkgray"))
# clara(dataset$dev_GAM_spat,5)
# clara(dataset$dev_GAM_spat,5)$clustering
#============================================================
# Including number of clusters in clustering method
#============================================================
dataset$classes3=clara(dataset$dev_GAM_spat,3)$clustering

scatterplot(classes3~dev_GAM_spat, reg.line=FALSE, smooth=FALSE, spread=FALSE, boxplots=FALSE, span=0.5, ellipse=FALSE, levels=c(.5, .9), data=dataset, xlab = "Deviances", ylab = "# of clusters")
scatterplot(Restaurant~rest_GAM_spat,reg.line=lm, data=dataset, spread = F,xlab = "Restaurant prediction - GAM with spatial component")
# 
# breaks=seq(1, 5, 1)
# labels=c("Low - over", "Low - under", "No change", "High - under", "High - over")
# np <- findInterval(l_m$classes3, breaks)
# colors <- c("red","pink", "white", "skyblue2", "blue")
# plot(l_m, col = colors)#colors[np] manually set the color for each region
# mtext("Map of restaurant potential", cex=1.5, side = 3, line = 1)
# legend("bottom", legend = labels, ncol=5, fill = colors, bty = "n",cex=0.7)

#============================================================
# The Empirical RUle
#============================================================
d = dataset$dev_GAM_spat

Ic <- d[which(d >= mean(d)-sd(d) & d <= mean(d)+sd(d))]
Ica <- d[which(d >= mean(d)-sd(d) & d <= sd(d))]
II_intv_left <- d[which(d >= mean(d)-2*sd(d) & d < mean(d)-sd(d))]
II_intv_right <- d[which(d <= mean(d)+2*sd(d) & d > mean(d)+ sd(d))]                                           
II_intv_righta <- d[which(d <= 2*sd(d) & d >  sd(d))]                                           

# II_dupl_left <- c(Ic,II_intv_left, II_intv_right)
# II_dupl_right <- c(Ic,II_intv_right, II_intv_left)
# x <- sapply(II_intv_left, function(x) sum(as.numeric(II_intv_left == x)))
# II_left <- II_dupl_left[x==1]
# y <- sapply(II_dupl_right, function(y) sum(as.numeric(II_dupl_right == y)))
# II_right <- II_dupl_right[y==1]

III_intv_left <- dataset$dev_GAM_spat[which(dataset$dev_GAM_spat < mean(d)-2*sd(d))] 
III_intv_right <- dataset$dev_GAM_spat[which(dataset$dev_GAM_spat > mean(d)+2*sd(d))]
III_intv_righta <- dataset$dev_GAM_spat[which(dataset$dev_GAM_spat > 2*sd(d))]

# III_dupl_left <- c(d, II_intv_left)
# III_dupl_right <- c(d, II_intv_right)
# q <- sapply(III_dupl_left, function(q) sum(as.numeric(III_dupl_left == q)))
# III_left <- III_dupl_left[q==1]
# w <- sapply(III_dupl_right, function(w) sum(as.numeric(III_dupl_right == w)))
# III_right <- III_dupl_right[w==1]

a= as.data.frame((d%in%Ica)*1)
b = as.data.frame((d%in%II_intv_left)*2)
c = as.data.frame((d%in%II_intv_righta)*4)
t = as.data.frame((d%in%III_intv_left)*3)
e = as.data.frame((d%in%III_intv_righta)*5)
z = as.data.frame(c(a,b))
colnames(z) = c("I","II")
z[z==0] <- NA
z <- within(z, I <- ifelse(is.na(I), II, I))
z[2]<-c
z[z==0] <- NA
z <- within(z, I <- ifelse(is.na(I), II, I))
z[2]<-t
z[z==0] <- NA
z <- within(z, I <- ifelse(is.na(I), II, I))
z[2]<-e
z <- as.data.frame(within(z, I <- ifelse(is.na(I), II, I)))
dataset$category = as.numeric(unlist(z[1]))
#============================================================
# The map of estimation
#============================================================
par(mfrow=c(1,1))
breaks=seq(1, 5, 1)
labels=c("No change", "Underestimated", "Extremely underestim", "Overestimated", "Extremely overestim" )
np <- findInterval(dataset$category, breaks)
colors <- c("white","orange", "red","skyblue", "blue")
plot(l_m, col = colors[np])#colors[np] manually set the color for each region
mtext("Map of Lisbon - Estimation", cex=1.1, side = 3, line = 1)
legend("bottomright", legend = labels, fill = colors, bty = "n",cex=0.7, inset = c(-0.1,0))
# map.scale(-98000,-96000,len=12,"Kilometers",2,0.5,sfcol='red')
north.arrow(xb = -95002.74, yb=-96000, len=0.05,col="cyan", lab = "???", cex.lab = 2)
# map.scale(-90000,-106000,len = 5,units = "Kilometers",ndivs = 2,subdiv = 4)
h = hist(d, col = "grey", main = "Histogram of residuals")
xfit = seq(min(d), max(d), length = 30)
yfit = dnorm(xfit, mean = mean(d), sd = sd(d))
yfit = yfit*diff(h$mids[1:2])*length(d)
lines(xfit,yfit, lw = 2, col = "red")

density = density(d)
plot(density, main = "Kernel density for deviation and cut offs", xlab = "Deviation", col = "red") +
polygon(density, col = "grey") + abline(v = c(sd(d), mean(d)-sd(d), 2*sd(d), mean(d)-2*sd(d)), col = "red", lty = "dashed")  + text(0.8,0.3,labels = "mean + sd = 0.86", cex = 0.8) + text(2.4,0.4,labels = "mean + 2 x sd = 2.02", cex = 0.8)

#============================================================
# Chebyshev???s Theorem
#============================================================
setwd("F:/Master thesis/Data/New_Lisboa") 
dataset$rest_GAM_spat = as.numeric(dataset$rest_GAM_spat)
l_m@data = as.data.frame(dataset)
writePolyShape(l_m,fn = "Lisboa_with_estimation_shift")
