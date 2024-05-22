library(oce)
library(rnaturalearth)
library(rnaturalearthdata)
library(ocedata)

data("coastlineWorld")
data("coastlineWorldFine")

sites <- read.csv("~/Documents/articles/functionalFeedingMode/map_coordinates.csv")
sites = as.data.frame(sites)
eco.col=rep('tan1',length(sites$site))
eco.col[sites$type=='lakes'] = 'brown'  
eco.col[sites$type=='streams'] = 'tan4'

sites$ID = letters[1:18]
sites$col = eco.col

####
## Figure A2
####

par(mfrow=c(1,2))

mapPlot(coastlineWorld,
        projection="+proj=robin",
        border=NA,
        grid=T,
        clip=T,
        col="lightgray")
mapPoints(sites$longitude,sites$latitude,cex=3,pch=21,col='whitesmoke',lwd=0.5,bg=sites$col)  
mapText(sites$longitude,sites$latitude,sites$ID,font=2,col='white')

mapPlot(coastlineWorldFine,
        #projection="+proj=moll",
        border=NA,
        grid=T,
        clip=T,
        col="lightgray",
        longitudelim = c(-10, 5), latitudelim =  c(60, 50))
mapPoints(sites$longitude,sites$latitude,cex=3,pch=21,col='whitesmoke',lwd=0.5,bg=sites$col)  
mapText(sites$longitude,sites$latitude,sites$ID,font=2,col='white')

