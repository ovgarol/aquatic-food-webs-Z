library(oce)
library(rnaturalearth)
library(rnaturalearthdata)
library(ocedata)

data("coastlineWorld")
data("coastlineWorldFine")

sites = read.csv("./data/map_coordinates.csv")
sites = as.data.frame(sites)
eco.col=rep('tan1',length(sites$site))
eco.col[sites$type=='lakes'] = 'brown'  
eco.col[sites$type=='streams'] = 'slategray'

sites$col = eco.col

####
## Figure A2
####

pdf('./figures/MAPS.pdf',family='Helvetica',10,4)
par(mfrow=c(1,2),family='Helvetica',bty='n')
layout(matrix(c(1,1,0,2), 2, 2, byrow = F),widths = c(2,1),heights = c(1,3))

mapPlot(coastlineWorld,
        projection="+proj=robin",
        border=NA,
        grid=T,
        clip=T,
        col="lightgray")
in.small.map = sites$ID%in%c('j','m','n','i','j','e','d')
mapPoints(sites$longitude[!in.small.map],sites$latitude[!in.small.map],cex=3,pch=21,col='whitesmoke',lwd=0.5,bg=sites$col[!in.small.map])  
mapText(sites$longitude[!in.small.map],sites$latitude[!in.small.map],sites$ID[!in.small.map],font=2,col='white')
title(main='Ecosystem locations',adj=0,line=0)

par(bty='o')
mapPlot(coastlineWorldFine,
        #projection="+proj=moll",
        border=NA,
        grid=T,
        clip=T,
        col="lightgray",
        longitudelim = c(-10, 5), latitudelim =  c(60, 50))
mapPoints(sites$longitude[in.small.map],sites$latitude[in.small.map],cex=3,pch=21,col='whitesmoke',lwd=0.5,bg=sites$col[in.small.map])  
mapText(sites$longitude[in.small.map],sites$latitude[in.small.map],sites$ID[in.small.map],font=2,col='white')

dev.off()

