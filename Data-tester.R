setwd('~/Documents/articles/functionalFeedingMode/Z/data')
database <- read.csv("./database.csv")
specialization <- read.csv("./database_specialization.csv")
specialization$alpha = specialization$a
ecosystem <- read.csv("./database_ecosystems.csv")
ecosystem_so <- read.csv("./database_ecosystems_size_model.csv")
coordinates <- read.csv("./map_coordinates.csv")
mc <- read.csv("./montecarlo-sampling-2885.csv")
database_calculated <- read.csv("./database_calculated.csv")
dinos <- read.csv("./Garcia-Oliva2022a.csv", comment.char="!")

#Fig2 uses database.csv
par(mfrow=c(2,3),pch=19)
for(i in unique(database$group)){
  temp = database[database$group==i,]
  plot(temp$esd,temp$opt,log='xy',main=i)
} 

#Fig3 uses database.csv
par(mfrow=c(2,2),pch=19)
plot(database$esd,database$opt,log='xy')

#Fig4 uses database_specialization.csv
par(mfrow=c(2,2),pch=19)
plot(specialization$s*specialization$stiffness,specialization$alpha*(specialization$alpha>0),log='') # Fig4a
plot(specialization$mean.D,specialization$stiffness,log='xy') # Fig4b
plot(specialization$mean.D,specialization$mean.m,log='x') # Fig4c
plot(specialization$esd,specialization$s*specialization$stiffness,log='x') # Fig4d

#Fig5 uses database_ecosystems.csv
par(mfrow=c(3,6))
for(i in unique(ecosystem$site)){
  db = ecosystem[ecosystem$site==i,]
  plot(db$esd,db$ops,col=c('tomato','dodgerblue')[1+db$represented],log='xy',main=i)
}
  
#Fig6 uses montecarlo-sampling-2885.csv
par(mfrow=c(1,1))
plot(mc$n.obs,mc$total.acc/100,log='x',ylim=c(0,1))

#EDF1 uses database_specialization.csv
par(mfrow=c(1,5))
for(i in unique(specialization$group)){
  temp = specialization[specialization$group==i,]
  plot(temp$s,temp$alpha*(temp$alpha>0),main=i,ylim=c(0,1.3),xlim=c(-5,5))
} 

#EDF2 uses database_calculated.csv
par(mfrow=c(1,2))
plot(database_calculated$opt,database_calculated$ops.2,log='yx')
plot(database_calculated$opt,database_calculated$ops.0,log='yx')

#EDF3 uses map_coordinates.csv
par(mfrow=c(1,1))
plot(coordinates$longitude,coordinates$latitude,type='p',xlim=c(-180,180),ylim=c(-90,90))
text(coordinates$longitude,coordinates$latitude,coordinates$site)

#EDF4 uses database_ecosystems_size_model.csv
par(mfrow=c(3,6))
for(i in unique(ecosystem_so$site)){
  db = ecosystem_so[ecosystem_so$site==i,]
  plot(db$esd,db$ops,col=c('tomato','dodgerblue')[1+db$represented],log='xy',main=i)
}

#EDF5 uses Garcia-Oliva2022a.csv
par(mfrow=c(1,3))
for(i in c(0,-1,1)){
  dinos.s = dinos[dinos$Specialization.factor..s==i,]
  plot(dinos.s$Equivalent.spherical.diameter..ESD..um.,dinos.s$Optimal.prey.size...OPS..um.,log='xy',ylim=c(3,100),xlim=c(3,100))
  
}

