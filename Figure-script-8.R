ww = 20
#set.seed(3985414->seed)
set.seed(abs(as.integer(1000*rnorm(1)))->seed)

calc.structure = function(db,ww){
  ##removing some species
  db.all = db
  if(length(db$autoID)<ww){return(NA)}
  db = db[db$autoID%in%sample(unique(db$autoID),ww),]
  if(ww == 0) db = db.all # if ww == 0 calculate the entire food web
  
  ## size calculations from mass
  esds = 2*(db$con.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  opss = 2*(db$res.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  
  ## checking if observations and feeding guilds matching
  represented = as.numeric(esds==-10) # Boolean vector 1:observation is included in model 0:observation is not included in model 
  touched = as.numeric(db.mm$min.esd == -1)  # Boolean vector 1:feeding guild is used in model 0:feeding guild is not used in model 
  for(i in 1:length(db.mm$tag))for(j in 1:length(esds)){
    if(!represented[j])if(esds[j]>db.mm$min.esd[i] & esds[j]<db.mm$max.esd[i])if(abs(Min.Spec(esds[j],db.mm$tag[i],db.mm)-log(opss[j]))<error.limit){
      represented[j] = 1
      touched[i] = 1
    }  
  }
  
  accuracy=100*sum(represented)/length(esds)
  
  
  ## calculating accuracy over the entire data set 
  ## size calculations from mass to entire data set
  esds = 2*(db.all$con.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  opss = 2*(db.all$res.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  represented = as.numeric(esds==-10) # Boolean vector 1:observation is included in model 0:observation is not included in model 
  
  for(i in 1:length(db.mm$tag))for(j in 1:length(esds)){
    if(touched[i])if(!represented[j])if(esds[j]>db.mm$min.esd[i] & esds[j]<db.mm$max.esd[i])if(abs(Min.Spec(esds[j],db.mm$tag[i],db.mm)-log(opss[j]))<error.limit){
      represented[j] = 1
    }  
  }
  
  total.accuracy=100*sum(represented)/length(esds)
  
  ## return to main loop
  return(data.frame(name=db$study.site[1],
                    type=unique(db$ecosystem.type)[1],
                    is.complete= (ww==0),
                    part.acc=accuracy, 
                    total.acc=total.accuracy, 
                    fw=as.character(paste(touched,collapse='')), 
                    n.species.pred=length(unique(db$con.taxonomy)), 
                    n.species.prey=length(unique(db$res.taxonomy)),
                    n.species.total=length(unique(c(db$res.taxonomy,db$con.taxonomy))),
                    n.obs=length(db$type)))
  
}


lower.boundary = function(x,y) {
  #result = aggregate(y ~ x, FUN = min)
  result = aggregate(y ~ x, FUN = function(x) quantile(x,0.05))
  result = approx(result$x,result$y,exp(seq(0,8,by=0.1)))
  return(result)
}

upper.boundary = function(x,y) {
  #result = aggregate(y ~ x, FUN = max)
  result = aggregate(y ~ x, FUN = function(x) quantile(x,0.95))
  result = approx(result$x,result$y,exp(seq(0,8,by=0.1)))
  return(result)
}


jaccard = function(A,B){
  A = as.character(unlist(strsplit(A,"")))
  B = as.character(unlist(strsplit(B,"")))
  a_11 = sum(A == 1 & B == 1) 
  b_01 = sum(A == 0 & B == 1)
  c_10 = sum(A == 1 & B == 0)
  similarity = a_11/(a_11 + b_01 + c_10)
  return(similarity)
}

import.files  = function(){
  file_list <- list.files(pattern = 'montecarlo-sampling-.*\\.csv')
  data_frames <- list()
  col_types <- c(fw = "character", default = NULL)
  
  for (file in file_list) {
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE, colClasses = col_types)
    data_frames[[file]] <- data
  }
  combined_df <- do.call(rbind, data_frames)
  #print(combined_df)
  row.names(combined_df) = 1:length(combined_df$name)
  # write.csv(combined_df, "combined_data.csv", row.names = FALSE)
  return(combined_df)  
}

###############
# to get new random sampling every time it is run 
list.samples = NULL

n = 1
n.samples = 20
for(idx in what.to.plot){
  ## data
  db = db0
  db = db[db$study.site == unique(db$study.site)[idx],] 
  db = db[db$con.movement.type=='swimming',] # restricting to swimming predators
  db = db[db$con.mass.mean.g.!=-999 & db$res.mass.mean.g.!=-999, ] # removing NA values
  db = db[!is.na(db$autoID),] # removing NA values
  db = db[!(db$res.taxonomy=='Calanus finnmarchicus'|db$res.mass.mean.g == 2.099800e-04),]
  
  ## continue if the site has no data
  if(length(db$autoID)==0){
    print(paste(idx, 'failed!!!'))
    next
  } 
  
  ## reference run
  list.samples= rbind(list.samples,calc.structure(db,0))   
  
  ## calculating the censored foodweb
  if(length(db$autoID)>n.samples)possible.w = sample(1:length(db$autoID),n.samples) ##100 samples
  else possible.w = 1:length(db$autoID) ##100 samples
  
  for(ww in possible.w)for(i in 1:2){
    list.samples = rbind(list.samples,calc.structure(db,ww)) 
    n = n+1
  }
  print(idx)
}

list.samples = na.omit(list.samples)
write.csv(list.samples, paste0("montecarlo-sampling-",seed,".csv"), row.names=FALSE)

#################

db = import.files()
db = db[!duplicated(db),]
db = na.omit(db)
#db = db[!db$type%in%c('lakes','streams'),]
n.db = length(db$name)

fw.l = unique(db$fw)[order(as.numeric(unique(db$fw)))]
fw.n = length(fw.l)
sim.matrix = matrix(nrow=fw.n,ncol=fw.n )
for(i in 1:fw.n)for(j in 1:fw.n)sim.matrix[i,j]=jaccard(fw.l[i],fw.l[j])

####
## Heatmap cluster analysis of foodwebs (not usede)
####
sim.matrix=data.frame(sim.matrix)
row.names(sim.matrix) = fw.l
colnames(sim.matrix) = fw.l
heatmap(as.matrix(sim.matrix))

db$fw.ref = NA
db$fw.sim = NA
for(i in 1:length(db$name)){
  fw.ref = db$fw[(db$name[i]==db$name)&db$is.complete]
  db$fw.ref[i] = fw.ref
  db$fw.sim[i] = jaccard(fw.ref,db$fw[i])
}

db$score = 0.5*(db$total.acc/100+db$fw.sim)
db$score0 = pmin(db$total.acc/100,db$fw.sim)

palette(alpha(c('brown','tan1','tan4'),0.15))
#unique(as.factor(db$type))

#############
# Figure 6
############

par(bty='o',mfrow=c(1,1),mai=c(.8,.8,.2,.2),oma=c(0,0,0,0))
#plot(total.acc~part.acc,data=db)
plot(fw.sim~n.obs,data=db,log='x',
     xaxt='n',#yaxt='n',
     col=alpha('lightgray',0.025),#as.factor(db$type),
     pch=19,xlim=c(1,3000),
     ylab='Jaccard index',
     xlab='number of observations')
title(main='food-web similarity',adj=0)
lines(approx(db$n.obs,db$fw.sim,exp(seq(0,8,by=0.5))),lwd=3)
abline(v=200,col='red')

smallest.values = lower.boundary(db$n.obs,db$fw.sim)
lines(supsmu(smallest.values$x,smallest.values$y),type='l')
largest.values = upper.boundary(db$n.obs,db$fw.sim)
lines(supsmu(largest.values$x,largest.values$y),type='l')

y=c(1e0,1e1,1e2,1e3)
axis(side=1, at=(y), las=1,cex.lab=1.)
y=c(2:9*1e0,2:9*1e1,2:9*1e2,2:9*1e3)
axis(side=1, at=(y), labels = NA, las=1,cex.lab=1.,cex=0.5,tck=-0.025)

