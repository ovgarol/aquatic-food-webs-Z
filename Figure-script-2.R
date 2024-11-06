matching.colors = c(alpha('tomato',0.15),alpha('royalblue',0.15),'gold2')
site.palette = c('tan1','brown','#708238')

#db.mm = read.csv("minimal_model_meca.csv",comment.char = '#') # to load alternative models
db.mm = db.sp # using the mechanistic food web model
#db.mm$tag = paste(db.mm$group,db.mm$type,sep=' - ')
db.mm$mean.D = 0.5*(log(db.mm$min.esd)+log(db.mm$max.esd)) # calculating median log size
for(i in 1:length(db.mm$tag)) db.mm$mean.OPS[i] = Min.Spec(exp(db.mm$mean.D[i]),db.mm$tag[i],db.mm)###TEMPORAL

my.pal = colorRampPalette(hcl.colors(5,palette='Zi'))(5)

par(bty='o',mfrow=c(3,6),mai=0.051+c(.2,.2,.2,.2),oma=c(5,5,1,0))

n = 1
error.limit = 2/sqrt(3) # width of feeding kernel

## metric values across food webs
n.species = numeric() # number of species
n.links = numeric() # number of links
n.obs = numeric() # number of observations
connectance = numeric() # n.links/n.obs^2
representativeness=numeric() # accuracy
n.mm = numeric() # number of feeding guilds included in the mechanistic model
eco.type = numeric() # ecosystem type
db.touched = data.frame() # list of including feeding guilds

#real.net = list() # list of real observed foodwebs
art.net = list() # list of artificial reduced foodwebs

############
## Figure 4: reconstruction of real food webs
############
what.to.plot = c(11,5,3,19,2,18,4,6,22,13,1,16,14,15,20,17,21,12) # sorted by ecosystem type
#what.to.plot =  c(20,17,12) # TO DELETE
coor = NULL

## reading and preparing GATEWAy database (https://doi.org/10.25829/iDiv.283-3-756).
## download from https://idata.idiv.de/ddm/Data/ShowData/283?version=3

db = read.csv("./external_data/283_2_FoodWebDataBase_2018_12_10.csv")
db = db[!is.na(db$autoID),] # removing empty lines
db = db[db$interaction.type%in%c('predacious','bactereeivorous'),] # selecting predacius interactions
db = db[db$ecosystem.type %in% c('marine','lakes','streams'),] # selecting aquatic ecosystems
db$type = paste(db$con.movement.type,db$con.metabolic.type,db$interaction.type)
db$cols = as.numeric.factor(db$con.metabolic.type)

esds = 2*(db$con.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
opss = 2*(db$res.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6

## Checking data
db0 = db # keeping a copy of filtered GATEWAy
length(unique(db0$foodweb.name[db$study.site%in%unique(db$study.site)[what.to.plot]]))
unique(db0$study.site[db$study.site%in%unique(db$study.site)[what.to.plot]])
length(db0$study.site[db$study.site%in%unique(db$study.site)[what.to.plot]])

##
total.obs  = 0 # observations counter
for(idx in what.to.plot){
  db = db0
  db = db[db$study.site == unique(db$study.site)[idx],] 
  db = db[db$con.movement.type=='swimming',] # restricting to swimming predators
  db = db[db$con.mass.mean.g.!=-999 & db$res.mass.mean.g.!=-999, ] # removing NA values
  db = db[!is.na(db$autoID),] # removing NA values
  db= db[!(db$res.taxonomy=='Calanus finnmarchicus'|db$res.mass.mean.g == 2.099800e-04),]
  
  ## continue if the site has no data
  if(length(db$autoID)==0){
    print(paste(idx, 'failed!!!'))
    next
  } 
  
  ## extracting data to plot the map
  eco.type[n] = unique(db$ecosystem.type)[1]
  coor[idx] = paste(db$study.site[1],db$ecosystem.type[1],db$latitude[1],db$longitude[1],sep=',') 
  
  ## size calculations from mass
  esds = 2*(db$con.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  opss = 2*(db$res.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  
  ## checking if observations and feeding guilds matching
  represented = as.numeric(esds==-1) # Boolean vector 1:observation is included in model 0:observation is not included in model 
  touched = as.numeric(db.mm$min.esd == -1)  # Boolean vector 1:feeding guild is used in model 0:feeding guild is not used in model 
  for(i in 1:length(db.mm$tag))for(j in 1:length(esds)){
    if(!represented[j])if(esds[j]>db.mm$min.esd[i] & esds[j]<db.mm$max.esd[i])if(abs(Min.Spec(esds[j],db.mm$tag[i],db.mm)-log(opss[j]))<error.limit){
      represented[j] = 1
      touched[i] = 1
    }  
  }
  
  ##plotting size-based model
  if(T){
    x.lim = range(esds,opss)#*c(0.8,1.2)
    y.lim = range(opss,esds)#*c(0.9,1.1)
    if(n==12){
      x.lim=range(esds,opss)*c(0.75,2.)
      y.lim =range(esds,opss)*c(0.75,2.)
    }
    
    plot(esds,opss,type='n',
         col=alpha('lightgray',0.95),
         bg=alpha(hcl.colors(10)[as.numeric.factor(db$foodweb.name)],0.95),
         #col=cols,
         lwd=0.1,
         pch=19,cex=.75,
         log='xy',
         xlim=x.lim,ylim=y.lim,
         xaxt='n',yaxt='n',
         xlab='',
         ylab=''
    )
    
    y=c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)
    y.lab = -7:1
    
    axis(side=1, at=(y), labels=y.lab,las=1,cex.axis=1.5)
    axis(side=2, at=(y), labels=y.lab,las=1,cex.axis=1.5)
    text(x.lim[1],y.lim[2],letters[n],adj=c(0,1),font=2,cex=1.5)
    eco.col=site.palette[1]
    if(eco.type[n]=='lakes') eco.col=site.palette[2]
    if(eco.type[n]=='streams') eco.col=site.palette[3]
    mtext(side=3,unique(db$study.site),outer=F,line=0,adj=0,font=2,col=eco.col,cex=1.25)
    
    abline(a=0,b=1,lty=2)
    
    ## to plot the allometric rule
    if(F){
      x = log(10^(-60:100 * 0.1)*1e6 )
      y = exp(-1.65+x-0.011*x**2)
      x = exp(x)
      lines(x,y)
      polygon(c(x,rev(x)),c(1/(2*error.limit)*y,rev(2*error.limit*y)),border='red',col=alpha('pink',0.2),lty=3) 
    }
    
    points(esds,opss,col=matching.colors[1+as.numeric(represented)],pch=19)
    if(F)for(i in 1:length(db.mm$tag)) if(touched[i]) minimal.model(db.mm$min.esd[i],db.mm$max.esd[i],db.mm$s[i],col='black',lty=1,lwd=2) # overlays the mechanistic model
    
  }
  
  n.mm[n] = length(db.mm[touched,]$tag)
  n.species[n] = length(unique(c(db$con.taxonomy,db$res.taxonomy)))
  n.links[n] = length(unique(paste(db$con.taxonomy,db$res.taxonomy)))
  n.obs[n] = length(db$autoID)
  connectance[n] = n.links[n] /n.species[n] **2
  representativeness[n]=100*sum(represented)/length(esds)
  total.obs = total.obs + length(db$autoID)
  
  ##printing data for table 2
  if(T){
    if(n==1)print(paste(c('Food-Web','#species','#links','#observations','#guilds','represented links'),collapse = ' & '))
    resul = paste(round(c(n.species[n],n.links[n],n.obs[n],n.mm[n],representativeness[n]),digits = 2),sep='&',collapse=' & ')
    print(paste(paste(unique(db$study.site),length(unique(db$foodweb.name))),resul,sep=' & '))
  }
  
  ## analysis of food web complexity
  #for(i in 1:length(db.mm$tag)) if(touched[i]) minimal.model(db.mm$min.esd[i],db.mm$max.esd[i],db.mm$s[i],col='black',lwd=2.5)
  db.touched = rbind(db.touched,touched)
  
  ## including text in each subplot
  if(T){
    if(n==1){
      text(x.lim[1],exp(log(y.lim[2])-0.15*diff(log(y.lim))),paste(round(representativeness[n],digits=1),'accuracy (%)'),adj=c(0,1),font=1,cex=1.5)
      text(x.lim[1],exp(log(y.lim[2])-0.3*diff(log(y.lim))),paste(n.mm[n], 'guilds'),adj=c(0,1),font=2,cex=1.5)
      
    }
    else{
      text(x.lim[1],exp(log(y.lim[2])-0.15*diff(log(y.lim))),paste(round(representativeness[n],digits=1)),adj=c(0,1),font=1,cex=1.5)
      text(x.lim[1],exp(log(y.lim[2])-0.3*diff(log(y.lim))),n.mm[n],adj=c(0,1),font=2,cex=1.5)
    } 
  }

  n=n+1
}

db = db0
mtext(side=1,TeX(paste0('log$_{10}$ ','predator size (m)')),outer=T,line=3,adj=0.5,cex=2)
mtext(side=2,TeX(paste0('log$_{10}$ ','prey size (m)')),outer=T,line=2,adj=0.5,las=0,cex=2)

## printing metric values across food webs
mean(representativeness)
sd(representativeness)

mean(n.species)
sd(n.species)

mean(n.links)
sd(n.links)

mean(n.obs)
sd(n.obs)

mean(connectance) 
sd(connectance) 

mean(n.mm) 
sd(n.mm) 
