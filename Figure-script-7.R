plot.it=T
iter.limit = 30 #minimum iteration for calculating the Trophic Level

#######
error.limit = 2/sqrt(3) # width of feeding kernel

calc.TP = function(net,is.exact=T){
  a = as.matrix(as_adjacency_matrix(net,sparse = F)) 
  a.list = list()
  n.mm = length(V(net))
  diagonal = diag(x = 1, n.mm, n.mm, names = F)

  ss = V(net)$x
  
  ##random TL
  n = 1
  if(!is.exact)for(k in 1:1000){
    a2=a
    for(i in 1:n.mm)for(j in 1:n.mm){
      if(a[j,i]>0) a2[j,i] = runif(1) >0.5
      #if(a[j,i]>0) a2[j,i] = runif(1,0.0,p[j,i]) #>0.5
      #if(a[j,i]>0) a2[j,i] = runif(1,0.01,0.01+0.99*exp(-1.5*(ss[i]-mean(ss[a[j,]>0]))**2)) >0.5
    }
    
    in.degree = degree(net,mode='in')
    in.degree[in.degree==0]=1
    
    a2 = scale(a2, center=FALSE, scale=in.degree+0e-10)

    if(det(diagonal-t(a2)) == 0) {print(k); next}
    
    TL = solve(diagonal-t(a2),rep(1,n.mm))
    a.list[[n]]=TL  
    n=n+1
    if(n>=iter.limit) break
  }
  
  TL = apply(simplify2array(na.omit(a.list)), 1, mean)
  TL.sd = apply(simplify2array(na.omit(a.list)), 1, sd)
  #TL.min = apply(simplify2array(a.list), 1, min)-min(TL)+1
  #TL.max = apply(simplify2array(a.list), 1, max)-min(TL)+1
  
  ##deterministic TL
  if(is.exact){
    in.degree = degree(net,mode='in')
    in.degree[in.degree==0]=1
    #a = scale(a, center=FALSE, scale=colSums(a)+1e-10)
    a = scale(a, center=FALSE, scale=in.degree+1e-10)
    TL = solve(diagonal-t(a),rep(1,n.mm))
    TL.sd = NA
  }
  
  return(data.frame(TL=TL+1,sd=TL.sd))
}

construct.fw = function(db){
  esds = 2*(db$con.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  opss = 2*(db$res.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  
  meta.data <- data.frame(
    size = log(c(opss,esds)),
    name = c(db$res.taxonomy,db$con.taxonomy)
  )
  
  relations = data.frame(
    from = db$res.taxonomy,
    to = db$con.taxonomy
  )
  net = graph_from_data_frame(relations,vertices=meta.data$names, directed = T)
  net = simplify(net, remove.multiple = T, remove.loops = T)
  V(net)$x = 0 * 1:length(V(net))
  
  meta.data=aggregate(meta.data,by=list(meta.data$name), mean)
  for(i in 1:length(V(net))) V(net)$x[i] = meta.data$size[which(V(net)$name[i]==meta.data$Group.1)][1] 
  
  TL = calc.TP(net,is.exact=F)
  
  posi.2 = TL$TL
  posi.1 = V(net)$x
  l = matrix(c(posi.1,posi.2),ncol=2)
  V(net)$y = posi.2
  
  if(plot.it){
    plot(net,
         vertex.color='dodgerblue',#hcl.colors(10)[posi.1/max(posi.1)*9+1],
         vertex.label=NA,#vertex.size=8*(1+as.numeric(n.mm<30)),
         edge.color=alpha('lightgray',0.5),
         edge.arrow.size=.25,edge.color='lightgray',
         add=F)#,layout=l)
    text(0,-0.5,paste('N =',length(V(net))),adj=c(0,0),font=2)
    text(0,-0.75,paste('C =',round(length(E(net))/length(V(net))**2,digits=2)),adj=c(0,0),font=2)
  }
  
  return(data.frame(name=V(net)$name, lesd=V(net)$x, TP=V(net)$y,sd=TL$sd))
}

fit.fw.model = function(db,w){

  if(F){
    ## creating data frame for mechanistically defined model
    db.sp = data.frame(min.esd=0,max.esd=0,s=0,tag=NA)
    ii = 1 # auxiliary iterator
    #if(is.na(w)) w = 3.5 # universal width for each feeding guild
    
    for(i in 1:ceiling(14/w)){ # four size classes
      for(j in -2:1){ # possible values of specialization 
        D.min = exp(ii+abs(j)+0*(j < -1)) # min predator size + offset to specialists
        D.max = exp(ii+w+abs(j)+0*(j < -1)) # max predator size + offset to specialists
        the.D = 0.5*(D.min+D.max) # mean predator size
        the.s = j*0.014*(log(the.D/1e0))**2 # Table 1, Eq(d) 
        #the.s = j*(0.15*(ii+w+1*(j < -1)+abs(j))) # alternative Eq
        the.s = j*(0.17*(log(the.D))) # linear Eq(d)
        
        if(j< -1 & (ii+w > log(1e4) )) next # ignoring additional small prey specialists in larger size classes
        db.sp = rbind(db.sp,c(D.min,D.max,the.s,paste0('tag',i,j)),stringsAsFactors = FALSE) # adding model to 
      }
      ii = ii+w
    }
    
    ## fixing typing in mechanistic model
    db.sp = db.mm# na.omit(db.sp)
    db.sp$min.esd = as.numeric(db.sp$min.esd)
    db.sp$max.esd = as.numeric(db.sp$max.esd)
    db.sp$s = as.numeric(db.sp$s)

    ## using db.sp
    #db.mm = db.sp # using the mechanistic food web model
    #db.mm$tag = paste(db.mm$group,db.mm$type,sep=' - ')
    db.mm$mean.D = 0.5*(log(db.mm$min.esd)+log(db.mm$max.esd)) # calculating median log size
    db.mm$mean.OPS = db.mm$mean.D-3
    for(i in 1:length(db.mm$tag)) db.mm$mean.OPS[i] = Min.Spec(db.mm$mean.D[i],db.mm$tag[i],db.mm)###TEMPORAL
    
  }
    

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
  
  ##modeled web
  if(T){
    touched = as.logical(touched)
    n.mm = length(db.mm[touched,]$tag)
    connections = matrix(1, nrow = n.mm, ncol = n.mm)
    for(i in 1:n.mm) for(j in 1:n.mm) connections[j,i] = db.mm[touched,]$mean.OPS[i]>log(db.mm[touched,]$min.esd[j]) & db.mm[touched,]$mean.OPS[i]<log(db.mm[touched,]$max.esd[j])
    links.mm = sum(connections)  
    net = graph_from_adjacency_matrix(connections, mode = "directed", diag = T)
    net = simplify(net, remove.multiple = T, remove.loops = T)
    V(net)$x = (db.mm$mean.D[touched==1])
    V(net)$name = (db.mm$tag[touched==1])
    #art.net[[n]] = net
  }
  
  
  TL = calc.TP(net,is.exact = F)
  
  posi.2 = TL$TL
  posi.1 = V(net)$x
  l = matrix(c(posi.1,posi.2),ncol=2)
  V(net)$y = posi.2
  
  
  if(plot.it){
    plot(net, 
         vertex.color='black',#,hcl.colors(10)[posi.1/max(posi.1)*9+1],
         vertex.label=NA,vertex.size=8*(1+as.numeric(n.mm<30)),
         edge.color=alpha('lightgray',0.5),
         edge.arrow.size=.25,edge.color='lightgray',add=F,layout=l)
    text(0,-0.5,paste('N =',length(V(net))),adj=c(0,0),font=2)
    text(0,-0.75,paste('C =',round(length(E(net))/length(V(net))**2,digits=2)),adj=c(0,0),font=2)
  }
  
  return(data.frame(name=V(net)$name, lesd=V(net)$x, TP=V(net)$y, sd=TL$sd))
  
}

get.broken.model = function(db,ww){
  
  if(F){
    ## creating data frame for mechanistically defined model
    db.sp = data.frame(min.esd=0,max.esd=0,s=0,tag=NA)
    ii = 1 # auxiliary iterator
    #if(is.na(w)) w = 3.5 # universal width for each feeding guild
    
    for(i in 1:ceiling(14/w)){ # four size classes
      for(j in -2:1){ # possible values of specialization 
        D.min = exp(ii+abs(j)+0*(j < -1)) # min predator size + offset to specialists
        D.max = exp(ii+w+abs(j)+0*(j < -1)) # max predator size + offset to specialists
        the.D = 0.5*(D.min+D.max) # mean predator size
        the.s = j*0.014*(log(the.D/1e0))**2 # Table 1, Eq(d) 
        #the.s = j*(0.15*(ii+w+1*(j < -1)+abs(j))) # alternative Eq
        the.s = j*(0.17*(log(the.D))) # linear Eq(d)
        
        if(j< -1 & (ii+w > log(1e4) )) next # ignoring additional small prey specialists in larger size classes
        db.sp = rbind(db.sp,c(D.min,D.max,the.s,paste0('tag',i,j)),stringsAsFactors = FALSE) # adding model to 
      }
      ii = ii+w
    }
    
    ## fixing typing in mechanistic model
    db.sp = db.mm# na.omit(db.sp)
    db.sp$min.esd = as.numeric(db.sp$min.esd)
    db.sp$max.esd = as.numeric(db.sp$max.esd)
    db.sp$s = as.numeric(db.sp$s)
    
    ## using db.sp
    #db.mm = db.sp # using the mechanistic food web model
    #db.mm$tag = paste(db.mm$group,db.mm$type,sep=' - ')
    db.mm$mean.D = 0.5*(log(db.mm$min.esd)+log(db.mm$max.esd)) # calculating median log size
    db.mm$mean.OPS = db.mm$mean.D-3
    for(i in 1:length(db.mm$tag)) db.mm$mean.OPS[i] = Min.Spec(db.mm$mean.D[i],db.mm$tag[i],db.mm)###TEMPORAL
    
  }
  
  ##removing some species
  if(length(unique(db$con.taxonomy))<ww){plot.new();return(NA)}
  db = db[!db$con.taxonomy%in%sample(unique(db$con.taxonomy),floor(ww/10*length(unique(db$con.taxonomy)))),]
   
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
  
  ##modeled web
  if(T){
    touched = as.logical(touched)
    n.mm = length(db.mm[touched,]$tag)
    connections = matrix(1, nrow = n.mm, ncol = n.mm)
    for(i in 1:n.mm) for(j in 1:n.mm) connections[j,i] = db.mm[touched,]$mean.OPS[i]>log(db.mm[touched,]$min.esd[j]) & db.mm[touched,]$mean.OPS[i]<log(db.mm[touched,]$max.esd[j])
    links.mm = sum(connections)  
    net = graph_from_adjacency_matrix(connections, mode = "directed", diag = T)
    net = simplify(net, remove.multiple = T, remove.loops = T)
    V(net)$x = (db.mm$mean.D[touched==1])
    V(net)$name = (db.mm$tag[touched==1])
    #art.net[[n]] = net
  }
  
  
  TL = calc.TP(net,is.exact = F)
  
  posi.2 = TL$TL
  posi.1 = V(net)$x
  l = matrix(c(posi.1,posi.2),ncol=2)
  V(net)$y = posi.2
  
  
  if(plot.it){
    plot(net, 
         vertex.color=hcl.colors(10)[2*ww],
         vertex.label=NA,vertex.size=10*(1+as.numeric(n.mm<30)),
         edge.color=alpha('lightgray',0.5),
         edge.arrow.size=.25,edge.color='lightgray',add=F,layout=l)
    text(0,-0.5,paste('N =',length(V(net))),adj=c(0,0),font=2)
    text(0,-0.75,paste('C =',round(length(E(net))/length(V(net))**2,digits=2)),adj=c(0,0),font=2)
  }
  
  return(data.frame(name=V(net)$name, lesd=V(net)$x, TP=V(net)$y, sd=TL$sd))
  
}


########
sum.art=NULL
sum.real=NULL

####
## Figure A6
####

par(bty='o',mfrow=c(3,6),mai=0.051+c(.2,.2,.2,.2),oma=c(5,5,1,0),las=1)

for(idx in (1:length(what.to.plot))){
  ## data
  db = db0
  db = db[db$study.site == unique(db$study.site)[what.to.plot[idx]],] 
  db = db[db$con.movement.type=='swimming',] # restricting to swimming predators
  db = db[db$con.mass.mean.g.!=-999 & db$res.mass.mean.g.!=-999, ] # removing NA values
  db = db[!is.na(db$autoID),] # removing NA values
  db = db[!(db$res.taxonomy=='Calanus finnmarchicus'|db$res.mass.mean.g == 2.099800e-04),]
  
  ## continue if the site has no data
  if(length(db$autoID)==0){
    print(paste(idx, 'failed!!!'))
    next
  } 
  
  ## extracting data to plot the map
  eco.type[n] = unique(db$ecosystem.type)[1]
  coor[idx] = paste(db$study.site[1],db$ecosystem.type[1],db$latitude[1],db$longitude[1],sep=',') 
  eco.col='tan1'
  if(eco.type[n]=='lakes') eco.col='brown'
  if(eco.type[n]=='streams') eco.col='slategray'
  
  real=construct.fw(db)
  
  ## to fix the site names
  title.site = unique(db$study.site)
  if(title.site=='NewZealandStreams') title.site = 'New Zealand Streams'
  if(title.site=='Eastern Weddell Sea Shelf') title.site = 'Eastern Weddell Sea'
  if(title.site=='Ythan Estuary, tidal Estuary of River Ythan, Forvie Nature Reserve') title.site = 'Ythan Estuary'
  if(title.site=='Puerto Rico-Virgin Islands (PRVI) shelf complex') title.site = 'PRVI shelf complex'
  if(idx ==14) title.site = paste(title.site,'(I)')
  if(idx ==15) title.site = paste(title.site,'(II)')
  mtext(side=3,title.site,outer=F,line=0,adj=0,font=2,col=eco.col,cex=1.25)

  art=fit.fw.model(db,3.4)
  #mtext(side=2,'0%',outer=F,line=0,adj=0,font=2,col='black',cex=1.25)
  
  #ss=get.broken.model(db,0.5)
  #mtext(side=2,'5%',outer=F,line=0,adj=1,font=2,col='black',cex=1.25)
  
  #ss=get.broken.model(db,1)
  #mtext(side=2,'10%',outer=F,line=0,adj=1,font=2,col='black',cex=1.25)
  
  #ss=get.broken.model(db,2)
  #mtext(side=2,'20%',outer=F,line=0,adj=1,font=2,col='black',cex=1.25)
  
  #ss=get.broken.model(db,8)
  #mtext(side=2,'80%',outer=F,line=0,adj=1,font=2,col='black',cex=1.25)
  
  esds = 2*(db$con.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  opss = 2*(db$res.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  
  if(T){
    plot(1,1,xlim=range(c(real$lesd,art$lesd)),ylim=range(c(real$TP,art$TP)),
         xaxt='n',yaxt='n',
         log='',col='gray',cex=0.5,pch=19,type='n')
    
    points(real$lesd,real$TP,bg=alpha('dodgerblue',0.25),pch=21,col='dodgerblue',cex=1)
    for(j in 1:length(real$TP)) lines(real$lesd[j]+c(0,0),real$TP[j]+real$sd[j]*c(-1,1),col=alpha('dodgerblue',0.1))
    
    for(j in 1:length(art$TP)) lines(art$lesd[j]+c(0,0),art$TP[j]+art$sd[j]*c(-1,1),col='black')
    for(j in 1:length(art$TP)) lines(art$lesd[j]+1.7*c(-1,1),art$TP[j]+c(0,0),col='black')
    
    points(art$lesd,art$TP,col='black',pch=15,cex=1.2)
    
    if(F){
      filt = art$TP>2
      if(sum(filt)>2) abline(s.real<-lm(TP~lesd,data=art[filt,],weights = 1/(1+art$sd[filt])),col='gray',lwd=3)
      else abline(s.real<-lm(TP~1,data=art[filt,],weights = 1/(1+art$sd[filt])),col='gray',lwd=3)
      filt = real$TP>2
      abline(s.art<-lm(TP~lesd,data=real[filt,],weights = 1/(1+real$sd[filt])),col='dodgerblue4')
      
      r.real=summary(s.real)
      r.art=summary(s.art)
      
      if(sum(art$TP>2)>2) sum.real[idx]=r.real$coefficients[2,1]
      else sum.real[idx]=0
      sum.art[idx]=r.art$coefficients[2,1]    
    }
    
    y=log(c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7))
    y.lab = -7:1
    
    axis(side=1, at=(y), labels=y.lab,las=1,cex.axis=1.5)
    axis(side=2, at=c(2,2.5,3,3.5,4),las=1,cex.axis=1.5)
  }
 
}

db = db0
mtext(side=1,TeX(paste0('log$_{10}$ ','predator size (m)')),outer=T,line=3,adj=0.5,cex=2)
mtext(side=2,'Trophic level',outer=T,line=2,adj=0.5,las=0,cex=2)

###############
