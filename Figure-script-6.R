error.limit = 2/sqrt(3) # width of feeding kernel

fit.net.model = function(db,w,do.real=F){
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
  db.sp = na.omit(db.sp)
  db.sp$min.esd = as.numeric(db.sp$min.esd)
  db.sp$max.esd = as.numeric(db.sp$max.esd)
  db.sp$s = as.numeric(db.sp$s)
  
  ## using db.sp
  db.mm = db.sp # using the mechanistic food web model
  #db.mm$tag = paste(db.mm$group,db.mm$type,sep=' - ')
  db.mm$mean.D = 0.5*(log(db.mm$min.esd)+log(db.mm$max.esd)) # calculating median log size
  db.mm$mean.OPS = db.mm$mean.D-3
  for(i in 1:length(db.mm$tag)) db.mm$mean.OPS[i] = Min.Spec(db.mm$mean.D[i],db.mm$tag[i],db.mm)###TEMPORAL
  
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
    net = simplify(net, remove.multiple = T, remove.loops = F)
    art.net[[n]] = net
  }
  
  ##real web
  if(do.real){
    meta.data <- data.frame(
      size = c(db$con.mass.mean.g.,db$res.mass.mean.g.)
    )
    
    relations = data.frame(
      from = db$res.taxonomy,
      to = db$con.taxonomy
    )
    
    a = matrix(nrow=length(meta.data$the.names),ncol=length(meta.data$the.names))
    #a = scale(a, center=FALSE, scale=colSums(a)+1e-9)
    net.0 = graph_from_data_frame(relations,vertices=meta.data$the.names, directed = T)
    net.0 = simplify(net.0, remove.multiple = T, remove.loops = F)
  }else{net.0=NA}
  
  return(list(net.0,net,accuracy))

}

calc.connectance = function(net){
  #net = simplify(net,remove.multiple = T,remove.loops = F)
  return(length(E(net))/length(V(net))^2)
} 

par(bty='o',mfrow=c(3,3),mai=0.051+c(.2,.2,.2,.2),oma=c(5,5,1,0))

#what.to.plot = 1:222
what.to.plot = c(11,5,3,19,2,18,4,6,22,13,1,16,14,15,20,17,21,12) # sorted by ecosystem type

list.w = list()
list.idx = list()
list.acc = list()
list.net = list()

i = 1
possible.w = seq(0.5,5,by=0.1) #c(1,1.5,2,2.5,3,3.5,4,4.5,5)
#possible.w = c(1,1.5,2,2.5,3,3.4,3.5,4,4.5,5.)#c(3.4,2.1) #shorter version

for(w in possible.w){
  for(idx in what.to.plot){
    db = db0
    #db = db[db$foodweb.name == unique(db$foodweb.name)[idx],] 
    db = db[db$study.site == unique(db$study.site)[idx],] 
    db = db[db$con.movement.type=='swimming',] # restricting to swimming predators
    db = db[db$con.mass.mean.g.!=-999 & db$res.mass.mean.g.!=-999, ] # removing NA values
    db = db[!is.na(db$autoID),] # removing NA values
    #db= db[!(db$res.taxonomy=='Calanus finnmarchicus'|db$res.mass.mean.g == 2.099800e-04),]
    
    net = fit.net.model(db,w)

    list.w[[i]] = w
    list.acc[[i]] = net[[3]]
    list.idx[[i]] = idx
    list.net[[i]] = net[[2]]
    i = i+1
  }
  
}


list.acc=unlist(list.acc)
list.w=unlist(list.w)
list.idx = unlist(list.idx)

real.net = list()
i = 1
for(idx in what.to.plot){
  db = db0
  #db = db[db$foodweb.name == unique(db$foodweb.name)[idx],] 
  db = db[db$study.site == unique(db$study.site)[idx],] 
  db = db[db$con.movement.type=='swimming',] # restricting to swimming predators
  db = db[db$con.mass.mean.g.!=-999 & db$res.mass.mean.g.!=-999, ] # removing NA values
  db = db[!is.na(db$autoID),] # removing NA values
  db= db[!(db$res.taxonomy=='Calanus finnmarchicus'|db$res.mass.mean.g == 2.099800e-04),]
  
  net = fit.net.model(db,3.5,do.real=T)

  real.net[[i]] = net[[1]]
  i = i+1
}

####
## beta as function of number of feeding guilds (not used)
####
par(mfrow=c(1,1),las=1)

acc.mean=aggregate(list.acc,by=list(list.w),mean)$x
acc.min=aggregate(list.acc,by=list(list.w),min)$x
acc.max=aggregate(list.acc,by=list(list.w),max)$x

tt = 1:length(list.acc)
for(i in 1:length(tt)) tt[i] = length(V(list.net[[i]]))
tg.mean=aggregate(tt,by=list(list.w),mean)$x
tg.max=aggregate(tt,by=list(list.w),max)$x

tt = 1:length(list.acc)
for(i in 1:length(tt)) tt[i] = calc.connectance(list.net[[i]])

vv = tt *0
cc = tt * 0

for(i in 1:length(tt)){
  net = list.net[[i]]
  vv[i] = length(V(net))
  cc[i] = calc.connectance(net)
}

beta = NULL
beta.sd = NULL
plot(cc~vv,log='xy',pch=1,col=alpha('lightgray',0.5))
for(w in possible.w){
  points(cc[list.w==w]~vv[list.w==w],pch=19)
  #abline(s0<-lm(log10(cc)[list.w==w]~log10(vv)[list.w==w]),col=w*2)
  s0<-lm(log(cc)[list.w==w][cc[list.w==w]>0]~log(vv)[list.w==w][cc[list.w==w]>0])
  beta= c(beta,s0$coefficients[2])  
  beta.sd= c(beta.sd,summary(s0)$coefficients[2,2])  
}

####
## Figure A5
####

par(mfrow=c(1,3),las=1)

plot((possible.w),possible.w,pch=19,type='n',ylim=c(0,1),xlab='',log='',yaxt='n')
title(main=TeX('number of trophic guilds $N$, and accuracy',bold=T),adj=0,line=0.66)
polygon(c(possible.w,rev(possible.w)),(c(acc.min,rev(acc.max))-60)/(100-60),col=alpha('royalblue',0.15),border = NA)
lines(possible.w,(tg.mean-min(tg.mean))/(max(tg.mean)+3-min(tg.mean)),col='black',lwd=3)
lines(possible.w,(tg.max-min(tg.mean))/(max(tg.mean)+3-min(tg.mean)),col='black',lwd=1)
lines(possible.w,(acc.mean-60)/(100-60),col='royalblue',lwd=3,type='l',pch='|')
axis(side=4, at=c(0,0.25,0.5,0.75,1),label=c(60,70,80,90,100),col='royalblue')
axis(side=2, at=c(0,0.25,0.5,0.75,1),label=c(round(min(tg.mean)),10,15,20,3+round(max(tg.mean))),col='black')

plot(0:5,1:6,ylim=c(-.8,-0.1),xlim=c(.5,5),type='n',
     xlab='feeding guild width, log-units',
     ylab='')
title(main=TeX('scaling exponent, $\\beta$',bold=T),adj=0,line=0.66)
#polygon(c(0,0,6,6),-0.35+2*0.11*c(1,-1,-1,1),col='gray95',border=NA)
polygon(c(0,0,6,6),-0.35+0.11*c(1,-1,-1,1),col='gray85',border=NA)
text(3,-0.3,TeX('real food web\n$\\beta \\pm error$',bold=T),cex=1.5,col='white')

abline(h=-0.35,lwd=3,col='white')
for(i in 1:length(beta))lines(possible.w[i]+c(0,0),beta[i]+c(-1,1)*beta.sd[i])
points(possible.w,beta,pch=19,type='p',cex=0.5)
#abline(v=3.4)

vv.real = 1:length(real.net) *0
cc.real = 1:length(real.net) * 0
for(i in 1:length(real.net)){
  net = real.net[[i]]
  vv.real[i] = length(V(net))
  cc.real[i] = calc.connectance(net)
}

plot(vv.real,cc.real,log='xy',ylim=c(1e-2,1),xlim=c(1,500),pch=1,col=alpha(eco.col,0.8),lwd=2)
title(main=TeX('connectance, $C = L/N^2$',bold=T),adj=0,line=0.66)
i = 1
for(w in 3.4*10){
  points(cc[10*list.w==w]~vv[10*list.w==w],pch=19,col=alpha(eco.col,0.8))
  i = i +1
}
abline(lm(log10(cc.real)~log10(vv.real)),lty=2)
cc.d = cc.real
vv.d = vv.real
s0=lm(log(cc.d)~log(vv.d))
summary(s0)

w=34
abline(lm(log10(cc[10*list.w==w])~log10(vv[10*list.w==w])))
cc.d = cc[10*list.w==w]
vv.d = vv[10*list.w==w]
s1=lm(log(cc.d)~log(vv.d))
summary(s1)
anova(s0,s1)

comparison = data.frame(vv=c(vv.d,vv.real),cc=c(cc.d,cc.real),group=c(rep('S',18),rep('O',18)))
s0 = lm(log(cc)~log(vv)*group-1,data=comparison)
summary(s0)

z.beta = unname(abs(beta+0.35)/0.35)

score = (z.beta+(100-acc.mean)/100+1*(log10(tg.mean)-.5))/3
score = (1-(tg.mean-min(tg.mean))/(max(tg.mean)-min(tg.mean)) + 0*(1-(z.beta-min(z.beta))/(max(z.beta)-min(z.beta))) + (acc.mean-min(acc.mean))/(max(acc.mean)-min(acc.mean)))/2

######
##  Plot quality score (not used figure)
#####
plot((possible.w),score,pch=19,type='n',ylim=c(0.4,1),log='',xlab='feeding guild log-size w')
title(main=TeX('quality scores, $Q$',bold=T),adj=0,line=0.66)
lines((possible.w),score,lwd=3,type='l')
abline(v=3.4)

##
  

