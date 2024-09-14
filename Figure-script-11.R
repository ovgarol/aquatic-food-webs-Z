library(mclust)
the.breaks = seq(-7,7,by=.25/sqrt(3))
db.all = db.all[!(db.all$class%in%c('Copepoda')&round(db.all$esd)==425),]
db.all$group[db.all$group=='Invertebrate'&db.all$opt<(20e-5*1e6)] = 'Invertebrate (II)'

par(mfrow=c(1,3),las=1)


res = log(db.all$opt)-(r.offset+log(db.all$esd)-gamma*log(db.all$esd)**2+db.all$mean.m)+r.offset
res = res[db.all$group%in%unique(db.all$group)[1:6]]
res0 = res
res = res0[res0>-3.5 & res0<3.5]
db.all0 = db.all[(res0>-3.5 & res0<3.5),]
BIC <- mclustBIC(res,G=1:5)
plot(BIC,main='')
title(main="Optimal number of components",adj=0,line=0.5)
text(1,-1640,'a',adj=c(0,1),cex=2,font=2)

summary(BIC)

hh=hist(res0,breaks=the.breaks,freq=F,border='NA',col=alpha('#f0d9f0',0.99),main="",xlab='residual')
title(main="histogram of residuals, r'",adj=0,line=0.5)
text(-6.5,0.3,'b',adj=c(0,1),cex=2,font=2)

he.breaks=hh$mids

nn.list=3



for(nn in nn.list){
  #nn = 3
  ff = 0 * the.breaks
  gmm_model <- Mclust(res,G=nn,modelNames = 'E')  
  summary(gmm_model)
  
  for(i in 1:nn)lines(lwd=2,the.breaks,gmm_model$parameters$pro[i]*dnorm(the.breaks,gmm_model$parameters$mean[i],sqrt(gmm_model$parameters$variance$sigmasq)))
  for(i in 1:nn) abline(v=gmm_model$parameters$mean[i])
  
  print(gmm_model$parameters$mean)
  
  cluster_assignments <- predict(gmm_model)$classification
  
  for(i in 1:nn)ff=ff+gmm_model$parameters$pro[i]*dnorm(the.breaks,gmm_model$parameters$mean[i],sqrt(gmm_model$parameters$variance$sigmasq))
  #lines(the.breaks,ff,col="red")
  print(c(nn,gmm_model$BIC))
  print(mean(gmm_model$uncertainty))
}


sqrt(gmm_model$parameters$variance$sigmasq)
3/sqrt(3)

if(F){
  plot(db.all0$esd,db.all0$opt,col=c('darkmagenta','gray','darkcyan')[cluster_assignments],log='xy',pch=19)
  plot(db.all0$esd,db.all0$opt,log='xy')
  plot(db.all$esd,db.all$opt,log='xy')
  #plot(gmm_model$parameters$mean,gmm_model$parameters$pro)
  boxplot(res~cluster_assignments)
  plot(log(db.all0$esd)-db.all0$mean.D+db.all0$mean.m,res,bg=db.all0$col, pch=c(25,23,24)[cluster_assignments],log='',col='white')
  
}

db.all0$col2 = cluster_assignments
#"Jellyfish""Invertebrate" "Fish" "Invertebrate (II)" "Plankton" "Mammals"
db.all0=db.all0[db.all0$group%in%unique(db.all0$group)[c(1,3,4,6)],]
#db.all0=db.all[db.all$group%in%unique(db.all$group)[c(1,3,5,6,4)],]

plot(log(db.all0$esd)-db.all0$mean.D,+log(db.all0$opt)-db.all0$mean.D-db.all0$mean.m,col=alpha(c('darkmagenta','gray','darkcyan'),0.5)[db.all0$col2],log='',pch=19,
     xlab=TeX('relative predator size, $D-\\bar{D}$'),
     ylab=TeX('relative OPS, $D_{opt}-\\bar{D}-m$'))
title(main="Z-pattern",adj=0,line=0.5)
text(-3,2.5,'c',adj=c(0,1),cex=2,font=2)

abline(a=r.offset,b=1-0*gamma*2)
unique(db.all0$col)
unique(db.all0$group)

lines(smooth.spline(spar=.99,log(db.all0$esd)-db.all0$mean.D,log(db.all0$opt)-db.all0$mean.D-db.all0$mean.m-0*r.offset),lwd=4,col='darkgray')

db.all1=db.all0[db.all0$col2==1,]
lines(smooth.spline(spar=.99,log(db.all1$esd)-db.all1$mean.D,log(db.all1$opt)-db.all1$mean.D-db.all1$mean.m-0*r.offset),lwd=4,col='darkmagenta')

db.all1=db.all0[db.all0$col2==3,]
lines(smooth.spline(spar=.99,log(db.all1$esd)-db.all1$mean.D,log(db.all1$opt)-db.all1$mean.D-db.all1$mean.m-0*r.offset),lwd=4,col='darkcyan')
abline(lm(log(db.all1$esd)-db.all1$mean.D,log(db.all1$opt)-db.all1$mean.D-db.all1$mean.m-r.offset))
