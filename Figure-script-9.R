
AIC = function(k,x,y){
  n = length(x)
  sigma.hat.2=sum((x-y)**2,na.rm=T)/(n-k-1)
  dd = 2*k*n/(n-k-1) + n*log(sigma.hat.2)
  return(dd)
}


##
# Indices for each model are
# ops.0 for specialization model
# ops.1 for no stiffness mean a is used instead
# ops.2 for no specialization no rotation
# ops.3 for size-only model 
# ops.4 for no m but yes s
##


db.all$mean.D = NA + db.all$esd
db.all$mean.m = NA + db.all$esd
db.all$s = NA + db.all$esd
db.all$a = NA + db.all$esd
db.all$ops.0 = NA + db.all$esd
db.all$ops.1 = NA + db.all$esd
db.all$ops.2 = NA + db.all$esd
db.all$ops.3 = NA + db.all$esd
db.all$ops.4 = NA + db.all$esd

for(j in 1:length(db.all$tag)){
  temp = spec[spec$tag==db.all$tag[j],]
  
  db.all$mean.D[j] = temp$mean.D[1]
  db.all$mean.m[j] = temp$mean.m[1]
  db.all$mean.m.sd[j] = mean(temp$mean.m.sd)
  db.all$s[j] = temp$s[1]
  db.all$a[j] = temp$a[1]
  db.all$stiffness[j] = temp$stiffness[1]
  db.all$ops.0[j] = exp(-0.33-r.offset+OPS(log(db.all$esd[j]),db.all$mean.D[j],db.all$mean.m[j],db.all$s[j],a=db.all$stiffness[j]))
  db.all$ops.1[j] = exp(-0.33-r.offset+OPS(log(db.all$esd[j]),db.all$mean.D[j],db.all$mean.m[j],db.all$s[j],a=mean(db.all$stiffness)))
  db.all$ops.2[j] = exp(-0.33-r.offset+OPS(log(db.all$esd[j]),db.all$mean.D[j],db.all$mean.m[j],0,a=0,gamma=gamma))
  db.all$ops.3[j] = exp(OPS(log(db.all$esd[j]),db.all$mean.D[j],0,0,a=0,gamma=gamma))
  db.all$ops.4[j] = exp(-0.33+OPS(log(db.all$esd[j]),db.all$mean.D[j],0,db.all$s[j],a=db.all$stiffness[j]))
  db.all$ops.5[j] = exp(-0.33-r.offset+OPS(log(db.all$esd[j]),db.all$mean.D[j],db.all$mean.m[j],sign(round(db.all$s))[j],a=mean(db.all$stiffness)))
  
}

a0=AIC((5*5)+3,log(db.all$opt),log(db.all$ops.0)) # specialization full
a1=AIC((4*5)+3+1,log(db.all$opt),log(db.all$ops.1)) # no-stiffness
a2=AIC((3*5)+3,log(db.all$opt),log(db.all$ops.2)) # no-specialization
a3=AIC(3,log(db.all$opt),log(db.all$ops.3)) # size-only
a4=AIC((4*5)+3+1,log(db.all$opt),log(db.all$ops.4)) # no m no Z displacement
a5=AIC((4*5)+3+1,log(db.all$opt),log(db.all$ops.5)) # fixed s, no rotations

par(mfrow=c(2,3))
plot(db.all$opt,db.all$ops.0,log='xy',xlab='',
     main = 'specialization = Full model',adj=0,
     ylab='',pch=19,col=alpha('gray',0.7),lwd=1,cex=0.5,yaxt='n',xaxt='n',
     xlim=size.range,ylim=size.range)
abline(a=0,b=1)
text(1,1e7,paste0('AIC = ',round(a0)),font=2,adj=c(0,1),cex=1.5)
text(1e7,1e0,'a',font=2,adj=c(1,0),cex=1.5)
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)
mtext(side=2,'calculated OPS (m)',las=0,line=-2,cex=0.8,outer=T)
mtext(side=1,'observed OPS (m)',las=0,line=-2,cex=0.8,outer=T)

plot(db.all$opt,db.all$ops.1,log='xy',xlab='',
     main = 'no stiffness = no Z-scaling',adj=0,
     ylab='',pch=19,col=alpha('gray',0.7),lwd=1,cex=0.5,yaxt='n',xaxt='n',
     xlim=size.range,ylim=size.range)
abline(a=0,b=1) 
text(1,1e7,round(a1),font=2,adj=c(0,1),cex=1.5)
text(1e7,1e0,'b',font=2,adj=c(1,0),cex=1.5)
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)

plot(db.all$opt,db.all$ops.5,log='xy',xlab='',
     main = 'no rotation',adj=0,
     ylab='',pch=19,col=alpha('gray',0.7),lwd=1,cex=0.5,yaxt='n',xaxt='n',
     xlim=size.range,ylim=size.range)
abline(a=0,b=1)
text(1,1e7,round(a5),font=2,adj=c(0,1),cex=1.5)
text(1e7,1e0,'c',font=2,adj=c(1,0),cex=1.5)
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)

plot(db.all$opt,db.all$ops.4,log='xy',xlab='',
     main = 'no m = no Z-displacement',adj=0,
     ylab='',pch=19,col=alpha('gray',0.7),lwd=1,cex=0.5,yaxt='n',xaxt='n',
     xlim=size.range,ylim=size.range)
abline(a=0,b=1)
text(1,1e7,round(a4),font=2,adj=c(0,1),cex=1.5)
text(1e7,1e0,'d',font=2,adj=c(1,0),cex=1.5)
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)

plot(db.all$opt,db.all$ops.2,log='xy',xlab='',
     main = 'no specialization = only PFG differences',adj=0,
     ylab='',pch=19,col=alpha('gray',0.7),lwd=1,cex=0.5,yaxt='n',xaxt='n',
     xlim=size.range,ylim=size.range)
abline(a=0,b=1)
text(1,1e7,round(a2),font=2,adj=c(0,1),cex=1.5)
text(1e7,1e0,'e',font=2,adj=c(1,0),cex=1.5)
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)

plot(db.all$opt,db.all$ops.3,log='xy',xlab='',
     main = 'size-only model',adj=0,
     ylab='',pch=19,col=alpha('gray',0.7),lwd=1,cex=0.5,yaxt='n',xaxt='n',
     xlim=size.range,ylim=size.range)
abline(a=0,b=1)
text(1,1e7,round(a3),font=2,adj=c(0,1),cex=1.5)
text(1e7,1e0,'f',font=2,adj=c(1,0),cex=1.5)
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)





