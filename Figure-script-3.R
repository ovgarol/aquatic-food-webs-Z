error.limit = 2/sqrt(3) # feeding kernell width
SP = .5/sqrt(3) # critical distance for vertical clustering

####
# FUNCTION DEFINITIONS
####

my.polygon = function(x,y,ef=ef0,has.borders=F,...){
  ef0=0.5 
  efy=0.025
  polygon(c((1-ef)*min(x),(1-ef)*min(x),(1+ef)*max(x),(1+ef)*max(x)),c((1-efy)*min(y),(1+efy)*max(y),(1+efy)*max(y),(1-efy)*min(y)),...)
  if(has.borders)polygon(c((1-ef)*min(x),(1-ef)*min(x),(1+ef)*max(x),(1+ef)*max(x)),c(min(y),max(y),max(y),min(y)),border='#525252')
} 

my.range = function(x)return(max(x)-min(x))

#ma = function(x,n){as.vector(filter(x, rep(1 / n, n), sides = 2))}
ma = function(...)runmean(...)

temp.pal = hcl.colors(5,palette='Zi')
temp.pal = c(temp.pal[1:2],temp.pal[2:5])

####
## Figure 2
####
pdf('./figures/BANDS.pdf',family='Helvetica',8,5)
par(mfrow=c(2,3),family='Helvetica',mai=.1+c(0.1,0.1,0.2,0.1),oma=c(3.5,3.5,.0,.0),las=1)
i = 1 # auxiliary iterator
## separating invertebrates in two groups
db.all$group[db.all$group=='Invertebrate'&db.all$opt<(20e-5*1e6)] = 'Invertebrate (II)'

for(kk in unique(db.all$group)[c(5,4,2,1,3,6)]){
  db.all.s = db.all[db.all$group==kk&!(db.all$class%in%c('Copepoda')&round(db.all$esd)==425),]
  db.all.s$opt.2=round(db.all.s$opt/10^floor(log10(db.all.s$opt)),digits = 2)*10^floor(log10(db.all.s$opt))

  the.breaks = seq(min(log(c(db.all.s$opt,db.all.s$ops.2)))-.10,max(log(c(db.all.s$opt,db.all.s$ops.2)))+.3,by =.5/sqrt(3))
  
  x.range = range(db.all.s$esd)*c(0.9,1.3)
  xhist=hist(log(db.all.s$opt),breaks=the.breaks,probability = T,border=NA,
             #col=alpha(the.spec.color[1:length(the.breaks)],0.15),
             col=alpha('#f0d9f0',0.75),
             bty='o',plot=F,
             main=paste('horizontal bands for',kk),
             xaxt='n',
             xlab=TeX('relative prey log-size D$_{opt}$-$\\bar{D}$'))
  
  db.all.s$opt.min = db.all.s$opt
  db.all.s$x = cut(log(db.all.s$opt),the.breaks,labels=F)

  plot(db.all.s$esd,db.all.s$opt,log='yx',
            col=alpha(temp.pal[i],0.75),
            pch=19,type='n',
            xlim=x.range,
            #ylim=exp(range(the.breaks)),
            xlab='',xaxt='n',
            ylab='',yaxt='n')
  text(min(db.all.s$es),max(db.all.s$opt),letters[i],font=2,cex=1.5,adj=c(0,1))
  
  title(main=c('unicellular','planktonic invertebrates','invertebrates','jellyfish','fish','mammals')[i],adj=0,line=0.33)
  lines(x.range,x.range,lty=2) # add 1:1 dashed line
  #points(db.all.s$esd,db.all.s$ops.0) # to plot using allometric model
  #abline(v=exp(db.all.s$mean.D)) # add line for mean size
  #abline(h=exp(db.all.s$mean.D+db.all.s$mean.m)) # add line for mean OPS

  ## plot allometric rule
  if(T){
    x = log(seq(0.65*min(db.all.s$esd),1.6*max(db.all.s$esd),length.out=20))#log(10^(-60:10 * 0.1)*1e6 )
    y = exp(x+mean(log(db.all.s$opt/db.all.s$esd)) - gamma*(x-mean(log(db.all.s$esd))-3)**2)
    x = exp(x)
    lines(x,y,col='gray',lwd=3)
    polygon(c(x,rev(x)),c(1/((1+SP)*error.limit)*y,rev((1+SP)*error.limit*y)),border=NA,col=alpha('gray',0.2),lty=1,lwd=.5) 
  }
  
  for(k in round(length(db.all.s$esd)/2):2){
    d = dist(log(subset(db.all.s,select=c(opt))), method = 'euc')
    hclust_avg = hclust(d, method = 'average')
    cut_avg = cutree(hclust_avg, k=k)
    db.all.s$kk = cut_avg
    for(jj in 1:length(db.all.s$name)) db.all.s$opt.c[jj] = mean(db.all.s$opt[db.all.s$kk==db.all.s$kk[jj]] )
    min.dist = min(dist(log(db.all.s$opt.c))[!dist(log(db.all.s$opt.c))%in%c(0)])
    if(min.dist>SP) break
  }
  
  
  #for(jj in 1:(k-1)) if(jj%in%unique(db.all.s$kk)) if(abs(mean(log(db.all.s$opt[db.all.s$kk == jj])) - mean(log(db.all.s$opt[db.all.s$kk == (jj+1)]))) < error.limit) db.all.s$kk2[db.all.s$kk == (jj+1)] = jj
  #for(jj in 1:k) if(((log(max(db.all.s$esd[cut_avg==jj]))-log(min(db.all.s$esd[cut_avg==jj])))>(1*error.limit))&sum(cut_avg==jj)>=5)my.polygon(range(db.all.s$esd[cut_avg==jj]),range(db.all.s$opt.min[cut_avg==jj]),col=alpha('#f0d9f0',0.5),border='#f099f0',lwd=1)
  
  ## height as range
  #for(jj in 1:k) if(sum(cut_avg==jj)>=5) if(my.range(log(db.all.s$esd[cut_avg==jj]))>(1.*error.limit)) if(sd(log(db.all.s$opt[cut_avg==jj]))<(1.*SP*error.limit)) my.polygon(range(db.all.s$esd[cut_avg==jj]),range(db.all.s$opt.min[cut_avg==jj]),col=alpha('green',0.25),border='gray',lwd=1)
  ## height as sd
  for(jj in 1:k) if((my.range(log(db.all.s$opt[cut_avg==jj]))<(1.96*SP*error.limit))&(my.range(db.all.s$esd[cut_avg==jj])>(1.*error.limit))&sum(cut_avg==jj)>=6)my.polygon(range(db.all.s$esd[cut_avg==jj]),exp(mean(log(db.all.s$opt.min[cut_avg==jj]))+c(-1,1)*sd(log(db.all.s$opt.min[cut_avg==jj])))+0*range(db.all.s$opt.min[cut_avg==jj]),col=alpha('#f0d9f0',0.9),border=NA,lwd=1)
  #for(jj in 1:k) my.polygon(range(db.all.s$esd[cut_avg==jj]),range(db.all.s$opt.min[cut_avg==jj]),col=alpha('blue',0.05),border=NA,lwd=.1)
  
  #abline(h=unique(db.all.s$opt.2),col=alpha('gray',0.5))
  
  opt.2=round(seq(min(db.all.s$opt),max(db.all.s$opt))/10^floor(log10(seq(min(db.all.s$opt),max(db.all.s$opt)))),digits = 0)*10^floor(log10(seq(min(db.all.s$opt),max(db.all.s$opt))))
  opt.2=unique(opt.2)
  ## to plot censored data
  #abline(h=unique(opt.2),col=alpha('gray',0.5))
  
  II = 1:k
  for(jj in 1:k)II[jj] = ((log(max(db.all.s$esd[cut_avg==jj]))-log(min(db.all.s$esd[cut_avg==jj])))>.5*error.limit)&sum(cut_avg==jj)>=5
  test.ds = db.all.s[db.all.s$kk%in%c((1:k)*II),]
  print(paste(kk,length(db.all.s$esd),'============================================'))
  print(summary(lm(log(opt)~as.factor(kk),data=test.ds)))
  
  
  for(j in 1:length(db.all.s$esd)) lines(db.all.s$esd[j]*c(1,1),c(db.all.s$opt.c[j],db.all.s$opt[j]),lwd=.5,col=alpha(temp.pal[i],1))
  points(db.all.s$esd,db.all.s$opt.c,bg=alpha(temp.pal[i],1), pch=21,type='p',col='white',cex=1.,lwd=0.25)
  points(db.all.s$esd,db.all.s$opt,col=alpha(temp.pal[i],1), pch='_',type='p',cex=1.,lwd=1)
  
  y=c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)
  y.lab = c(TeX('10$^{-7}$'),TeX('10$^{-6}$'),TeX('10$^{-5}$'),TeX('10$^{-4}$'),TeX('10$^{-3}$'),TeX('10$^{-2}$'),TeX('10$^{-1}$'),TeX('10$^{0}$'),TeX('10$^{1}$'))
  axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
  axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)
  
  i = i+1
  }

mtext(side=1,'predator size (m)',outer=T,line=2,cex=1.)
mtext(side=2,'optimal prey size (m)',outer=T,las=0,line=2,cex=1.)

dev.off()
