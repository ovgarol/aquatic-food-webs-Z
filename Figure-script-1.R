setwd('/home/og/Documents/articles/functionalFeedingMode/')

library(scales)
library(latex2exp)

####
# DATA INPUT
####
## read db.all directly
db.all = read.csv("repo/database.csv")
db.all$tag = as.factor(paste(db.all$group,'-',db.all$type))

## read database of limits of each feeding group constructed by data
db.mm = read.csv("repo/minimal_model.csv",comment.char = '#')
db.mm$tag = paste(db.mm$group,db.mm$type,sep=' - ')
db.mm$mean.D = 0.5*(log(db.mm$min.esd)+log(db.mm$max.esd))


####
# FUNCTION DEFINITIONS
####

as.numeric.factor = function(...) as.numeric(droplevels(as.factor(...)))

OPS = function(x,x.bar,m.bar,s,a=1,gamma=0.011){
  alpha = exp(-1*(a*s)^2)
  return(alpha*x + r.offset - gamma*(x-x.bar-3)^2 + m.bar + s + (1.-alpha)*x.bar)
}

my.lines = function(x,y,lwd=1,...){
  lines(x,y,col='white',lwd=2*lwd)
  lines(x,y,lwd=lwd,...)
}

col.white = 'white'
col.black = 'black'

median=function(...) stats::median(...)

abs = function(x)  base::abs(x)


bounding.box = function(esd.min,esd.max,kw,mean.m,s,mean.D,stiffness,col='gray',...){
  qq.min<-1e0*exp(exp(-s**2)*(log(esd.min)-mean.D)+mean.D+mean.m+s/stiffness-0*gamma*(log(esd.min)-mean.D+3)**2-0.45)
  qq.max<-1e0*exp(exp(-s**2)*(log(esd.max)-mean.D)+mean.D+mean.m+s/stiffness-0*gamma*(log(esd.max)-mean.D+3)**2-0.45)
  polygon(c(esd.min,esd.max,esd.max,esd.min), c(qq.min/kw,qq.max/kw,qq.max*kw,qq.min*kw),border=alpha(col,1),lwd=1,col=alpha(col,0.15),...)
  #lines(c(esd.min,esd.max),c(qq.min,qq.max))
  #points(exp(0.5*(log(esd.min)+log(esd.max))), exp(0.5*(log(qq.min)+log(qq.max))),pch='*',cex=3)
  return(lm(log(c(qq.min,qq.max))~log(c(esd.min,esd.max))))
}


minimal.model = function(esd.min,esd.max,s,type=1,...){
  esd.bar = 0.5*(log(esd.min)+log(esd.max))
  m.bar = -exp(-0.77+0.14*esd.bar) # Eq(c)
  #m.bar = 0.13-0.23*esd.bar # Eq(c)
  a = exp(-2.92 + 0.26*esd.bar) # Eq(b)
  rr = bounding.box(esd.min,esd.max,2*2/sqrt(3),m.bar,s*sqrt(a),esd.bar,sqrt(a),...)
  return(rr)
} 

#####
# CHECKING CONSISTENCY   
#####

#my.pal = colorRampPalette(c("orange","limegreen", "yellow3","tomato","dodgerblue")) # alternative colorpalette
my.pal = colorRampPalette(hcl.colors(5,palette='Zi')[c(4,2,3,5,1)])
palette(my.pal(5))

## counting guilds, groups and tags
aggregate(db.all$name,by=list(db.all$tag), FUN = length)
aggregate(db.all$name,by=list(db.all$group), FUN = length)
aggregate(db.all,by=list(db.all$tag), FUN = mean)

levels(as.factor(db.all$type))
levels(as.factor(db.all$group))
unique(db.all$tag)


#####
# SPECIALIZATION CALCULATIONS
#####

## calculation of structural coefficient \gamma Eq(SM1)
s = lm(log(opt)-log(esd)~I(log(esd)^2),data=db.all)
summary(s)
gamma = -s$coefficients[2]

displ = 0 # displacement in specialization trait (not used)

## calculation of PFG and guild specific traits
spec = NULL
spec$tag = unique(db.all$tag)
dummy.list = 0 * 1:length(spec$tag)
spec$group =  NA + dummy.list
spec$esd = NA + dummy.list
spec$esd.min = NA + dummy.list
spec$esd.max = NA + dummy.list
spec$opt = NA + dummy.list
spec$mean.D = NA + dummy.list
spec$mean.m = NA + dummy.list
spec$mean.m.sd = NA + dummy.list
spec$m = NA + dummy.list
spec$a = NA + dummy.list
spec$s =0  + dummy.list
spec$a.sd = NA + dummy.list
spec$s.sd = NA + dummy.list
spec$stiffness = NA + dummy.list
spec = as.data.frame(spec)

## iteration over all guilds
i = 1 # auxiliar iterator
for(j in spec$tag){
  if(j =='NA - NA') next # ignoring empty lines
  temp = db.all[db.all$tag==j,] # temp is a subset of all predators in the guild
  spec$group[i] = temp$group[1] # asignation of PFG
  spec$esd[i] = exp(mean(log(temp$esd))) # mean predator size
  spec$opt[i] = exp(mean(log(temp$opt))) # mean optimal prey size
  
  spec$esd.min[i] = min(temp$esd,na.rm=t) # min predator size
  spec$esd.max[i] = max(temp$esd,na.rm=t) # max predator size
  
  spec$opt.min[i] = min(temp$opt,na.rm=t) # min optimal prey size
  spec$opt.max[i] = max(temp$opt,na.rm=t) # max optimal prey size
  
  s = lm(log(opt)~log(esd),data=temp) # regression Eq(SM5)
  #s = lmodel2(log(opt)~log(esd),data=temp,nperm=99) # major axis regression (not used)
  spec$m[i] = s$coefficients[1] # m'_{jk} in Eq(SM5)
  spec$a[i] = s$coefficients[2] # \alpha_{jk} in Eq(SM5)
  
  ## assignation of uncertainties in alpha and m' based on standard error of OLS regression needed for Eq(SM10)
  rrr = summary(s)
  if(length(rrr$coefficients)<5) next
  if(length(temp$tag)>2) spec$a.sd[i] = rrr$coefficients[2,2] # standard error in alpha
  else spec$a.sd[i] = 0
  if(length(temp$tag)>2) spec$s.sd[i] = rrr$coefficients[1,2] # standard error in m'
  else spec$s.sd[i] = 0
  if(length(temp$tag)>2) spec$r2[i] = rrr$r.squared # r2 coefficient for regression
  else spec$s.sd[i] = 0
  
  temp = db.all[db.all$group==temp$group[1],] # selecting all members of a PFG
  spec$mean.D[i] = (mean(log(temp$esd))) -> D.mean # calculation of mean PFG log size Eq(2)
  #s = lm(log(opt)~log(esd),data=temp) # ignoring the structural coefficient \gamma (not used)
  s = lm(log(opt)-log(esd)+gamma*(log(esd)-D.mean-3)^2~1,data=temp) # regression to calculate m_k Eq(SM4)
  
  spec$mean.m[i] = s$coefficients[1] # assigning the value of m_k
  rrr = summary(s)
  spec$mean.m.sd[i] = rrr$coefficients[1,2] # standard error in m_k
  
  i = i+1 # auxiliar iteration 
}

spec$s.sd = sqrt(spec$s.sd) # required for unit consistency
r.offset = mean(log(db.all$opt/db.all$esd),na.rm=T) # defining the prototyping scaling r as the mean log-size ratio
spec$s = spec$m-spec$mean.m-(1-spec$a)*spec$mean.D+gamma*(log(spec$esd)-spec$mean.D-3)**2+0.33 # calculation of specialization Eq(SM6)
spec$s.sd = 0.5*spec$a.sd*spec$mean.D + 0.5*spec$s.sd # error in specialization Eq(SM10)

the.col = (as.numeric(as.factor(spec$group)))*as.numeric(spec$a.sd>0) # vector with the colors of each PFG

############
### Figure A1
############

par(mfrow=c(1,5))
par(bty='o',mai=0.5*c(.8,.8,.2,.2),oma=1+c(0,0,0,0),las=1)
L = seq(-5,5,length.out=300) # range of possible values of specific specialization

## calculation of stiffness
## iteration over all PFGs
jj = 1 # auxiliary iterator
for(i in unique(spec$group)[c(4,2,1,3,5)]){#c(1,3,2,5,4)
  if(is.na(i)) next # ignoring empty lines
  
  plot(a*(a>-0)~I(s*stiffness),data=spec,type='n',xlim=5*c(-1,1),log='',
       ylim=c(0,1.2),
       adj=0,
       xaxt='n',yaxt='n',
       xlab='',
       ylab='',
       main=c('unicellular','invertebrates','jellyfish','fish','mammals')[jj]
  )
  text(-5,1.2,letters[jj],adj=c(0,1),font=2,cex=1.5)
  
  filt = spec$group==i # Selecting PFG
  
  ## calculation of stiffness 
  temp.a = 0.01 * 1:1000 # dummy list of possible vales of a
  temp.b = temp.a # dummy list of errors y 
  for(j in 1:1000) temp.b[j] = sum((exp(-temp.a[j] * abs(spec$s[filt])^2) - spec$a[filt])^2,na.rm = T) # calculating y Eq(SM8)
  spec[filt,]$stiffness = sqrt(temp.a[which.min(temp.b)]) # assigning stiffness that minimizes y
  
  my.lines(L,1*exp(-spec$stiffness[filt][1]**2*abs(L+displ)^2),lwd=3,col='black') # plotting theoretical line
  text(0,1.2,paste('a = ',round(spec$stiffness[filt][1]**2,digits=2)),adj=c(0,1.)) # printing result
  
  for(i in 1:length(spec$tag)) lines(spec$s[filt][i]+c(0,0),spec$a[filt][i]*(spec$a[filt][i]>-100)+0.99*spec$a.sd[filt][i]*c(-1,1)*as.numeric(spec$a[filt][i]>-100),col=the.col[filt][i],lwd=1.5) # error bar in y
  for(i in 1:length(spec$tag)) lines(spec$s[filt][i]+0.5*spec$s.sd[filt][i]*c(1,-1)*as.numeric(spec$a[filt][i]>-100),spec$a[filt][i]*(spec$a[filt][i]>-100)+c(0,0),col=the.col[filt][i],lwd=1.5) # error bar in x
  
  points(a*(a>-0)~I(s),data=spec[filt,],pch=19,col=the.col[filt],cex=1.5,lwd=.5) # plotting data points
  
  axis(side=2,at=c(0,0.5,1))
  axis(side=1,at=-c(-4,0,4))
  
  jj=jj+1 # auxiliary iterator
} 

mtext('specific specialization S',side=1,line=-1,cex=0.75,las=0.5,adj=0.5,outer=T)
mtext(TeX('scaling exponent $\\alpha$'),side=2,line=-1,cex=0.75,las=0,outer=T)



############
### Figure A4 error in OPS calculation
############

##
# Indices for each model are
# ops.0 for specialization model
# ops.2 for size-only model
##


db.all$mean.D = NA + db.all$esd
db.all$mean.m = NA + db.all$esd
db.all$s = NA + db.all$esd
db.all$a = NA + db.all$esd
db.all$ops.0 = NA + db.all$esd
db.all$ops.2 = NA + db.all$esd

for(j in 1:length(db.all$tag)){
  temp = spec[spec$tag==db.all$tag[j],]
  
  db.all$mean.D[j] = temp$mean.D[1]
  db.all$mean.m[j] = temp$mean.m[1]
  db.all$mean.m.sd[j] = mean(temp$mean.m.sd)
  db.all$s[j] = temp$s[1]
  db.all$a[j] = temp$a[1]
  db.all$stiffness[j] = temp$stiffness[1]
  db.all$ops.0[j] = exp(-0.33-r.offset+OPS(log(db.all$esd[j]),db.all$mean.D[j],db.all$mean.m[j],db.all$s[j],a=db.all$stiffness[j]))
  db.all$ops.2[j] = exp(OPS(log(db.all$esd[j]),db.all$mean.D[j],-0.,0,a=db.all$stiffness[j],gamma=gamma))
}

the.col=as.numeric(as.factor(db.all$group))#*as.numeric(spec$a.sd>0)

###
# comparison size-only and specialization
###

par(mfrow=c(2,3),oma=c(2,2,1,1),lend=1)
layout(matrix(c(1, 1, 2, 3), nrow=2, byrow=TRUE),widths=c(1,1,1,1),heights = c(2,2))

if(T){
  abs(1*(db.all$opt-db.all$ops.0)/(db.all$opt))->s0
  abs(1*(db.all$opt-db.all$ops.2)/(db.all$opt))->s1

  the.names = c('fish','invertebrates','jellies','mammals','plankton')
  
  plot(1,xlim=c(-0.2,5.2),ylim=1*c(0.,4),type='n',xaxt='n',yaxt='l',xlab='',ylab='',log='',
       main=c('error comparison'),adj=0)
  abline(h=0,lty=1)
  text(-0.3,4,'a',adj=c(0,1),font=2,cex=1.5)
  
  dis=0.1
  i=0+dis
  x=s0
  y = quantile(x,probs=c(.0,.0,.5,.8,.9),na.rm = T)
  lines(i+c(0,0),c(0,y[4]),lwd=15)
  
  i=0-dis
  x=s1
  y = quantile(x,probs=c(.0,.0,.5,.8,.9),na.rm = T)
  lines(i+c(0,0),c(0,y[4]),lwd=15,col='gray')
  
  if(T){
    for(i0 in 1:5){
      i = i0+dis
      filt=c('Plankton','Invertebrate','Jellyfish','Fish','Mammals')[i0]# unique(db.all$group)[c(1,2,5,3,4)][i]
      
      x= s0[db.all$group==filt]
      y = quantile(x,probs=c(.0,.0,.5,.8,.9),na.rm = T)
      lines(i+c(0,0),c(0,y[4]),lwd=15,col=c(5,2,3,1,4)[i0])
    }
    
    for(i0 in 1:5){
      i=i0-dis
      filt=c('Plankton','Invertebrate','Jellyfish','Fish','Mammals')[i0]# unique(db.all$group)[c(1,2,5,3,4)][i]
      x= s1[db.all$group==filt]
      y = quantile(x,probs=c(.0,.0,.5,.8,.9),na.rm = T)
      lines(i+c(0,0),c(0,y[4]),lwd=15,col='gray')
      print(y)
      print(paste(filt,y[4]))
      
    }
  }
  
  axis(side=1,at=0:5,labels = c('all','unicellular','invertebrates','jellyfish','fish','mammals'),las=0)
  mtext(TeX('relative error in OPS (%)'),side=2,line=3,cex=0.75,las=0)
  
  #axis(side=4,at=c(-2,-1,0,1,2),las=1,cex=0.3)
  #axis(side=2,at=log(c(0.1,0.2,0.5,1,2,5,10)),labels=c(0.1,0.2,0.5,1,2,5,10),las=1)
  abline(h=900,lty=5,col='gray',lwd=5)
  abline(h=900,lty=1,col='white',lwd=4)
  axis.break(axis=2,breakpos=900,pos=NULL,bgcol="white",breakcol="black",
             style="slash",brw=0.02)
}

y=c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)
y.lab = c(' ',TeX('10$^{-6}$'),' ',TeX('10$^{-4}$'),' ',TeX('10$^{-2}$'),' ',TeX('10$^{0}$'),' ')
size.range = c(2e-1,3e6)
size.range = c(5e-1,1e6)

plot(db.all$opt,db.all$ops.2,log='xy',xlab='',
     main = 'size-only model',adj=0,
     ylab='',pch=19,col=the.col,lwd=1,cex=0.5,yaxt='n',xaxt='n',
     xlim=size.range,ylim=size.range)
text(min(size.range),max(size.range),'b',adj=c(0,1),font=2,cex=1.5)
abline(b=1,a=0,lty=2)
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)
mtext(side=2,'calculated OPS (m)',las=0,line=3,cex=0.8)
mtext(side=1,'observed OPS (m)',las=0,line=0,cex=0.8,outer=T)

plot(db.all$opt,db.all$ops.0,log='xy',xlab='',ylab='',
     main = 'specialization model',adj=0,
     pch=19,col=the.col,lwd=1,cex=0.5,yaxt='n',xaxt='n',
     xlim=size.range,ylim=size.range)
text(min(size.range),max(size.range),'c',adj=c(0,1),font=2,cex=1.5)
abline(b=1,a=0,lty=2)
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=NA,las=1,cex.lab=1.)



####
## Figure 2
####

par(mfrow=c(1,3))
par(bty='o',mfrow=c(2,2),mai=c(.8,.8,.2,.2),oma=c(0,0,0,0))
size.range = c(1e0,10e6) # predator size range in microns

## Fig 2a: plot of allometric rule and observed data

db.all$col = as.numeric.factor(db.all$group) # colored by PFG
db.all$ss = db.all$s*db.all$stiffness # normalized specialization

plot(1,xlim=range(db.all$esd),ylim=range(db.all$opt),log='xy',type='n',
     xaxt='n',yaxt='n',
     xlab='predator size (m)',
     ylab='optimal prey size (m)',
     main='allometric rule and data')

lines(size.range,size.range,lty=2)
the.col = as.numeric(as.factor(db.all$group))
the.shape = c(25,23,24)[1 + as.numeric(db.all$ss>-0.33) + as.numeric(db.all$ss >0.33)] # shaped as generalists, small or large specialists
points(opt~esd,data=db.all,pch=the.shape,bg=palette()[the.col],lwd=1,cex=1*(1+0*as.numeric(the.shape>23)),col='white')

y=c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)
y.lab = c(' ',TeX('10$^{-6}$'),' ',TeX('10$^{-4}$'),' ',TeX('10$^{-2}$'),' ',TeX('10$^{0}$'),' ')
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)

v = log(1e0):log(2e7)
my.lines(exp(v),exp(OPS(v,0,0,0,a=0)),col='black',lwd=3) # plotting allometric rule


## Fig 2b: plot of feeding guilds

plot(1,xlim=range(db.all$esd),ylim=range(db.all$opt),log='xy',type='n',
     xaxt='n',yaxt='n',
     xlab='predator size (m)',
     ylab='optimal prey size (m)',
     main='feeding guilds')

lines(size.range,size.range,lty=2)
the.shape = c(25,23,24)[1 + as.numeric(db.all$ss>-0.33) + as.numeric(db.all$ss >0.33)]
#points(opt~esd,data=db.all,pch=the.shape,bg=the.col,lwd=1,cex=1*(1+0*as.numeric(the.shape>23)),col='white') # to include individual data points
#my.lines(exp(v),exp(OPS(v,0,0,0,a=0)),col='black',lwd=1) # to include allometric rule

axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)

## plotting the feeding guilds based on observations
for(i in 1:length(db.mm$group)) if(db.mm$lty[i]==1) minimal.model(db.mm$min.esd[i],db.mm$max.esd[i],db.mm$s[i],col=palette()[c(5,2,3,1,4)][db.mm$color[i]],lty=db.mm$lty[i],lwd=1.5)

## Fig 2c: Z-pattern to be included later

plot.new()

## Fig 2d: mechanistic food web model

plot(1,xlim=range(db.all$esd),ylim=range(db.all$opt),log='xy',type='n',
     xaxt='n',yaxt='n',
     xlab='predator size (m)',
     ylab='optimal prey size (m)',
     main='food web model')
lines(size.range,size.range,lty=2)

#for(i in 1.5+c(0,2,4))for(j in -5:3) minimal.model(10^(i-1),10^(i+1),j*i/5) # old calculation (REMOVE)
#for(i in 1+c(0,2,4))for(j in -5:2) minimal.model(10^(i-2),10^(i+2),j*i/5,col='gray')

## creating data frame for mechanistically defined model
db.sp = data.frame(min.esd=0,max.esd=0,s=0,tag=NA)

ii = 1 # auxiliary iterator
w = 3.5 # universal width for each feeding guild

if(T)for(i in 1:4){ # four size classes
  for(j in -2:1){ # possible values of specialization 
    D.min = exp(ii+abs(j)+0*(j < -1)) # min predator size + offset to specialists
    D.max = exp(ii+w+abs(j)+0*(j < -1)) # max predator size + offset to specialists
    the.D = 0.5*(D.min+D.max) # mean predator size
    the.s = j*0.014*(log(the.D/1e0))**2 # Table 1, Eq(d) 
    #the.s = j*(0.15*(ii+w+1*(j < -1)+abs(j))) # alternative Eq
    the.s = j*(0.17*(log(the.D))) # linear Eq(d)
    
    if(j< -1 & (ii+w > log(1e4) )) next # ignoring additional small prey specialists in larger size classes
    minimal.model(D.min,D.max,the.s,col=c('magenta','darkmagenta','darkgray','darkcyan')[3+j]) # plotting mechanistic model
    db.sp = rbind(db.sp,c(D.min,D.max,the.s,paste0('tag',i,j)),stringsAsFactors = FALSE) # adding model to 
  }
  ii = ii+w
}

## fixing typing in mechanistic model
db.sp = na.omit(db.sp)
db.sp$min.esd = as.numeric(db.sp$min.esd)
db.sp$max.esd = as.numeric(db.sp$max.esd)
db.sp$s = as.numeric(db.sp$s)

axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)


##################
# Figure 3
##################

par(mfrow=c(2,2),las=1)
par(bty='o',mfrow=c(2,2),mai=c(.8,.8,.2,.2),oma=c(0,0,0,0))
high.color = 'red'

##aggregation of PFGs
db.group = aggregate(db.all,by=list(db.all$group),mean) # PFG means
db.sd = aggregate(db.all,by=list(db.all$group),function(...)if(is.numeric(...))sd(...))  # PFG standard deviations as error
db.min = aggregate(db.all,by=list(db.all$group),function(...)if(is.numeric(...))min(...)) # PFG min values
db.max = aggregate(db.all,by=list(db.all$group),function(...)if(is.numeric(...))max(...)) # PFG max values


## Fig 3a: alpha - s 

plot(a*(a>-0)~I(s*stiffness),data=spec,type='n',adj=0.5,
     xlim=4*c(-1,1),ylim=c(0,1.15),
     log='',
     xaxt='n',
     xlab='specialization s',
     ylab=TeX('scaling exponent $\\alpha$'),
     main='variation in OPS scaling')

axis(side=1,at=-c(-4,-2,0,2,4))

L = seq(-10,10,length.out=1000)
my.lines(L,1*exp(-1.0*abs(L+displ)^2),lwd=5,col='black')
the.col=(as.numeric(as.factor(spec$group)))*as.numeric(spec$a.sd>0)

for(i in 1:length(spec$tag)) lines(spec$s[i]*spec$stiffness[i]+c(0,0),spec$a[i]*(spec$a[i]>0)+0.99*spec$a.sd[i]*c(-1,1)*as.numeric(spec$a[i]>-0),col=the.col[i],lwd=1.5)
for(i in 1:length(spec$tag)) lines(spec$s[i]*spec$stiffness[i]+0.5*spec$s.sd[i]*c(1,-1)*as.numeric(spec$a[i]>0),spec$a[i]*(spec$a[i]>0)+c(0,0),col=the.col[i],lwd=1.5)
points(a*(a>-0)~I(s*stiffness),data=spec,pch=19,col=the.col*as.numeric(spec$a.sd>0)*as.numeric(spec$a>-1000),cex=1.5,lwd=.5)
labl = spec

text(-4,1.15,'a',font=2,cex=1.5,adj=c(0,1))

col.list = c(5,2,3,1,4)
tag.list = c('unicellular','invertebrates','jellyfish','fish','mammals')
legend("topright",legend=tag.list,
       horiz=F,
       lwd=0.5,
       lty=NA,
       pt.cex=1.5,bty='n',
       col = col.list,
       pch=19
) 

pseudo.r2 = 1-sum((spec$a*(spec$a>-0) - exp(-(spec$s*spec$stiffness)**2))**2)/sum((spec$a*(spec$a>-0) - mean(spec$a*(spec$a>-0)))**2) # calculation pseudo R2
text(-3,0.8,TeX(paste('$r^2$ = ',round(pseudo.r2,digits=2)) ))

## Fig. 2b: stiffness - mean PFG body size

plot(db.group$esd,db.group$stiffness**2,log='yx',type='n',xaxt='n',
     xlim=c(1e0,1e7),
     ylim=c(0.05,5),
     main = 'effect of specialization in the OPS',
     ylab = 'stiffness a',
     xlab=TeX('mean PFG body size $\\bar{D}$ (m)'))

y=c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)
y.lab = c('',TeX('10$^{-6}$'),'',TeX('10$^{-4}$'),'',TeX('10$^{-2}$'),'',TeX('10$^{0}$'),'')
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
text(1e0,5,'b',font=2,cex=1.5,adj=c(0,1))
#text(db.group$esd,db.group$stiffness**2,db.group$stiffness**2,font=2,cex=1.5,adj=c(0,1))

for(i in 1:length(db.sd$X))lines(c(db.min$esd[i],db.max$esd[i]),db.group$stiffness[i]**2+c(0,0),col=db.group$col[i],lwd=1.5)
for(i in 1:length(db.sd$X))lines(exp(db.group$mean.D[i])+c(0,0),db.group$stiffness[i]**2+0.5*db.group$stiffness[i]**2*c(-1,1),col=db.group$col[i],lwd=1.5)
points(exp(db.group$mean.D),db.group$stiffness**2,col=db.group$col,cex=1.5,pch=15,lwd=.5)


abline(s<-lm(log10(stiffness**2)~log10(exp(mean.D)),data=db.group),col='black',lwd=5)
s = lm(log(stiffness**2)~mean.D,data=db.group) # regression of stiffness and D_k
ss = summary(s)
coef = round(s$coefficients, digits=2)
r2 = round(ss$r.squared,digits = 2)

text(1e0,3, paste('a =',round(exp(coef[1]),digits = 2),'exp(',coef[2],'D)') ,adj=c(0,1))
text(1e0,2, TeX(paste('$r^2$ = ',r2) ),adj=c(0,1))


## Fig. 2c: mean feeding mode and mean PFG body size

plot(db.group$esd,db.group$mean.m,log='x',xaxt='n',type='n',
     xlim=c(1e0,1e7),
     ylim=c(-4,0),
     main='variation in predator-prey ratio',
     ylab='mean feeding mode m',
     xlab=TeX('mean PFG body size $\\bar{D}$ (m)'))

y=c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)
y.lab = c('',TeX('10$^{-6}$'),'',TeX('10$^{-4}$'),'',TeX('10$^{-2}$'),'',TeX('10$^{0}$'),'')
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)

abline(s<-lm(mean.m~log10(exp(mean.D)),data=db.group),lwd=5,col='black')
s = lm(mean.m~mean.D,data=db.group)
ss = summary(s)

for(i in 1:length(db.sd$X))lines(c(db.min$esd[i],db.max$esd[i]),db.group$mean.m[i]+c(0,0),col=db.group$col[i],lwd=1.5)
for(i in 1:length(db.sd$X))lines(exp(db.group$mean.D[i]+c(0,0)),db.group$mean.m[i]+db.group$mean.m.sd[i]*c(1,-1),col=db.group$col[i],lwd=1.5)
points(exp(db.group$mean.D),db.group$mean.m,cex=1.5,pch=15,col=db.group$col,lwd=.5)

coef = round(s$coefficients, digits=2)
r2 = round(ss$r.squared,digits = 2)
text(1e0,-2.75, paste('m  =',coef[1],'+',coef[2],'D') ,adj=c(0,1))
text(1e0,-3, TeX(paste('$r^2$ = ',r2) ),adj=c(0,1))
text(1e0,0,'c',font=2,cex=1.5,adj=c(0,1))

## alternative calculation as log-log regression
s = lm(log(-mean.m)~log(exp(mean.D)),data=db.group)
#lines(xx,-exp(-0.71+0.13*log(xx)))

## Fig. 2d: specialization and mean guild body size

plot(0, xlim=c(1e0,1e7),
     ylab='specialization s',
     xlab='mean guild body size (m)',xaxt='n',
     ylim=c(-4.5,4.5),type='n',log='x')

y=c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)
y.lab = c('',TeX('10$^{-6}$'),'',TeX('10$^{-4}$'),'',TeX('10$^{-2}$'),'',TeX('10$^{0}$'),'')
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)

db.mm = spec
db.mm$ref.esd = exp(0.5*(log(db.mm$esd.min)+log(db.mm$esd.max)))
db.mm$s = db.mm$s*db.mm$stiffness 
db.mm$color = as.numeric.factor(db.mm$group)
text(1e0,4.5,'d',font=2,cex=1.5,adj=c(0,1))

s = lm(abs(s)~log(ref.esd)-1,data=db.mm[c(1,2,3,6,9,10,11,12,14,15),]) # linear model (not used)
summary(s)
#abline(s)

ss = summary(s<-lm(abs(s)~I(log(ref.esd)**2)-1,data=db.mm[c(1,2,3,6,9,10,11,12,14,15),]))
coef = round(s$coefficients, digits=3)
r2 = round(ss$r.squared,digits = 3)
text(1e0,-2, paste('|s| =',coef[1],'D^2') ,adj=c(0,1))
text(1e0,-3, TeX(paste('$r^2$ = ',r2) ),adj=c(0,1))

xx = 10^(0.1*(-3:74))
lines(xx,-coef*(log(xx))**2,lwd=5,col='black')
lines(xx,coef*(log(xx))**2,lwd=5,col='black')
polygon(c(xx,rev(xx)),c(-coef*(log(xx))**2,-5+0*xx),border=NA,col=alpha('gray',0.25))
polygon(c(xx,rev(xx)),c(coef*(log(xx))**2,5+0*xx),border=NA,col=alpha('gray',0.25))

for(i in 1:length(db.mm$group))lines(c(db.mm$esd.min[i],db.mm$esd.max[i]),db.mm$s[i]+c(0,0),col=palette()[db.mm$color[i]],lwd=1.5)
for(i in 1:length(db.mm$group))lines(db.mm$ref.esd[i]+c(0,0),db.mm$s[i]+0.5*c(-1,1)*db.mm$s.sd[i]*as.numeric(spec$a[i]>-0.5),col=palette()[db.mm$color[i]],lwd=1.5)
for(i in 1:length(db.mm$group))points(db.mm$ref.esd[i],db.mm$s[i],col=palette()[db.mm$color[i]],lwd=.5,pch=19,cex=1.5)


