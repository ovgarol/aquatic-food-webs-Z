porta2019 = function(name){
  limit1 = read.csv(paste0("./physio_limits/",name), header=FALSE, comment.char="#")
  kk = order(limit1$V1)
  limit1 = limit1[kk,] 
  limit1$V1=2*(limit1$V1*2e-3/1e3/(4*pi/3))**(1./3.)*1e6
  limit1$V2=2*(limit1$V2*2e-3/1e3/(4*pi/3))**(1./3.)*1e6
  return(limit1)
}

####
## Figure A3
####
par(mfrow=c(1,1),family='Helvetica')
plot(1,xlim=range(db.all$esd),ylim=range(db.all$opt),log='xy',type='n',
     xaxt='n',yaxt='n',
     xlab='predator size (m)',ylab='optimal prey size (m)',
     main='allometric rule and data')
the.col=as.numeric(as.factor(db.all$group))
the.shape = c(25,23,24)[1 + as.numeric(db.all$ss>-0.33) + as.numeric(db.all$ss >0.33)]

physio.col= c('lightgray','cyan','gray','darkmagenta','magenta','white')
if(T)for(i in 2:6){
  limit1 = porta2019(paste0('limit',i))
  #lines(limit1,lwd=2)
  s<-lm(log10(V2)~log10(V1),data=limit1)
  #abline(s,lwd=3)
  x0 = c(1e-1,1e8)
  x1 = 10^(predict(s,newdata = data.frame(V1=x0)))
  polygon(c(x0,rev(x0)),c(x1,c(0,0)+1e-7),col=physio.col[i],border=NA)
  s<-lm(log(V2)~log(V1),data=limit1)
  
  print(s$coefficients)
}

lines(size.range,size.range,lty=2,col='black')
v=1:(17**2) * 0.5
my.lines(exp(v),exp(OPS(v,0,0,0,a=0)),col='black',lwd=3)
abline(a=r.offset-0.33,b=1,col='gold',lwd=3)

#points(opt~esd,data=db.all,pch=the.shape,bg=palette()[the.col],lwd=1,cex=2*(1+0*as.numeric(the.shape>23)),col='white')
#points(opt~esd,data=db.all,pch=1,bg='gray',lwd=1,cex=0.1,col='black')

#polygon(1e6*c(10e-5,10e-5,1e-3,1e-3),1e4*c(1e-4,1e-1,1e-1,1e-4),lwd=3,angle=45,fillOddEven=T)

y=c(1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7)
y.lab = c(' ',TeX('10$^{-6}$'),' ',TeX('10$^{-4}$'),' ',TeX('10$^{-2}$'),' ',TeX('10$^{0}$'),' ')
axis(side=1, at=(y), labels=y.lab,las=1,cex.lab=1.)
axis(side=2, at=(y), labels=y.lab,las=1,cex.lab=1.)

## to plot data
if(F){
  db = read.csv("./external_data/283_2_FoodWebDataBase_2018_12_10.csv")
  db = db[!is.na(db$autoID),]
  db = db[db$interaction.type%in%c('predacious','herbivorous'),]
  #db = db[db$interaction.classification=='nibi',]
  db = db[db$ecosystem.type %in% c('marine','lakes','streams'),]
  db$type = paste(db$con.movement.type,db$con.metabolic.type,db$interaction.type)
  db$cols = 1+as.numeric(db$interaction.classification=='nibi')
  db = db[db$con.movement.type=='swimming',]
  
  esds = 2*(db$con.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  opss = 2*(db$res.mass.mean.g.*1e-3/1e3/(4*pi/3))**(1./3.)*1e6
  points(esds,opss,cex=0.01+0.2*db$cols,pch=19,col=db$cols)  

}


ii = 1
w = 3.4

## to plot model
if(F)for(i in 1:4){
  w =3.5# 0.5*exp(1)*(3.5-.8975*log(abs(0.5*(ii+w)-9)))
  for(j in -2:1){
    the.s = j*(-0.3 + 0.2*(ii+w+1*(j < -1)+abs(j))) #j*(ii+0.5*w)/5
    if(ii>log(1e4)) the.s = sqrt(4)*j
    if(j< -1 & (ii+w > log(1e4) )) next
    minimal.model(exp(ii+abs(j)+0*(j < -1)),exp(ii+w+abs(j)+1*(j < -1)),the.s,col=c('darkmagenta','darkmagenta','darkgray','darkcyan')[3+j])
    print(c(exp(ii+abs(j)+0*(j < -1)),exp(ii+w+abs(j)+1*(j < -1)),the.s))
  }
  
  ii = ii+w
}




