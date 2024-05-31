library(scales)
library(stringr)

db.zoo = read.csv("./data/Hansen1994.csv")
db.zoo$type = db.zoo$group
db.zoo$group = 'Plankton'
db.zoo = subset(db.zoo,select = c(specie,esd,opt,type,group))
db.zoo$source = 'Hansen 1994'

###

db.dino = read.csv("./data/Garcia-Oliva2022a.csv")
db.dino$type = 'Dinoflagellates'
db.dino$type2 = as.numeric(db.dino$Specialization.factor..s)
db.dino$specie = db.dino$Species.reported.name
db.dino$opt = as.numeric(db.dino$Optimal.prey.size...OPS..um.)
db.dino$esd = as.numeric(db.dino$Equivalent.spherical.diameter..ESD..um.)
db.dino$group = 'Plankton'
db.dino$col = as.numeric(as.factor(db.dino$type)) 
db.dino = subset(db.dino,select = c(specie,esd,opt,type,group))
db.dino$source = 'Garcia-Oliva 2022'

###

db = read.delim("./data/Barnes2008.txt",stringsAsFactors=F,na.strings="n/a")

spec = NULL
spec$specie = levels(as.factor(db$Predator))
spec$type = levels(as.factor(db$Predator))
spec$esd = 0 * 1:length(spec$specie)
spec$opt = 0 * 1:length(spec$specie)
spec$n.values = 0 * 1:length(spec$specie)

i = 1
for(j in spec$specie){
  temp = db[db$Predator==j,]
  spec$type[i] = temp$Type.of.feeding.interaction[1]
  spec$esd[i] = 2*((median(temp$Predator.mass)*1e-3/1e3/(4*pi/3))**(1./3.))*1e6 # in micrometers
  spec$opt[i] = 2*((median(temp$Prey.mass)*1e-3/1e3/(4*pi/3))**(1./3.))*1e6 # in micrometers
  spec$n.values[i] = length(temp$Predator.mass)
  i = i+1
  if(T)if(j=='Clupea harengus'){
    print(temp$Type.of.feeding.interaction)
    break
  } 
}

spec = as.data.frame(spec)
spec = spec[spec$n.values>5,]
spec$group = 'Fish' 
db.fish = subset(spec,select=c(specie,esd,opt,type,group))
db.fish$source = 'Barnes 2008'

###

db.mammal = read.csv("./data/Tucker2014.csv")

db.mammal = db.mammal[db.mammal$Environment == 'Marine',]
db.mammal = db.mammal[db.mammal$Diet == 'Carnivore',]
db.mammal$esd = 2*(10^db.mammal$Mass..log10.kg.*1e-3*(3/4/pi))^(1/3) * 1e6 # in micrometers
db.mammal$opt = 2*(10^as.numeric(db.mammal$Prey.mass..log10.kg.)*1e-3*(3/4/pi))^(1/3) *1e6 # in micrometers
db.mammal$specie = as.factor(gsub("_", " ", db.mammal$Taxon))
db.mammal$type = db.mammal$Niche
db.mammal$group = 'Mammals'

db.mammal = subset(db.mammal,select = c(specie,esd,opt,type,group))
db.mammal$source = 'Tucker 2014'

###

db = read.csv("./data/Brose2019.csv")
db = db[db$con.taxonomy.level=='species',]
db = db[db$res.taxonomy.level=='species',]
db = db[!is.na(db$autoID),]
db = db[db$interaction.type=='predacious',]
db = db[db$ecosystem.type %in% c('marine','lakes','streams'),]
## to restrict the type of movement
#db = db[db$con.movement %in%c('swimming','other'),]

db = db[db$con.metabolic.type=='invertebrate',]
db$type = paste(db$con.movement.type,db$interaction.type)

db = db[db$type!='NA NA NA',]
unique(db$interaction.type)

spec = NULL
spec$specie = levels(as.factor(db$con.taxonomy))
spec$type = levels(as.factor(db$con.taxonomy))
spec$n.values = 0 * 1:length(spec$Predator)
spec$esd = 0 * 1:length(spec$Predator)
spec$opt = 0 * 1:length(spec$Predator)

i = 1
for(j in spec$specie){
  temp = db[db$con.taxonomy==j,]
  spec$type[i] = temp$type[1]
  spec$esd[i] = 2*(median(temp$con.mass.mean.g.)*1e-3/1e3/(4*pi/3))**(1./3.)*1e6 # in micrometers
  spec$opt[i] = 2*(median(temp$res.mass.mean.g.)*1e-3/1e3/(4*pi/3))**(1./3.)*1e6 # in micrometers
  spec$n.values[i] = length(temp$con.mass.mean.g.)
  i = i+1
}

spec$group = 'Invertebrate'
spec = as.data.frame(spec)
spec = spec[spec$n.values>5,]

db.inv = subset(spec,select=c(specie,esd,opt,type,group))
db.inv$source = 'Brose 2019'

levels(as.factor(db.inv$type))
write.csv(spec,'many-taxa.csv')

###

db.otr = read.csv("./data/Fuchs2010.csv",comment.char="#")
db.otr = as.data.frame(db.otr)

db.otr$esd = db.otr$multiplier*db.otr$size*1e6
db.otr$opt = db.otr$multiplier*db.otr$opt*1e6
db.otr$type = db.otr$group
db.otr$group = db.otr$group.2
db.otr = subset(db.otr,select=c(specie,esd,opt,type,group))
db.otr$source = 'Fuchs 2010'

###

rm(db,spec,temp,i,j)

###

db.all = rbind(db.zoo,db.dino,db.otr,db.inv,db.fish,db.mammal)
db.all = unique(db.all)
db.all = db.all[complete.cases(db.all),]
db.all$name = db.all$specie
db.all$specie = word(paste(db.all$specie, 'sp.'), 1,2, sep=" ")
unique(db.all$group)

###

db.names = read.csv("./data/WORMS_names_taxonomy.csv",comment.char="#")
db.names = subset(db.names,select=c(name,kingdom,phylum,class,order,family,genus))

db.final = merge(db.all,db.names,by='name',all.x=T,all.y=F)

db.all=db.final

###
### my classification
if(T){
  db.all$group[db.all$group=='inv'] ='Invertebrate'
  db.all$type[db.all$type%in%c('Mammal-fish')] = 'Mammal'
  db.all$type[db.all$type%in%c('Squid-fish')] = 'Fish'
  db.all$type[db.all$type%in%c('Doliolid')] = 'Salp'
  db.all$group[db.all$type%in%c('swimming ectotherm vertebrate')] = 'Zz'
  db.all$group[db.all$type%in%c('swimming endotherm vertebrate')] = 'Zz'
  db.all$type[db.all$type%in%c('swimming ectotherm vertebrate')] = 'EctoV'
  db.all$type[db.all$type%in%c('swimming endotherm vertebrate')] = 'EndoV'
  db.all$type[db.all$type%in%c('predacious/piscivorous')] = 'predacious'
  db.all$type[db.all$type%in%c('other')] = 'predacious'
  db.all$type[db.all$type%in%c('Dinoflagellate')] = 'dinos'

  db.all$type[grep('Daphnia',db.all$specie)] = 'Cladoceran'
  db.all$group[db.all$type%in%c('Cladoceran','Rotifer')] = 'Invertebrate'
  db.all$type[db.all$type%in%c('Cladoceran','Rotifer')] = 'sessile herbivorous'
  
  ### CHANGE COPEPODS
  if(F)for(i in c('Copepod','Meroplankton','Chaetognath')){
    db.all$group[db.all$type==i] = 'Invertebrate'
    db.all$type[db.all$type==i] = 'swimming predacious'
  }
  
  db.final = db.all
}


###
### taxonomical classification
if(T){

  #ciliates
  db.final$group[db.final$phylum=='Ciliophora'] = 'Plankton'
  db.final$type[db.final$phylum== 'Ciliophora'] = 'Ciliates'
  
  #copepods
  db.final$group[db.final$type=='Copepod'] = 'Invertebrate'#'Plankton'
  db.final$type[db.final$type=='Copepod'] = 'swimming predacious'#'Copepoda'
  
  db.final$group[db.final$class=='Copepoda'] = 'Invertebrate'#'Plankton'
  db.final$type[db.final$class=='Copepoda'] = 'swimming predacious'#'Copepoda'
  
  #Meroplankton
  db.final$group[db.final$type=='Meroplankton'] = 'Invertebrate'#'Plankton'
  db.final$type[db.final$type=='Meroplankton'] ='Meroplankton'
  
  #mammals
  db.final$group[db.final$class=='Mammalia'] = 'Mammals'
  
  #chaetognath
  db.final$group[db.final$phylum=='Chaetognatha'] = 'Invertebrate'
  db.final$type[db.final$phylum== 'Chaetognatha'] = 'swimming predacious'
  
  #fish
  db.final$group[db.final$class=='Teleostei'] = 'Fish'
  db.final$group[db.final$class=='Elasmobranchii'] = 'Fish'
  
  #dinoflagellates
  db.final$group[db.final$class=='Dinophyceae'] = 'Plankton'
  db.final$type[db.final$class=='Dinophyceae'] = 'Dinoflagellates'
  
  #ctenophora
  db.final$group[db.final$phylum=='Ctenophora'] = 'Jellyfish'
  db.final$type[db.final$phylum== 'Ctenophora'] = 'Ctenophora'
  
  #Scyphozoa
  db.final$group[db.final$class=='Scyphozoa'] = 'Jellyfish'
  db.final$type[db.final$class=='Scyphozoa'] = 'Scyphozoa'
  
  #Siphonophorae
  db.final$group[db.final$order=='Siphonophorae'] = 'Jellyfish'
  db.final$type[db.final$order=='Siphonophorae'] = 'Siphonophorae'
  
  db.final$group[db.final$type=='Siphonophora'] = 'Jellyfish'
  db.final$type[db.final$type=='Siphonophora'] = 'Siphonophorae'
  
  #Anthozoa = anemoni
  db.final$group[db.final$class=='Anthozoa'] = 'Invertebrate'
  db.final$type[db.final$class=='Anthozoa'] = 'sessile predacious'
  
  #Tubularia = anemoni-like
  db.final$group[db.final$family=='Tubulariidae'] = 'Invertebrate'
  db.final$type[db.final$family=='Tubulariidae'] = 'sessile predacious'
  
  #Branchiopoda = cladoceran-likes
  #db.final$group[db.final$class=='Branchiopoda'] = 'Invertebrate'
  #db.final$type[db.final$class=='Branchiopoda'] = 'sessile herbivorous'
  
  #cladoceras = cladoceran-likes
  db.final$group[db.final$type=='Cladoceran'] = 'Invertebrate'
  db.final$type[db.final$type=='Cladoceran'] = 'sessile herbivorous'
  
  #Rotifera
  db.final$group[db.final$phylum=='Rotifera'] = 'Invertebrate'
  db.final$type[db.final$phylum=='Rotifera'] = 'sessile herbivorous'
  
  db.final$group[db.final$type=='Rotifer'] = 'Invertebrate'
  db.final$type[db.final$type=='Rotifer'] = 'sessile herbivorous'
  
  #Thaliacea = salps and dolioids
  db.final$group[db.final$class=='Thaliacea'] = 'Jellyfish'
  db.final$type[db.final$class=='Thaliacea'] = 'Salp'
  
  #Flagellate
  db.final$group[db.final$type=='Flagellate'] = 'Plankton'
  db.final$type[db.final$type=='Flagellate'] = 'Flagellate'
  
  ##default to invertebrate
  db.final$group[db.final$type=='sessile herbivorous'] = 'Invertebrate'
  
  #db.final$group[is.na(db.final$group)] = 'Invertebrate'
  db.final$type[db.final$group=='Invertebrate' & db.final$type=='predacious'] = 'swimming predacious'
  
}

###
###

db.all = unique(db.final)
db.all$tag = as.factor(paste(db.all$group,'-',db.all$type))

unique(db.all$group)
unique(db.all$type)
unique(db.all$tag)

db.all = db.all[!db.all$type%in%c('walking predacious','sessile predacious'),]
#db.all = db.all[db.all$specie != 'Noctiluca scintillans',]

db.all$type[db.all$type=='Scyphozoa'] = 'Siphonophorae'
db.all$type[db.all$type=='Ctenophora'] = 'Siphonophorae'
db.all$type[db.all$type=='Squid'] = 'Fish'
db.all$type[db.all$type=='Ciliate'] = 'Ciliates'
db.all$type[db.all$type=='Siphonophorae'] = 'tentacled'

db.all$tag = as.factor(paste(db.all$group,'-',db.all$type))
db.all$X = 1:length(db.all$esd)
write.csv(db.all,file='db_all_dirty2.test.csv',row.names = F)
