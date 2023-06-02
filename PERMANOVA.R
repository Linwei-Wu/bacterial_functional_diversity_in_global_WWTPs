

##0. read data of functional composition
fcomp_bin<-read.csv("All CWM traits_PCAid.csv",row.names = 1)   ##community-level traits based on CWM genomic traits
fcomp_orf<-read.csv("All CWM traits_PCAid_ORF_based.csv",row.names = 1)  ##community-level traits based on direct gene annotations
env<-read.csv("GWMC_env.csv",row.names = 1) # environmental data


##1. calculate beta distances
library(vegan)  #install.packages("vegan")
comm1<-scale(fcomp_bin[,-(1:2)],center =F)    #cwm functional traits from bin data
comm2<-scale(fcomp_orf[,-(1:2)],center =F)    #read-based functional groups

beta.weight1=vegdist(comm1)
beta.weight2=vegdist(comm2)


##2.PERMANOVA (adonis) on multi-variables
#e.g., adonis of continent, climate type and activated sludge

#fcomp based on bin data
bray<-as.matrix(beta.weight1)

env1<-env[match(row.names(bray),row.names(env)),] #match
sum(colnames(bray)==row.names(env1))   #check
grp=env1[,c(2,9,11)] #the three categorical variables
head(grp)

id<-complete.cases(grp) #indicate which samples have no missing values

adonis.bray<-adonis2(as.dist(bray[id,id])~Continent*Climate.type2*activated.sludge.type,data=grp[id,])

sink(file="adonis.Bray.fcomBin.summary.txt")
adonis.bray
sink()

#fcomp based on direct gene annotation data

bray<-as.matrix(beta.weight2)

env1<-env[match(row.names(bray),row.names(env)),] #match
sum(colnames(bray)==row.names(env1))
grp=env1[,c(2,9,11)] #the three categorical variables
head(grp)

id<-complete.cases(grp) #indicate which samples have no missing values

adonis.bray<-adonis2(as.dist(bray[id,id])~Continent*Climate.type2*activated.sludge.type,data=grp[id,])

sink(file="adonis.Bray.fcomORF.summary.txt")
adonis.bray
sink()


