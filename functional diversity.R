
### 0. Prepare the 'species abundance table'
bincov<-read.csv("bin_coverage.csv",row.names = 1) #the raw coverage data for genomes
bincov.rnd<-round(bincov)

library(vegan) 
test<-rowSums(bincov.rnd)
hist(test)
sort(test)[22:23]

#remove 10% samples with small coverages
comm<-bincov.rnd[test>345,]
res.dep = min(rowSums(comm)) #resample size (min total coverage)
comm.res = rrarefy(comm, sample = res.dep)  #resample
comm.res = comm.res[, colSums(comm.res) > 0] #remove genomes with 0 abundance
dim(comm.res);sum(comm.res[1,]);sum(comm.res[11,]); #check

env2<-env[match(row.names(comm.res),row.names(env)),]
sum(row.names(env2)==row.names(comm.res))
rch<-rowSums(comm.res>0)  ## richness of MAGs
plot(env2$HQ.bases, rch)
summary(lm(rch~env2$HQ.bases))
#the sequence depth now only explain 1.3% of the variations of MAG richness. p=0.056 not significant
write.csv(comm.res,"bin_cov_resampled.csv")



### 1. functional diversity

library(FD)

bininfo<-read.csv("bin_info_traits.csv",row.names = 1,check.names = F)
comm_bin_cov<-read.csv("bin_cov_resampled.csv",row.names = 1)
sum(row.names(bininfo)==colnames(comm_bin_cov)) # check if the species matched

traits<-bininfo[,4:31]   
bin.ab2<-comm_bin_cov[rowSums(comm_bin_cov)>0,]
HB_FD<-dbFD(x=traits[, -c(3,5,21,25:27)], a=bin.ab2,m=5) ##exclude irep as it contain many NA;# exclude rare eggnog categories (almost 0) and [S] functions unknown
#select 5 axes based on the distribution of eigen values
saveRDS(HB_FD, file="All_FD.RData")  ##save the functional diversity data

#################### related figures ############
###1. functional diversity box plot
FD<-readRDS("All_FD.RData")
dat<-data.frame(FRic=FD$FRic,FDis=FD$FDis)

env<-read.csv("GWMC_env.csv",row.names = 1)
env2<-env[match(row.names(dat),row.names(env)),]   #match sample names
sum(row.names(env2)==row.names(dat))   

library(ggplot2)
dat=data.frame(dat,Continent=env2$Continent)
ggplot(dat,aes(x=Continent,y=FRic,colour=Continent))+
  geom_boxplot(outlier.alpha =0)+
  geom_point(position=position_jitterdodge(0.7),alpha=0.6)+
  scale_colour_manual(values=c("#4DBBD5FF","#E64B35FF","#7E6148FF","#00A087FF","#3C5488FF","#F39B7FFF"))+
  xlab(NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

a1<-aov(dat$FRic ~ dat$Continent)
summary(a1)
#Df   Sum Sq Mean Sq F value Pr(>F)
#dat$Continent   5   880788  176158   1.292  0.269
#Residuals     198 26987198  136299 
posthoc <- TukeyHSD(x=a1, 'dat$Continent', conf.level=0.95)
posthoc

ggplot(dat,aes(x=Continent,y=FDis,colour=Continent))+
  geom_boxplot(outlier.alpha =0)+
  geom_point(position=position_jitterdodge(0.7),alpha=0.6)+
  scale_colour_manual(values=c("#4DBBD5FF","#E64B35FF","#7E6148FF","#00A087FF","#3C5488FF","#F39B7FFF"))+
  xlab(NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

a1<-aov(dat$FDis ~ dat$Continent)
summary(a1)
#Df Sum Sq Mean Sq F value Pr(>F)  
#dat$Continent   5  1.842  0.3683   2.958 0.0134 *
#  Residuals     198 24.657  0.1245 
posthoc <- TukeyHSD(x=a1, 'dat$Continent', conf.level=0.95)
posthoc
#South America-Europe        -0.381621466 -0.7069666 -0.056276379 0.0112877
#South America-North America -0.303557320 -0.5602088 -0.046905812 0.0102983

###2. functional diversity vs absolute latitude
fundiv<-readRDS("All_FD.RData")
fundiv2<-data.frame(FRic=fundiv$FRic,FDis=fundiv$FDis)
env<-read.csv("GWMC_env.csv",row.names = 1)
env2<-env[match(row.names(fundiv2),row.names(env)),]   #match sample names of fundiv
sum(row.names(fundiv2)==row.names(env2))  #check if matched

library(ggplot2)
library(RColorBrewer)

dat<-data.frame(fundiv2,env2) 
hemi<-ifelse(dat$Latitude.Google>0,"North","South")
dat2<-data.frame(dat,hemi=hemi)
dat2$temp=dat2$Air.temperature.Mean.annual

ggplot(dat2,aes(x=lat.abs,y=FRic))+
  geom_point(aes(shape=hemi,fill=temp))+
  #geom_smooth(method = "lm",formula =y ~ poly(x, 2),col="gray20")+
  scale_shape_manual(values = c(21,22))+
  scale_fill_gradientn(colours=rev(brewer.pal(9,"RdBu")))+
  xlab("Absolute latitude") +
  ylab("Functional richness") +
  theme_bw() +
  #theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(dat2,aes(x=lat.abs,y=FDis))+
  geom_point(aes(shape=hemi,fill=temp))+
  #geom_smooth(method = "lm",formula =y ~ poly(x, 2),col="gray20")+
  scale_shape_manual(values = c(21,22))+
  scale_fill_gradientn(colours=rev(brewer.pal(9,"RdBu")))+
  xlab("Absolute latitude") +
  ylab("Functional dispersion") +
  theme_bw() +
  #theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
