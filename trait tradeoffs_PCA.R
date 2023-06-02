

library(FD)

##1. PCA of genome traits
bininfo<-read.csv("bin_info_traits.csv",row.names = 1,check.names = F)

bins<-bininfo
traits<-bins[,4:31]  

test<-dudi.pca(traits[,-c(3,5,21,25:27)])  ##exclude irep as it contain many NA;# exclude rare eggnog categories (almost 0) and [S] functions unknown
#select 5 axes based on the distribution of eigen values
HB<-dudi.pca(df = traits[, -c(3,5,21,25:27)], scannf = FALSE, nf = 5)
scatter(HB)
#HB$co
saveRDS(HB, file="All_gtraits_PCA.RData")
#test<-readRDS("All_gtraits_PCA.RData") #this is how to load the data


###2. PCA of community-level traits
##2.1 functional composition based on community-weighted means of genome traits
comm_bin_cov<-read.csv("bin_cov_resampled.csv",row.names = 1)  # the resampled coverage of genomes
bin.ab<-comm_bin_cov[,match(row.names(traits),colnames(comm_bin_cov))]  #match the species names
bin.ab2<-bin.ab[rowSums(bin.ab)>0,] # remove samples with no species
cwm <- functcomp(traits[, -c(3,5,21,25:27)],as.matrix(bin.ab2))  # community-weighted means of genome traits

library(FactoMineR)
bacteria.pca = PCA(cwm, scale.unit=TRUE, ncp=2, graph=T)
saveRDS(bacteria.pca, file="All_funcomp_PCA.RData")

sample.coord<-(bacteria.pca$ind)$coord  # sample coordinates in the functional space
write.csv(data.frame(sample.coord,cwm),"All CWM traits_PCAid.csv")

##2.2 functional composition based on gene annotations..
contig.traits<-read.csv("community_traits.csv",row.names = 1)

contig.traits2<-data.frame(contig.traits[,-c(5,21,25:27)]) ## exclude rare eggnog categories (almost 0) and [S] functions unknown

library(FactoMineR)
bacteria.pca2 = PCA(contig.traits2, scale.unit=TRUE, ncp=2, graph=T)
saveRDS(bacteria.pca2, file="All_funcomp_PCA_ORF_based.RData")

sample.coord<-(bacteria.pca2$ind)$coord
write.csv(data.frame(sample.coord,contig.traits2),"All CWM traits_PCAid_ORF_based.csv")




################ related figures #############

### 1. genome traits PCA
gtraits.pca<-readRDS("All_gtraits_PCA.RData")
gtraits.pca$eig/sum(gtraits.pca$eig) #explained proportion
# 0.234476414 0.155647294

#principal axes of traits
ptraits<-gtraits.pca$c1

l.x <- ptraits[,1]
l.y <- ptraits[,2]
# Draw arrows
plot(c(-0.4,0.5),c(-0.5,0.5),xlab="PCA1, 23.4%",ylab="PCA2, 15.6%")
arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="steelblue", length=0.06, lwd=1.2,angle = 18)

# Label position
labels<-row.names(ptraits)
test<-sapply(strsplit(labels,"[.][.]"), function(x) x[1])
test2<-paste0("[",test[3:21],"]")
test[3:21]<-test2
labels_simple<-test

l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")

# Variable labels
text(l.x, l.y, labels=labels_simple, pos=l.pos,cex=0.8)

###2. genome PCA scatter plot 
bin_co<-gtraits.pca$li
bininfo<-read.csv("bin_info_traits.csv",row.names = 1)
sum(row.names(bin_co)==row.names(bininfo)) #check if matched

#get the phylum and class of bins
phylum<-sapply(strsplit(bininfo$gtdb.new,";"), function(x) x[2])
class<-sapply(strsplit(bininfo$gtdb.new,";"), function(x) x[3])
phylum<-gsub("p__","",phylum)
class<-gsub("c__","",class)

group<-sapply(1:length(phylum),function(i){
  ifelse(phylum[i]=="Proteobacteria",class[i],phylum[i])
})

group.count<-c(table(group))
group.sort<-sort(group.count,decreasing = T)
top10<-names(group.sort[1:10])
group2<-sapply(1:length(group), function(i){
  ifelse(group[i] %in% top10, group[i], "other")
})

cols=c("Bacteroidota"="#e41a1c",
       "Gammaproteobacteria"="#874f6f",
       "Actinobacteriota"="#449b75",
       "Myxococcota"="#3881b0",
       "Chloroflexota"="#f17eb4",
       "Alphaproteobacteria"="#cb8cad",
       "Acidobacteriota"="#ffe528",
       "Patescibacteria"="#ffa10d",
       "Planctomycetota"="#e1c62f",
       "Verrucomicrobiota"="#b16c29",
       "other"="#999999"
)

dat<-data.frame(bin_co,group=group2)

library(ggplot2)

ggplot(dat,aes(x=Axis1,y=Axis2,colour=group))+
  geom_point(size=0.7)+
  xlab("PCA1, 23.4%")+  #get the proportion by checking the eigenvalues file
  ylab("PCA2, 15.6%")+
  scale_colour_manual(values = cols)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###3. density distribution of PC1 scores
ggplot(dat, aes(Axis1, y = stat(density), color = group, fill = group )) +
  geom_density(alpha = 0)+
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  xlab("PCA1 score")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


###4. ARG counts across MAG groups
sum(row.names(bininfo)==row.names(dat))
arg.mean=t(sapply(split(bininfo$arg_tot_counts,dat$group), function(x){
  result=c(mean(x,na.rm=T),sd(x,na.rm=T),length(!is.na(x)),sd(x,na.rm=T)/sqrt(length(!is.na(x))))
  names(result)=c("mean","sd","N","se")
  result
}))

#write.csv(arg.mean,"mean arg_taxa groups.csv")

##anova and tukey test
a1<-aov(bininfo$arg_tot_counts ~ dat$group)
summary(a1)
#Df Sum Sq Mean Sq F value Pr(>F)    
#dat$group     10  10992  1099.2   104.8 <2e-16 ***
#  Residuals   2633  27620    10.5  

posthoc <- TukeyHSD(x=a1, 'dat$group', conf.level=0.95)
posthoc$`dat$group`
write.csv(posthoc$`dat$group`,"tax group arg Tukey result.csv")

##bar plot
genome.mean<-data.frame(X=row.names(arg.mean),  arg.mean)
test<-sort(genome.mean$mean,index.return=T,decreasing = T)
genome.mean<-genome.mean[test$ix,]
genome.mean$X=factor(genome.mean$X,levels=unique(genome.mean$X))

genome.mean<-genome.mean[genome.mean$X!='other',]
library(ggplot2)
ggplot() + 
  geom_bar(data=genome.mean, aes(x=X, y=mean),stat="identity", width=0.7,fill="Goldenrod") +
  geom_errorbar(data=genome.mean, aes(x=X, y=mean, ymin=mean, ymax=mean+sd), width=.3)+
  ylim(0,22)+
  xlab(NULL)+
  ylab("ARG count")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


###5. ARG vs genome size 

ggplot(bininfo, aes(x=genome.size, y=arg_tot_counts)) + 
  geom_point(shape=21, color="white",fill="#386cb0",size=2)+
  geom_smooth(method=lm, se=FALSE, color="black")+
  xlab("Genome size")+
  ylab("ARG abundance")+
  theme_classic()

summary(lm(arg_tot_counts~genome.size,data = bininfo))
#Residual standard error: 3.25 on 2642 degrees of freedom
#Multiple R-squared:  0.2774,	Adjusted R-squared:  0.2771 
#F-statistic:  1014 on 1 and 2642 DF,  p-value: < 2.2e-16

##ARG abundance vs [Q]secondary metabolites
ggplot(bininfo, aes(x=Q..Secondary.metabolites.biosynthesis..transport.and.catabolism., y=arg_tot_counts)) + 
  geom_point(shape=21, color="white",fill="#386cb0",size=2)+
  geom_smooth(method=lm, se=FALSE, color="black")+
  xlab("Relative abundance of [Q]")+
  ylab("ARG abundance")+
  theme_classic()

summary(lm(arg_tot_counts~Q..Secondary.metabolites.biosynthesis..transport.and.catabolism.,data = bininfo))
#Residual standard error: 3.746 on 2642 degrees of freedom
#Multiple R-squared:   0.04,	Adjusted R-squared:  0.03964 
#F-statistic: 110.1 on 1 and 2642 DF,  p-value: < 2.2e-16


###5. ARG vs J
ggplot(bininfo, aes(x=J..Translation..ribosomal.structure.and.biogenesis., y=arg_tot_counts))+ 
  geom_point(shape=21, color="white",fill="#ff7f00",size=2)+
  geom_smooth(method=lm, se=FALSE, color="black")+
  xlab("Relative abundance of [J]")+
  ylab("ARG abundance")+
  theme_classic()

summary(lm(arg_tot_counts~J..Translation..ribosomal.structure.and.biogenesis.,data = bininfo))
#Residual standard error: 3.343 on 2642 degrees of freedom
#Multiple R-squared:  0.2351,	Adjusted R-squared:  0.2348 
#F-statistic: 812.1 on 1 and 2642 DF,  p-value: < 2.2e-16

###6. ARG vs iRep
ggplot(bininfo, aes(x=irep, y=arg_tot_counts))+ 
  geom_point(shape=21, color="white",fill="#ff7f00",size=2)+
  geom_smooth(method=lm, se=FALSE, color="black")+
  xlab("MAG iRep")+
  ylab("ARG abundance")+
  theme_classic()

summary(lm(arg_tot_counts~irep,data = bininfo))
#Residual standard error: 3.873 on 2078 degrees of freedom
#(564 observations deleted due to missingness)
#Multiple R-squared:  0.01452,	Adjusted R-squared:  0.01404 
#F-statistic: 30.61 on 1 and 2078 DF,  p-value: 3.553e-08


###community-level traits PCA plots

###1. CWM traits PCA (MAG-based)


cwm_mag.pca<-readRDS("All_funcomp_PCA.RData")
cwm_mag.pca$eig  #explained proportion
# 26.2670909 19.7618067
#principal axes of traits
ptraits<-(cwm_mag.pca$var)$coord

l.x <- ptraits[,1]
l.y <- ptraits[,2]
# Draw arrows
plot(c(-1.0,1.0),c(-1.0,1.0),xlab="PCA1, 26.3%",ylab="PCA2, 19.8%")
arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="steelblue", length=0.06, lwd=1.2,angle = 18)

# Label position
labels<-row.names(ptraits)
test<-sapply(strsplit(labels,"[.][.]"), function(x) x[1])
test2<-paste0("[",test[3:21],"]")
test[3:21]<-test2
labels_simple<-test

l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")

# Variable labels
text(l.x, l.y, labels=labels_simple, pos=l.pos,cex=0.8)

###2. CWM PCA scatter plot (MAG-based)
env<-read.csv("GWMC_env.csv",row.names = 1)
sample_co<-(cwm_mag.pca$ind)$coord
env2<-env[match(row.names(sample_co),row.names(env)),]   #match sample names
sum(row.names(sample_co)==row.names(env2)) #check if matched

dat<-data.frame(sample_co,group=env2$Continent)

library(ggplot2)

ggplot(dat,aes(x=Dim.1,y=Dim.2,colour=group))+
  geom_point(size=1)+
  xlab("PCA1, 26.3%")+  
  ylab("PCA2, 19.8%")+
  scale_colour_manual(values=c("#4DBBD5FF","#E64B35FF","#7E6148FF","#00A087FF","#3C5488FF","#F39B7FFF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###4. Community traits PCA (gene annotation-based)
cwm_read.pca<-readRDS("All_funcomp_PCA_ORF_based.RData")
cwm_read.pca$eig  #explained proportion
# 23.3070095  18.2290359
#principal axes of traits
ptraits<-(cwm_read.pca$var)$coord

l.x <- ptraits[,1]
l.y <- ptraits[,2]
# Draw arrows
plot(c(-1.0,1.0),c(-1.0,1.0),xlab="PCA1, 23.3%",ylab="PCA2, 18.2%")
arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="steelblue", length=0.06, lwd=1.2,angle = 18)

# Label position
labels<-row.names(ptraits)
test<-labels
test2<-paste0("[",test[4:22],"]")
test[4:22]<-test2
labels_simple<-test

l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")

# Variable labels
text(l.x, l.y, labels=labels_simple, pos=l.pos,cex=0.8)

###5. CWM PCA scatter plot (read-based)
sample_co<-(cwm_read.pca$ind)$coord
sum(row.names(sample_co)==row.names(env)) #check if matched


dat<-data.frame(sample_co,group=env$Continent)

ggplot(dat,aes(x=Dim.1,y=Dim.2,colour=group))+
  geom_point(size=1)+
  xlab("PCA1, 23.3%")+  
  ylab("PCA2, 18.2%")+
  scale_colour_manual(values=c("#4DBBD5FF","#E64B35FF","#7E6148FF","#00A087FF","#3C5488FF","#F39B7FFF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

