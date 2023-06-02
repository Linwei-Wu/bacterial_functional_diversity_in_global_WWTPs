
###prepare the data
env<-read.csv("GWMC_env.csv",row.names = 1)
fundiv<-readRDS("All_FD.RData")
funcomp<-read.csv("All CWM traits_PCAid_ORF_based.csv",row.names = 1)
fundiv2<-data.frame(FRic=fundiv$FRic,FDis=fundiv$FDis)

sum(row.names(env)==row.names(funcomp))
fundiv3<-fundiv2[match(row.names(env),row.names(fundiv2)),]
sum(row.names(env)==row.names(fundiv3))

dat<-data.frame(Fun.PC1=funcomp$Dim.1,fundiv3,ARG.abundance=funcomp$arg_tot,env[,16:86])

###SEM
library(lavaan)
dat.scaled<-data.frame(scale(dat,center = F))


model <- '
# regressions

Fun.PC1~Influent.BOD..mg.L.+F.M
F.M~Influent.BOD..mg.L.+population
FDis~F.M+pH
pH~F.M
FRic~population+pH
ARG.abundance~Fun.PC1+population
BOD.removal.percentage~Fun.PC1+Influent.BOD..mg.L.+population+F.M
NH4.Nremoval.percentage~Fun.PC1+FRic+FDis+F.M
FDis~~FRic
Fun.PC1~~FRic
Fun.PC1~~FDis
'
Fit <- sem(model, missing="ml",data=dat.scaled,fixed.x =F)

summary(Fit, rsquare=T, standardized=T,fit.measures=TRUE)

# residuals(Fit, type="cor")
# modificationIndices(Fit,standardized=F)

 
 sink(file="SEM.txt")
 summary(Fit, rsquare=T, standardized=T,fit.measures=TRUE)
 sink()
 

 
 