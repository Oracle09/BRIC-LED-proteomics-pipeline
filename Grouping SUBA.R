temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

###to count number of localization in the each datasets and merge result together

library(dplyr)
a<-table(SUBA$Location.Consensus)
a<-data.frame(a)
setnames(a,"Freq","RootSOLUp")
SUBA<-merge(SUBA,a, by="Var1", all= TRUE)

write.csv(SUBA,"SUBA.csv")

library(data.table)
library(dplyr)
b<-`SUBACON Up space Shoot SOL.csv`%>%
  group_by(`Location.Consensus`)%>%
  summarise(Protein.Properties)

c<-split(b,f=b$Location.Consensus)

##write for a similar condition and combine the localization

d<-`SUBACON Up space shootMEM.csv`%>%
  group_by(`Location.Consensus`)%>%
  summarise(Protein.Properties)

e<-split(b,f=b$Location.Consensus)
f<-c$cytosol$Protein.Properties
g<-e$cytosol$Protein.Properties
Cytosol<-cbind(f,g)
CytosolShootUP<-Cytosol[,-1]

t<-cbind(CytosolRootDOWN,CytosolRootUP,CytosolShootUP,CytosolShootDOWN)
write.csv(t,"cystosol.csv")

#########################################################################################

b<-`SUBACON Down space Shoot SOL.csv`
group_by(`Location.Consensus`)%>%
  summarise(Protein.Properties)

c<-split(b,f=b$Location.Consensus)

##write for a similar condition and combine the localization

d<-`SUBACON Down space shoot MEM.csv`%>%
  group_by(`Location.Consensus`)%>%
  summarise(Protein.Properties)

e<-split(b,f=b$Location.Consensus)
f<-c$plastid$Protein.Properties
g<-e$plastid$Protein.Properties
Plastid<-cbind(f,g)
PlastidShootDown<-Plastid[,-1]



###########################################################
#FOR THE RNA ANALYSIS: FOLLOW THIS FOR PROTEIN ALSO.
#Download from SUBA5, choose single consensus
#Then load into R script
up<-up%>%    ##up is the data frame 
  + group_by(up$Location)%>%
  + summarise(count=n())
up



