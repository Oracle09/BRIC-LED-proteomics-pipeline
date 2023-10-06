
### 7/26/2022
### create a list of TAIR ID, Description, SYMBOLS AND LFC VALUE
###Make tables for publication

#change name of enrich column
colnames(enrich)[1]<- "TAIR..ID"
a<-merge(RootMD.csv,enrich, by="TAIR..ID")

#REMOVE COLUMNS
a<-a[,-(2:3)]
a<-a[,-(3:12)]
a<-a[,-(6:17)]
a<-a[,-(7:13)]
a<-a[,-(4)]


#merge with shoot data
colnames(Root)[3]<-"TAIR..ID"
b<-merge(a,ShootMD.csv,by="TAIR..ID")
b<-b[,-(6:17)]
b<-b[,-(10:21)]
b<-b[,-(6)]
b<-b[,-(7)]
write.csv(b,"Rootdownshootup.csv")


###Make Table for Root and Shoot Photosynthetic proteins
library(dplyr)
library(tidyr)
a<-RootM.csv%>%
  filter(grepl('pet',Description,ignore.case = TRUE))
b<-ShootM.csv%>%
  filter(grepl('pet',Description,ignore.case = TRUE))
a<-a[,-(5:14)]
a<-a[,-(8:19)]
a<-a[,-(1:2)]
a<-a[,-(4)]

b<-b[,-(5:14)]
b<-b[,-(8:19)]
b<-b[,-(1:2)]
b<-b[,-(4)]

##merge both
c<-left_join(a,b, by="Description")
#eliminate TAIR ID
c<-c[,-5]

write.csv(c,"Transport.csv", row.names = FALSE)

###For transport proteins

library(data.table)

c<-merge(ShootM.csv,shootT.csv,by="TAIR.ID")
c<-c[,-(5:14)]
c<-c[,-(6:19)]
c<-c[,-(2:3)]


shootR<-c

c<-merge(shootR,rootR, by="TAIR.ID", all=TRUE)
## coalesce the two Description
c$Description.x<-ifelse(is.na(c$Description.x),c$Description.y, c$Description.x)
c<-c[,-4]


###For transport proteins graph
a<-rootUP[,-(1:3)]
a<-a[,-(3:4)]

b<-rootD[,-(1:3)]
b<-b[,-(3:4)]
d<-merge(a,b, by= "Pathway", all=TRUE)
e<-merge(c,d, by = "Pathway", all=TRUE)
write.csv(e, "transportgraph.csv")

###Classify transport proteins into class
###split columns in a data frame as separate data frame
library(data.table)
c<-split.default(a,names(a))
cation<-c$Cation.transport

setnames(c$Cation.transport, "c$Cation.transport" ,"x")

c$Protein.localization<-setdiff(c$Protein.localization,c$Macromolecule.localization)
CationT<-merge(b,c$Cation.transport, by = "x", all.y = TRUE)

c$Cation.transport<-c$Cation.transport[-(18:294),]
c$Cation.transport<-as.data.frame(c$Cation.transport)
setnames(c$Cation.transport, "c$Cation.transport" ,"x")

e<-merge(b$cation,by= "x" , all=T)

cation$y<-c(1:17)


###CELL WALL merge
a<-left_join(Cellwall,rootM, by="TAIR.ID")

b<-left_join(Cellwall,shootS,by="TAIR.ID")
b<-b[,-(2:3)]
b<-b[,-(3:12)]
b<-b[,-(6:17)]
b<-b[,-(7:13)]
b<-b[,-(4)]
b<-b[,-4]
e<-merge(b,c, by="TAIR.ID")
e$Description.x<-ifelse(is.na(e$Description.x),e$Description.y, e$Description.x)


###Remove duplicate rows in e
library(dplyr)
f<-e%>%distinct(e$TAIR.ID, .keep_all=TRUE)
g<-left_join(a,f, by="TAIR.ID", keep_all=TRUE)  ###joing shoot and root together
write.csv(g, "Cell wall protein.csv")


###for redox

#split redox columns
a<-split.default(redox,names(redox))
rD<-a$ROOT.DOWN
rU<-a$ROOT.UP
sD<-a$SHOOT.DOWN
sU<-a$SHOOT.UP
setnames(sD, "SHOOT.DOWN" ,"TAIR.ID")   ##repeat for all

RedoxSDS<-left_join(sD,shootS, by="TAIR.ID")

RedoxSD<-RedoxSD[,-(4:31)]


RedoxSDS<-RedoxSDS[,-(2:3)]
RedoxSDS<-RedoxSDS[,-(3:12)]
RedoxSDS<-RedoxSDS[,-(6:17)]
RedoxSDS<-RedoxSDS[,-(7:13)]
RedoxSDS<-RedoxSDS[,-(4)]
RedoxSDS<-RedoxSDS[,-4]

A<-merge(RedoxSU, RedoxSUS, by="TAIR.ID", all = TRUE)  ##merge soluble and membrane for RootUP
A$Description.x.x<-ifelse(is.na(A$Description.x.x),A$Description.x.y, A$Description.x.x)
A<-A[,-(4)]



RedoxRootUP<-A
RedoxRD
RedoxShootUP<-A
RedoxShootDOWN<-A

###merge all together to form a comprehensive list
A<-merge(RedoxShootUP,RedoxShootDOWN, by="TAIR.ID", all = TRUE)
ROOT<-A
SHOOT<-A
B<-merge(ROOT,SHOOT, by="TAIR.ID", all=TRUE)

write.csv(B, "Redox_full.csv")
