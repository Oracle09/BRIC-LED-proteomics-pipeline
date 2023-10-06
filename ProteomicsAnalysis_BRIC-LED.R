#Script origin
#Author: Colin Kruse
#Modified by: Gbolaga Olanrewaju
#Collaborator: Gbolaga (Wyatt Lab)
#Description: Import and initial annotation of protein data from BRIC-LED studies
#Anticipation work: all sheets from the reports provided by Nebraska collaborators were saved as .csv files

library(data.table)


# Import ------------------------------------------------------------------

setwd("C:/Users/Gbolaga/Desktop/BRIC LED RESULTS/RInputs/")


#Import all reports 
#Best import for multiple csv import:
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))


#Colin's import script:
for(reports in list.files(pattern = ".csv")){
  assign(strsplit(reports,split = "\\.")[[1]][1],fread(reports,skip = 1))}


#Import accession IDs from UniProt files (https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/)
IDs<-fread("ARATH_3702_idmapping.dat")
table(IDs$`UniProtKB-ID`)

RelevantAccessions<-c(RootMEM.csv$Accession,RootSOL.csv$Accession,ShootMEM.csv$Accession,ShootSOL.csv$Accession)  
RelevantAccessions<-unique(RelevantAccessions)

#Link UniprotID to ensemble genome and co
EnsemblGenomeIDs<-IDs[IDs$`UniProtKB-ID` == "EnsemblGenome",]
AraportIDs<-IDs[IDs$`UniProtKB-ID` == "Araport",]
ExtraIDs<-IDs[IDs$`UniProtKB-ID` == "Gene_OrderedLocusName",]
ManualAdd<-as.data.frame(c("A0A2P2CLH9","A0A2P2CLF6","A0A2P2CLF9","A0A2P2CLF3"))
colnames(ManualAdd)<-"P48347"   #change column name
#ManualAdd$P48347<-c("A0A2P2CLH9","A0A2P2CLF6","A0A2P2CLF9","A0A2P2CLF3")
ManualAdd$'UniProtKB-ID'<-"ManualAddition"
ManualAdd$'14310_ARATH'<-c("ATMG00640", "ATMG00070", "ATCG00120", "ATMG00090")

TairIDs<-rbind(EnsemblGenomeIDs,AraportIDs,ExtraIDs,ManualAdd)

TairIDs<-TairIDs[!duplicated(TairIDs$P48347),]   #remove duplicated value
  
#Checking for fully matched set of accessions 
summary(RelevantAccessions %in% TairIDs$P48347)

#The Stuff below was used to identify the additional fields 
#I needed to keep and the what was missing altogether ID wise
absence<-RelevantAccessions[!(RelevantAccessions %in% TairIDs$P48347)]
absenceIDs<-IDs[IDs$P48347 %in% absence,]
#Absent IDs manually found, these were used above
#"A0A2P2CLH9" ATMG00640
#"A0A2P2CLF6" ATMG00070
#"A0A2P2CLF9" ATCG00120
#"A0A2P2CLF3" ATMG00090, I had to blast this one to verify

head(TairIDs)
head (RelevantAccessions)
colnames(TairIDs)[1]<-"Accession"    # change TAIR ID name to V2, same as accession name in file





RootM<-merge(RootMEM.csv,TairIDs,by ="Accession")
RootS<-merge(Membrane,TairIDs,by ="Accession")
ShootS<-merge(Membrane,TairIDs,by ="Accession")
ShootM<-merge(Membrane,TairIDs,by ="Accession")  #all =TRUE will retain values in both ways

write.csv(Membrane, file = "Membrane", sep = ",")




#Compare the columns
compare<-read.csv("C:/Users/Gbolaga/Desktop/BRIC LED RESULTS/Mapped whole protein/Compare_all.csv")

#Find similar proteins between dataset
library(dplyr)
getwd()
compare<-read.csv("compare.csv")
head(compare)
a<-intersect(compare$RootSearth,compare$ShootSearth)  #shared protein in 2 column
length(a)<- length(compare$X)  #force equate the length of column to fill in NA values
compare$a<-a
library(data.table)
setnames(compare,"a","RootSolEarthAl") #to change a single column name in a dataframe
compare<-subset(compare,select= -sharedMemSp)  #to remove column
#find distinct membrane protein in space 
a<-setdiff(compare$RootSearth,compare$ShootSearth) #to find elements unique to [1]alone


#Gene expression volcano plot

vol<-read.csv("Membrane")
names(vol)<-vol[1,]  #####make first row the header
vol<-vol[-1,]  ####delete first row
library(ggplot2)
library(ggrepel)
library(dplyr)
RootM<-read.csv("C:/Users/Gbolaga/Desktop/BRIC LED RESULTS/VOLCANO PLOT/RootMEM",header=TRUE)   ###load full MembraneEM dataset

dim(RootMEM.csv)   #find length and column


###remove first row and make row 2 the header
names(Membrane)<-Membrane[1,]
Membrane<-Membrane[-1,]


##check if field is in numeric
RootMEM.csv$`P-Value`<-as.numeric(RootMEM.csv$`P-Value`)
Membrane$``Adj. P-Value``<-as.numeric(Membrane$`Adj. P-Value`)
Membrane$`Log2 ratio`<-as.numeric(Membrane$`Log2 ratio`)


# Set threshold for diffrentially expressed
ShootSOL.csv$Diffexpressed<-"NO"
ShootSOL.csv$Diffexpressed[ShootSOL.csv$Log2.ratio <=0.2 & ShootSOL.csv$Adj..P.Value<0.05] <- "DOWN"
ShootSOL.csv$Diffexpressed[ShootSOL.csv$Log2.ratio >=0.2 & ShootSOL.csv$Adj..P.Value<0.05] <- "UP"
head(ShootSOL.csv)
ShootSOL.csv$delabel<-NA   #To trick ggplot which requires that we add label to our plot

#Generate volcano plot
plot<-ggplot(data=ShootSOL.csv, aes(x=`Log2.ratio`, y= -log10(`Adj..P.Value`),col=Diffexpressed, label=delabel))+
  geom_point()+ylab("-Log10 (p-value)") +xlab("Log2 Fold Change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line=element_line(colour = "black"))+
  geom_text_repel()+ xlim (-4,5.5) + ylim (0,3.5) +
  scale_color_manual(values=c('red','black','blue'))+
  theme(text=element_text(size = 8))

plot + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12,face="bold")) +
  theme(legend.key.size = unit(0.9, 'cm'))+
  theme(legend.text = element_text(size=12))



#Create plot for all membrane protein in plant

# merge root and shoot data and eliminate replicate
library(data.table)
MEMPro<-rbindlist(list(rootM,ShootMEM.csv),use.names=FALSE)[,lapply(.SD,mean),list(`Log2 ratio`,`Adj. P-Value`,`P-Value`)] 

# for the whole plant. Dataset contains combined root and shoot
Membrane<-read.csv("C:/Users/Gbolaga/Desktop/BRIC LED RESULTS/MembraneP.csv")   ###Abandon this...not biologically feasible



install.packages("rio")    ####This converts all xls in a folder to csv
suppressMessages(require(rio))
setwd("C:/Users/Gbolaga/Desktop/BRIC LED RESULTS/SUBA/")
xls<-dir(pattern = ".xlsx")
created<-mapply(convert, xls, gsub("xlsx", "csv", xls))


##load all  csv files in the folder
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

###to count number of localization in the each datasets and merge result together

library(dplyr)
a<-table(SUBA$Location.Consensus)
a<-data.frame(a)
setnames(a,"Freq","RootSOLUp")
SUBA<-merge(SUBA,a, by="Var1", all= TRUE)

write.csv(SUBA,"SUBA.csv")

#to create a correlation plot between Root and shoot proteins
#select just specific column
library(dplyr)
a<-RootM.csv%>%select("TAIR.ID","Description","Log2.ratio","Adj..P.Value")
b<-ShootM.csv%>%select("TAIR.ID","Description","Log2.ratio","Adj..P.Value")
CorMEM<-merge(a,b, by="TAIR.ID")

library(data.table)
setnames(CorMEM,"Log2.ratio.x","LFC RootMEM")   #change column names
setnames(CorMEM,"Log2.ratio.y","LFC ShootMEM")
suppressMessages(require(ggplot2))

#Make scatter plot
MEM<-ggplot(CorMEM,aes(x=CorMEM$'LFC RootMEM',y=CorMEM$'LFC ShootMEM'))+geom_point()+
  ylab("Log2 Fold change (Shoot proteins)")+xlab("Log2 Fold change (Root proteins)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line=element_line(colour = "black"))

#add lines
MEM+geom_hline(yintercept =0, color= "red")+geom_vline(xintercept =0, color= "red")    # I completed it on excel  

write.csv(CorMEM,"Membrane correlation.csv")



# ## CLUSTERING WITH K-MEAN

##SCRIPT TITLE: K-mean clustering
#AUTHOR: Gbolaga Olanrewaju
#Date: 07/20/2022

#Pre: Find average of ground and space sample, save as csv or do it in script




library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

#load csv data
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

### remove columns and find average of space and ground
##Change name for each
ShootS<-ShootS.csv[,-(24:29)]
ShootS<-ShootS[,-(1:2)]
ShootS<-ShootS[,-(2:15)]


#average of each
ShootS$AvgGround<-rowMeans(ShootS[,c(2:4)])
ShootS$AvgSpace<-rowMeans(ShootS[,c(5:7)])

#Remove individual value column
ShootS<-ShootS[,-(2:7)]

#Remove all TAIR ID

#Clustering
RootM<-na.omit(RootM)
RootMEM<-kmeans(rootM, centers = 4, nstart = 25)
RootSOL<-kmeans(RootS, centers = 2, nstart = 25)
ShootMEM<-kmeans(ShootM, centers = 4, nstart = 25)
ShootSOL<-kmeans(ShootS, centers = 2, nstart = 25)

library(gridExtra)
p1 <- fviz_cluster(RootMEM, geom = "point", data = rootM) + ggtitle("Root MEM")+ xlab("Ground Control") +ylab("Spaceflight")
p2 <- fviz_cluster(RootSOL, geom = "point",  data = RootS) + ggtitle("Root SOL")+ xlab("Ground Control") +ylab("Spaceflight")
p3 <- fviz_cluster(ShootMEM, geom = "point",  data = ShootM) + ggtitle("Shoot MEM")+ xlab("Ground Control") +ylab("Spaceflight")
p4 <- fviz_cluster(ShootSOL, geom = "point",  data = ShootS) + ggtitle("Shoot SOL")+ xlab("Ground Control") +ylab("Spaceflight")

grid.arrange(p1, p2, p3, p4, nrow = 2)



##
rootM<-RootM[,-1]
fviz_cluster(K2,data=rootM) +ggtitle("Root Membrane protein") +
  xlab("Ground Control") +ylab("Spaceflight")



### 7/26/2022
### create a list of TAIR ID, DESCRIPTION, SYMBOLS AND LFC VALUE


