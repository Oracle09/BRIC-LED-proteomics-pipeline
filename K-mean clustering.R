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
ShootM<-ShootM.csv[,-(24:29)]
ShootM<-ShootM[,-(1:2)]
ShootM<-ShootM[,-(2:15)]


#average of each
ShootM$AvgGround<-rowMeans(ShootM[,c(2:4)])
ShootM$AvgSpace<-rowMeans(ShootM[,c(5:7)])

#Remove individual value column
ShootM<-ShootM[,-(2:7)]

#Remove all TAIR ID
ShootM<-ShootM[,-1]

#Check methods to determine appropriate number of cluster

library(readr)
scaled_data<-as.matrix(scale(ShootM))
## ELBOW METHOD:::
k.max<-15    #plot graph between cluster of 2 to 15
data<-scaled_data
wss<-sapply(1:k.max,function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
#plot graph
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
fviz_nbclust(ShootM, kmeans, method = "wss")
##SILHOUETTE METHOD

library(purrr)
avg_sil <- function(k) {
  km.res <- kmeans(ShootM, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(ShootM))
  mean(ss[, 3])}
k.values <- 2:15
avg_sil_values <- map_dbl(k.values, avg_sil)
plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")


fviz_nbclust(ShootM, kmeans, method = "silhouette")

##GAP STATISTICS METHOD
gap_stat <- clusGap(RootM, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)

print(gap_stat)
fviz_gap_stat(gap_stat)



#Clustering
ShootSOL<-na.omit(ShootS)
ShootMEM<-kmeans(ShootM, centers =3, nstart = 25)
RootMEM<-kmeans(RootM, centers = 3, nstart = 25)
RootSOL<-kmeans(RootS, centers = 2, nstart = 25)
ShootSOL<-kmeans(ShootS, centers = 2, nstart = 25)

library(gridExtra)
p1 <- fviz_cluster(RootMEM, geom = "point", data = RootM) + ggtitle("Root MEM")+ xlab("Ground Control") +ylab("Spaceflight")
p2 <- fviz_cluster(RootSOL, geom = "point",  data = RootS) + ggtitle("Root SOL")+ xlab("Ground Control") +ylab("Spaceflight")
p3 <- fviz_cluster(ShootMEM, geom = "point",  data = ShootM) + ggtitle("Shoot MEM")+ xlab("Ground Control") +ylab("Spaceflight")
p4 <- fviz_cluster(ShootSOL, geom = "point",  data = ShootS) + ggtitle("Shoot SOL")+ xlab("Ground Control") +ylab("Spaceflight")

grid.arrange(p1,p2,p3,p4, nrow = 2)



##
ShootM<-ShootM[,-1]
fviz_cluster(K2,data=ShootM) +ggtitle("Root Membrane protein") +
  xlab("Ground Control") +ylab("Spaceflight")

 
###Word cloud

BiocManager::install("GOsummaries")
library(GOsummaries) 
library(wordcloud)
library(tm)
library(slam)




knew<-kmeans(ShootM$p_value, centers =3, nstart = 25)

#bind cluster to description
a<-cbind(knew$cluster,ShootM)


##sort by cluster level
library(dplyr)
b<-a%>%
  group_by(`knew$cluster`)%>%
  summarise(description,p_value)
head(b)



###split according to clusters
c<-split(b,f=b$`knew$cluster`)

#set c as dataframes

c1<-as.data.frame(c$`1`)
c2<-as.data.frame(c$`2`)
c3<-as.data.frame(c$`3`)

#Wordcloud
library(ggwordcloud)
library(ggplot2)
ggplot(c2, aes(label = c2$description, size=c2$p_value,color = c2$p_value))+
  geom_text_wordcloud(eccentricity = .35)+
  xlab("Cluster 2")+
  scale_size_area(max_size = 4.8)+
  theme_minimal()+
  scale_color_gradient(low = "grey", high = "black")
 




