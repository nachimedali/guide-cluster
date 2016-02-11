# Mohamed Ali Nachi
# 07/02/2016
# Desktop Genetics

# Importing Database

db<-read.csv("C:/Users/Dally/Desktop/Deskgen/data/data.tsv",sep="\t")

# Converting database to str

str(db)

#First Step : Calculating medians for t1 = 7 Days, t2 = 14 Days. 
D7_t1_without_drugs<-c()
D14_t2_without_drugs<-c()
D7_t1_with_drugs<-c()
D14_t2_with_drugs<-c()

for(i in 1:nrow(db))
  D7_t1_without_drugs[i]<-(db[i,5]+db[i,6])/2

for(i in 1:nrow(db))
  D14_t2_without_drugs[i]<-(db[i,7]+db[i,8])/2

for(i in 1:nrow(db))
D7_t1_with_drugs[i]<-(db[i,9]+db[i,10])/2

for(i in 1:nrow(db))
D14_t2_with_drugs[i]<-(db[i,11]+db[i,12])/2

# Generating new_db containing genes_id, and medians of time
 new_db<-data.frame(D7_t1_without_drugs,D14_t2_without_drugs,D7_t1_with_drugs,D14_t2_with_drugs)
 
# Generating a dendrogram to deduce the number of clusters that need to be used with hierarchical clustering
# Kindly not that we used only 20000 samples of the database due to our computer performance, to use this 
# algorithm for the hole db, we need only to delete the next two lines and replace "data3_new" by new_db

 tri_alea<-sample(1:nrow(new_db),20000,replace=F)
 data3_new<-new_db[tri_alea,]
 class_new<-scale(data3_new, center = TRUE, scale = TRUE)
 dc<-dist(class_new, method ="euclidean", diag=FALSE, upper=FALSE)
 hier<-hclust(dc,"ward.D")
 
# After deducing that clusters number is 4, we fix this number to be able to cluster our data with a hiearchical 
# clustering approach
 nb_clusters<-4
 X11()
 plot(hier,hang=-1)
 cl4<-cutree(hier,nb_clusters)
 rect.hclust(hier,k=nb_clusters)
 
 # number of element per cluster
 a<-table(cl4)
 hier$order
 a1<-c()
 for(k in 1:length(a)) a1[k]<-a[[k]]

 #Render the number of genes per cluster
a1
 
 
 # Caracterize every cluster with genes that belong to it
 uu2<-data.frame(data3_new,cluster=as.data.frame(cl4))
 aa<-list(data.frame()) 
 for(k in 1:length(uu2[,1])) aa[[k]]<-uu2[uu2$cl4==k,]

 # To show gene sequences that belong to a cluster
 
 View(aa[[4]])
 
 #Retrieving genes id that are present in every cluster and save the data in csv file
 as.numeric(rownames(aa[[4]]))
 write.table(as.numeric(rownames(aa[[1]])), "clu1.csv", sep="\t", row.names=FALSE, col.names=FALSE)
 write.table(as.numeric(rownames(aa[[2]])), "clu2.csv", sep="\t", row.names=FALSE, col.names=FALSE)
 write.table(as.numeric(rownames(aa[[3]])), "clu3.csv", sep="\t", row.names=FALSE, col.names=FALSE)
 write.table(as.numeric(rownames(aa[[4]])), "clu4.csv", sep="\t", row.names=FALSE, col.names=FALSE)
 
##Defining id for every cluster
 lis<-list()
 for (i in 1:nb_clusters)
  lis[[i]]<-as.numeric(rownames(aa[[i]]))

# Retrieving random vector spacer 
 vect<-as.character(db$spacer_seq[tri_alea])
 

 # Getting sequences genes of every gene in cluster i
 
liste_spacer_seq<-list()
for (i in 1:nb_clusters)
liste_spacer_seq[[i]]<-as.character(db$spacer_seq[as.numeric(rownames(aa[[i]]))])

### Viewing genes sequences
liste_spacer_seq[[1]]

########################### Second Step, PCA approach
#spacer_seq==> stat descp
####create function dissociate spacer...
dissos_spacer_seq<-function(indice,num_cluster)
{
  a<-as.character(liste_spacer_seq[[num_cluster]][indice])
  b<-unlist(strsplit(a,""))
  c<-table(b)
return(c)
}

### Running the function for the 100th spacer_seq in cluster 1
dissos_spacer_seq(num_cluster = 1,indice = 100)

cvb<-data3_new

library(ade4) 
library(FactoClass) 
library(FactoMineR)

AFC <- dudi.coa(df=cvb, scannf=FALSE, nf=ncol(cvb)) 
X11()
plot.dudi(AFC) 
distMat <- dist.dudi(AFC, amongrow=TRUE) 
CAH <- ward.cluster(distMat, peso = apply(X=cvb, MARGIN=1, FUN=sum) , plots = TRUE, h.clust = 1) 
cvb$clusters <- cutree(tree = CAH, k = 4) 
s.class(cstar=1,addaxes=TRUE, grid=TRUE, axesell=TRUE, 
        dfxy=AFC$li, fac=as.factor(cvb$clusters), col=1:4, 
        label=c(1:4), csub=1.2, possub="bottomright") 
plot(as.dendrogram(CAH), leaflab = "none") 
X11()
plot.dudi(AFC, Tcol = TRUE, Trow = FALSE)  
s.class(cstar=1,addaxes=TRUE, grid=FALSE, axesell=TRUE, 
        dfxy=AFC$li, fac=as.factor(cvb$clusters), col=1:4, 
        label=c(1:4), csub=1.2, possub="bottomright", add=TRUE)  
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R") 
ordreClasses <- unique(cvb$clusters[CAH$order]) 
A2Rplot(x = CAH, k = 4, boxes = FALSE,  col.up = "gray50", col.down = c(1:4)[ordreClasses], show.labels = FALSE, main = "Dendrogramme")

### Input used is the genes sequences

determ_per_couples_per_cluster<-function(num_cluster)
{
data<-as.data.frame(liste_spacer_seq[[num_cluster]])
colnames(data)<-"cluster"

data$cluster<-as.character(data$cluster)

x<-as.character(data$cluster[1])

####if the data is from R (ME)
mydecoup=function(x){
  t=unlist(strsplit(x,split =""))
  return(c(t))
}

###recup dissoc seq_genes elements into list
var<-data$cluster
list_l<-list()
for(i in 1:nrow(data))
  list_l[[i]]<-var[i]
list_l_dissoc<-lapply(list_l,mydecoup)
############################

a<-list_l_dissoc[[20]]
recup_couple_genes<-function(a)
{
  k<-list()
  for(i in 1:(length(a)-1))
    k[[i]]<-paste(a[i],a[i+1],sep=",")
  return(unlist(k))
}
####test the function
recup_couple_genes(a)

###the list list_l_dissoc_per_couples contains the couples for each seq of genes
#using the ith cluster !
list_l_dissoc_per_couples<-lapply(list_l_dissoc,recup_couple_genes)
##################Description of couples that occurs in the ith cluster 
liste_of_all_couples<-unlist(list_l_dissoc_per_couples)
A<-rev(sort(table(liste_of_all_couples)))
###using percent... 

return(round((A/sum(A))*100,3))
#end
}

for(i in 1:nb_clusters)
print(determ_per_couples_per_cluster(num_cluster = i))

##### Splitting letters in the same cluster

determ_per_elements_per_cluster<-function(num_cluster)
{
data<-as.data.frame(liste_spacer_seq[[num_cluster]])
colnames(data)<-"cluster"

data$cluster<-as.character(data$cluster)

x<-as.character(data$cluster[1])
####if the data is from R (ME)
mydecoup=function(x){
  t=unlist(strsplit(x,split =""))
  return(c(t))
}

###recup dissoc seq_genes elements into list
var<-data$cluster
list_l<-list()
for(i in 1:nrow(data))
  list_l[[i]]<-var[i]
list_l_dissoc<-lapply(list_l,mydecoup)

return(round(rev(sort(table(unlist(list_l_dissoc)))),3)*100/sum(table(unlist(list_l_dissoc))))
}


for(i in 1:nb_clusters)
print(determ_per_elements_per_cluster(num_cluster = i))