#il est possible sur R de lire des fichiers json, les cluster sont une concatenation de fichiers json dans un fichier java script
#il faut donc transformer ces fichiers pour pouvoir les traiter

library(jsonlite)
library(igraph)

#avant de lancer ces lignes de code, ouvrir les fichiers du cluster, ouvrir data ? l'aide de notepad
#puis ouvrir https://www.convertsimple.com/convert-javascript-to-json/ et coller dans l'onglet adapt? la partie Node pour t?l?chager un fichier json et faire de m?me pour la partie Link
#placer les fichiers json dans le chemin indiqu? et les renommer en Link et Node

setwd("C:/Users/magot2/Desktop/Baptiste")


Link <- fromJSON("linksData.json")
Node <- fromJSON('nodesData.json')


#conditions


seuil_score<-49
seuil_rcov<-69
seuil_fcov<-69
seuil_fmatch<-9
seuil_rmatch<-9

lien<-Link[which((Link$score>=seuil_score)& (Link$fcov>=seuil_fcov|Link$rcov>=seuil_rcov)&(Link$fmatch>=seuil_fmatch|Link$rmatch>=seuil_rmatch)),]


#10 lien
g<-graph_from_data_frame(lien[,c(2,3)],Node$id,directed=F)


liste<-degree(g)

for (i in 1:length(liste)){
  if (liste[i]>10){
    num<-as.numeric(names(liste[i]))#numero de la molécule (id)
    nb<-length(which(lien$source==num))+length(which(lien$target==num)) #nbr de lien
    l<-c(which(lien$source==num),which(lien$target==num)) #liste de lien
    sort.list(lien[l,10],decreasing =TRUE)[11:nb] #lien ordonné par score décroissant
    l[sort.list(lien[l,10],decreasing =TRUE)[11:nb]]
    lien<-lien[-c(l[sort.list(lien[l,10],decreasing =TRUE)[11:nb]]),]
    g<-graph_from_data_frame(lien[,c(2,3)],Node$id,directed=F)
    liste<-degree(g)
  }
}

g<-graph_from_data_frame(lien[,c(2,3)],Node$id,directed=F)
which(degree(g)>10) # doit afficher "numeric(0)" pour indiquer que aucun sommet ne possède plus de 10 liaisons

#100 par clusters





gclus<-clusters(g)#recherche composante connexe du graphe

gclus$no # nbr de cluster

gclus$csize #nbr d'individus dans chaque cluster

Node$cluster<-gclus$membership

l<-V(g)$name[which(gclus$membership == 3)]  #individus dans le cluster 1

lien_clus1<-as.data.frame(c())
for (i in 1:length(l)){
  num<-l[i]
  ll<-c(which(lien$source==num))
  lien_clus1<-rbind(lien_clus1,lien[ll,])
}

g_clus1<-graph_from_data_frame(lien_clus1[,c(2,3)],V(g)$name[which(gclus$membership == 3)],directed=F)

clus1<-cluster_louvain(g_clus1,weights=lien_clus1$score,resolution=1.1)

#faire varier le paramêtre de resolution pour faire varier le nomrbe de sommets par groupe
#pour en avoir moins de 100

barplot(table(factor(clus1$membership, levels = 0:10)),ylim=c(0,60))

plot(g_clus1, vertex.label=NA, vertex.size=5,vertex.color=clus1$membership+1)

#renommage dans le fichier
Node[Node$cluster==3,12]<-clus1$membership+gclus$no+1

barplot(table(factor(Node$cluster, levels = 0:500)),ylim=c(0,100))#vérification nb dans chaque clusters <100

table(factor(Node$cluster))

Node <- as.data.frame(Node)
dim(Node)[1]




#probleme pour enregister en csv pour la colonne 1
#des listes dans un data frame
liste1<-c()
liste2<-c()
liste3<-c()
liste4<-c()
liste5<-c()

for (i in 1:dim(Node)[1]){
  liste1<-c(liste1,Node$areas[[i]][1])
  liste2<-c(liste2,Node$areas[[i]][2])
  liste3<-c(liste3,Node$areas[[i]][3])
  liste4<-c(liste4,Node$areas[[i]][4])
  liste5<-c(liste5,Node$areas[[i]][5])
  
}
Node$areas1<-liste1
Node$areas2<-liste2
Node$areas3<-liste3
Node$areas4<-liste4
Node$areas5<-liste5
Node<-Node[,-c(11)]


#suppression des liens entres les clusters
l_supp<-c()
for (i in 1:dim(lien)[1]){
  if (Node[which(Node$id==lien[i,2]),11]!=Node[which(Node$id==lien[i,3]),11]){
    l_supp<-c(l_supp,i)
  }
}

lien<-lien[-l_supp,]

g<-graph_from_data_frame(lien[,c(2,3)],directed=F)
gclus<-clusters(g)#recherche composante connexe du graphe

gclus$no # nbr de cluster

gclus$csize #nbr d'individus dans chaque cluster

#merge deux tableaux


# tab<-lien[,c(1,3,2,4:14)]
# tab[1]<-tab[1]+87518
# colnames(lien)
# 
# colnames(tab)<-colnames(lien)
# 
# lien<-rbind(lien,tab)
# 




data_cytoscape<-merge(lien,Node,by.y="id",by.x="source")





write.csv(Node,"noeud.csv")
write.csv(lien,"liaison.csv")
write.csv(data_cytoscape[,-c(2,4,5)],"data_cytoscape.csv")



#plot(g<-graph_from_data_frame(lien[,c(2,3)],Node$id,directed=F),vertex.label=NA, vertex.size=1,vertex.color=Node$cluster)
#illisible

