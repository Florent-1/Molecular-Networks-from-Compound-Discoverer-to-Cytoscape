######## Script for the conversion of Molecular Network from Compound Discovered to Cytoscape #############
#   By Florent Magot    #

#### Script beginning  ####

#Clean the workspace 
rm(list=ls()) 
dev.off

#Install the packages if there are not installed yet
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("stringr")

#Charge the corresponding packages at the beginning of each run
library(dplyr)
library(tidyr)
library(stringr)

#Upload all the data
setwd('C:/Users/magot/Desktop/Papier Romain 2022-03-15/V2')
Clusters = read.table("Clusters.txt", header = FALSE, dec = ".", sep= ",")
SingleNodes = read.table("SingleNodes.txt", header = FALSE, dec = ".", sep= ",")
MetaDataCompounds = read.csv("MetaDataCompounds3.csv", header = TRUE, dec = ".", sep= ",")

#Visualization of the data charged
#head(Clusters)
#summary(Clusters)
#str(Clusters)

#head(SingleNodes)
#summary(SingleNodes)
#str(SingleNodes)

#head(MetaDataCompounds)
#summary(MetaDataCompounds)
#str(MetaDataCompounds)

#For Clusters "id:0" have the name "﻿id:0"
#For SingleNodes "id:XY" have the name "﻿id:XY" (be carful le name of the id change according of the initial data)

Clusters[1,1] <- "id:0"
SingleNodes[1,1] <- "id:12986"


#### Transforming the "Clusters" data frame ####


#Name the colones
NamesClusters <- c("ID_Clusters", "ID", "ID_TARGET", "ID_BIS", "ID_TARGET_BIS",
                   "NAME_TRANSFORMATION", "MASS_TRANSFORMATION", "FORMULA_TRANSFORMATION",
                   "ERROR_TRANSFORMATION", "SCORE", "FCOV", "RCOV", "FMATCH", "RMATCH")

colnames(Clusters) <- NamesClusters

#Generate a new data frame, with the doubled colones removed
Clusters2 <- select(Clusters,-c("ID_BIS", "ID_TARGET_BIS"))

#Remove the text from colones: "SCORE", "FCOV", "RCOV", "FMATCH", "RMATCH"
Clusters3 <- Clusters2
Clusters3$SCORE <- gsub(".*:","",Clusters3$SCORE)
Clusters3$FCOV <- gsub(".*:","",Clusters3$FCOV)
Clusters3$RCOV <- gsub(".*:","",Clusters3$RCOV)
Clusters3$FMATCH <- gsub(".*:","",Clusters3$FMATCH)
Clusters3$RMATCH <- gsub(".*:","",Clusters3$RMATCH)

#Change the nature of the colones "SCORE", "FCOV", "RCOV", "FMATCH", "RMATCH" in numerical
Clusters3$SCORE <- as.numeric(Clusters3$SCORE)
Clusters3$FCOV <- as.numeric(Clusters3$FCOV)
Clusters3$RCOV <- as.numeric(Clusters3$RCOV)
Clusters3$FMATCH <- as.numeric(Clusters3$FMATCH)
Clusters3$RMATCH <- as.numeric(Clusters3$RMATCH)

#Change the name "target" by "source"
Clusters3$ID_TARGET <- str_replace_all(Clusters3$ID_TARGET, "target","source")


####  Selection of the edges ####


#Selection based on the score
Clusters4 <- Clusters3[Clusters3$SCORE >= 50 , ]

#Selction based on the coverage "RCOV" and "FCOV" (R stands for reverse and F for forward)
Clusters4 <- Clusters4[Clusters4$RCOV >= 70 | Clusters4$FCOV >= 70 , ]

#Selection based on the numbers of fragments "FMATCH" and "RMATCH" (R stands for reverse and F for forward)
Clusters4 <- Clusters4[Clusters4$FMATCH >= 11 | Clusters4$RMATCH >= 11, ]

#Remove the column "ID_CLUSTERS", that we are not using
Clusters4$ID_Clusters <-NULL


#### Transforming the "SingleNodes" data frame ####


#head(SingleNodes)
#NamesSingleNodes
#length(NamesSingleNodes)

#Name the colones
NamesSingleNodes <- c("ID", "Nom_de_molecules", "CalcMW","FormulaCD", "Error", "MaxArea", "RT", "NbreFragments",
                      "Structure_SN", "Confidence", "Area_G1", "Area_G2", "Area_G3","Area_G4", "Area_G5",
                      "Area_G6", "Area_G7", "Area_G8","Area_G9", "Area_G10", "Area_G11", "Area_G12","Area_G13",
                      "Area_G14", "Area_G15","Area_G16", "Area_G17", "Area_G18", "Area_G19","Area_G20")

colnames(SingleNodes) <- NamesSingleNodes

SingleNodes2 <- SingleNodes

#Change the name "id" by "source"
SingleNodes2$ID <- str_replace_all(SingleNodes2$ID, "id","source")

#Remove the characters
SingleNodes2$Nom_de_molecules <- gsub(".*:","",SingleNodes2$Nom_de_molecules)
SingleNodes2$FormulaCD <- gsub(".*:","",SingleNodes2$FormulaCD)
SingleNodes2$Error <- gsub(".*:","",SingleNodes2$Error)
SingleNodes2$MaxArea <- gsub(".*:","",SingleNodes2$MaxArea)
SingleNodes2$RT <- gsub(".*:","",SingleNodes2$RT)
SingleNodes2$NbreFragments <- gsub(".*:","",SingleNodes2$NbreFragments)
SingleNodes2$Structure_SN <- gsub(".*:","",SingleNodes2$Structure_SN)
SingleNodes2$Confidence <- gsub(".*:","",SingleNodes2$Confidence)
SingleNodes2$Area_G1 <- gsub("[[]","",SingleNodes2$Area_G1)
SingleNodes2$Area_G1 <- gsub(".*:","",SingleNodes2$Area_G1)
SingleNodes2$CalcMW <- gsub(".*:","",SingleNodes2$CalcMW)
SingleNodes2$Area_G20 <- gsub("]","",SingleNodes2$Area_G20)

#Change to numerical format
SingleNodes2$CalcMW <- as.numeric(SingleNodes2$CalcMW)
SingleNodes2$CalcMW <- format(round(SingleNodes2$CalcMW, 5), nsmall = 5)
SingleNodes2$CalcMW <- as.numeric(SingleNodes2$CalcMW)

SingleNodes2$Error <- as.numeric(SingleNodes2$Error)
SingleNodes2$MaxArea <- as.numeric(SingleNodes2$MaxArea)
SingleNodes2$RT <- as.numeric(SingleNodes2$RT)
SingleNodes2$NbreFragments <- as.numeric(SingleNodes2$NbreFragments)

#Change the format of the data
SingleNodes2$Structure_SN <- as.logical(SingleNodes2$Structure_SN)
#SingleNodes2$Confidence <- as.factor(as.numeric(SingleNodes2$Confidence))
SingleNodes2$Confidence <- as.numeric(SingleNodes2$Confidence)
#unique(SingleNodes2$Confidence)

SingleNodes2$Area_G1 <- as.integer(SingleNodes2$Area_G1)
SingleNodes2$Area_G20 <- as.integer(SingleNodes2$Area_G20)

#Add one colunm with a round value with 4 digits for the molecular weight
SingleNodes2$RoundMW <- round(SingleNodes2$CalcMW,4)


#### Merge the "SingleNodes2" with the "MetaDataCompounds" table ####


SingleNodes3 <- merge(SingleNodes2, MetaDataCompounds, by = c("CalcMW", "RT"), all.x = TRUE)


#### Merge the "SingleNodes2" with "Clusters4"  ####


Cytoscape <- merge(SingleNodes3, Clusters4, by = c("ID"), all.x = TRUE)

Cytoscape$Nom_de_molecules <- str_replace_all(Cytoscape$Nom_de_molecules, ",","-")

Cytoscape[is.na(Cytoscape)] <- 0


############## Save as a .csv file ##################


write.csv(Cytoscape, "RM_Cytoscape.csv")

