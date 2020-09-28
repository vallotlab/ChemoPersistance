#Global environment variable

working_directory = "/media/pacome/LaCie/InstitutCurie/Documents/GitLab/ChemoPersistance/"

pcaText <- TRUE
annotText <- "Sample"
hcText <- "Name"  ## column used in hierarchical clustering # change names from Sample_x to actual sample name
centering <- c("none","mean","median")[1]

##Hierarchical clustering
distHC <- c("distPearson","distCosine","euclidean","maximum","manhattan","canberra","binary","minkowski")[1]
methHC <- c("ward","ward.D","ward.D2","single","complete","average")[2]

##ConsClust
repsCC <- 1000
pItemCC <- 0.8
pFeatureCC <- 1
clusterAlgCC <- c("hc","pam","km","kmdist")[1]
distCC <- c("pearson","distCosine","euclidean","manhattan")[1]
innerLinkageCC <- c("ward","ward.D","ward.D2","single","complete","average")[2]
finalLinkageCC <- c("ward","ward.D","ward.D2","single","complete","average")[2]

## Heatmap
chRangeHM <-FALSE # Should be set to TRUE for expression data, FALSE for methylation data
hmColors <- colorRampPalette(c("royalblue","white","indianred1"))(256)
corColors <- colorRampPalette(c("indianred","white","forestgreen"))(256)


#metadataFile <- file.path(ProjectDir, "Scripts", expType, paste0("metadata_",expType,"_", annoType, ".txt") )
maxKHC <- 10
maxKCC <- 10

