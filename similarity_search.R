# Installation of R via RStudio;
# Installation of ChemminerR in RStudio: install.packages("BiocManager", dependencies=TRUE)
# BiocManager::install("ChemmineR")
# BiocManager::install("ChemmineOB")
# Installation of gplots in RStudio: install.packages("gplots")

######################
#Similarity searching#
######################

#Load libraries:
library(ChemmineR)
library(gplots)
library(stats)

# load query compounds from SDF file
sdf <- read.SDFset("C:\\Users\\woutv\\Desktop\\ARDI_chirality_BRENKS_NIH_PAINS_ZINC_filtered_ENAMINE.sdf")

# assign unique ID's
sdfid(sdf)
unique_ids <- makeUnique(sdfid(sdf))
cid(sdf) <- unique_ids

# atom pair descriptors generation
ap <- sdf2ap(sdf)
ap

# generate fingerprints from atom pair descriptors
fp <- desc2fp(ap)
fp

# distance matrix
disMA <- cmp.cluster(db=ap, cutoff=0, save.distances="distmat.rda", quiet=TRUE) 
load("distmat.rda") 

# clustering (cluster size is the amount of compounds and count is the amount of clusters with this size)
clusters <- cmp.cluster(db=ap, cutoff = c(0.7, 0.8, 0.9), quiet = TRUE)
cluster.sizestat(clusters, cluster.result = 1)
cluster.sizestat(clusters, cluster.result = 2)
cluster.sizestat(clusters, cluster.result = 3)

# similarity matrix and clustering
simMA <- sapply(cid(sdf), function(x) fpSim(fp[x], fp, sorted=FALSE))
hc <- hclust(as.dist(1-simMA), method="average")

# dendrogram (export 20000 width, otherwise not properly visible labels)
dendro <- as.dendrogram(hc)
plot(dendro, edgePar=list(col=4, lwd=1), horiz=TRUE, axes = FALSE)

# heatmap using '1-distmat' or 'simMA' (export 20000 width, otherwise not properly visible labels)
heatmap.2(simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), 
          col=colorpanel(40, "red", "white", "blue"), 
          density.info="none", trace="none")

####################
#Compound selection#
####################

# clustering
k <- 10
kmeans_result <- kmeans(simMA, centers = k, nstart = 30)
cluster_labels <- kmeans_result$cluster

# selection
num_compounds_to_select <- 2
selected_compounds <- list()

for (i in 1:k) {
  cluster_indices <- which(kmeans_result$cluster == i)
  cat("cluster_indices length:", length(cluster_indices), "\n")
  if (length(cluster_indices) > 1) {
    cluster_compounds <- sdf[cluster_indices]
    similarity_matrix <- simMA[cluster_indices, cluster_indices]
    cat("similarity_matrix dimensions:", dim(similarity_matrix), "\n")
    selected_indices <- integer(num_compounds_to_select)
    # calculate dissimilarity matrix between compounds in the cluster
    cluster_dissimilarity <- 1 - similarity_matrix
    # calculate the sum of dissimilarities for each compound
    compound_dissimilarity <- rowSums(cluster_dissimilarity)
    # find the two compounds with the highest combined dissimilarity
    selected_indices[1:2] <- order(compound_dissimilarity, decreasing = TRUE)[1:2]
    # add the selected compounds to the list
    selected_compounds[[i]] <- sdfid(cluster_compounds[selected_indices])
  } else {
    # if there is only one compound in the cluster, select it
    selected_compounds[[i]] <- sdfid(sdf[cluster_indices])
  }
}

# combine the selected compounds from each cluster
selected_compounds <- unlist(selected_compounds)
selected_df <- data.frame(compounds = selected_compounds, constant = "selected")
write.table(selected_df, file = "selected_20_compounds.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
