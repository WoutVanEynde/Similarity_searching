# Installation of R via RStudio;
# Installation of ChemminerR in RStudio: install.packages("BiocManager", dependencies=TRUE)
# Installation of packages in terminal: sudo apt-get install librsvg2-dev libeigen3-dev
# Installation of openbabel (not with sudo apt-get install!)
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
sdf <- read.SDFset("/home/wout/AR/6_similarity_search/ARDI_chirality_BRENKS_PAINS_filtered_ENAMINE.sdf")

# assign unique ID's
sdfid(sdf)
unique_ids <- makeUnique(sdfid(sdf))
cid(sdf) <- unique_ids

# MACCS fingerprints generation
fp <- fingerprintOB(sdf, "MACCS")
fp

# similarity matrix from MACCS fingerprints
simMA <- sapply(cid(fp), function(x) fpSim(fp[x], fp, method="Tanimoto", sorted=FALSE)) 
disMA <- 1-simMA
hc <- hclust(as.dist(disMA), method="single")
hc[["labels"]] <- cid(fp)

# dendrogram (export 20000 width, otherwise not properly visible labels)
dendro <- as.dendrogram(hc)
plot(dendro, edgePar=list(col=4, lwd=2), horiz=TRUE, axes = FALSE)

# heatmap (export 20000 width, otherwise not properly visible labels)
heatmap.2(simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), 
          col=colorpanel(40, "red", "white", "blue"), 
          density.info="none", trace="none", labRow = cid(fp), labCol = cid(fp))

####################
#Compound selection#
####################

# clustering
k <- 20
kmeans_result <- kmeans(simMA, centers = k, nstart = 100)
table(kmeans_result$cluster)
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
    list(compound_dissimilarity)
    # find the compound with the highest combined dissimilarity
    selected_indices[1] <- order(compound_dissimilarity, decreasing = TRUE)[1]
    selected_indices[1]
    # find second compound with highest dissimilarity compared to first one
    second_compound_dissimilarity <- as.vector(cluster_dissimilarity[,selected_indices[1:1]])
    list(second_compound_dissimilarity)
    selected_indices[2] <- order(second_compound_dissimilarity, decreasing = TRUE)[1]
    selected_indices[2]
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
write.table(selected_df, file = "selected_40_compounds.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
