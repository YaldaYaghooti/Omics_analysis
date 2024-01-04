##hierarchical clustering for gene module detection in each experimental group

#read data
#define directory to the logCPM file. You can see the code for creating the logCPM file in "limma_trend" in this repository.
library(readxl)
logCPM <- read_excel(logCPM_path)
logCPM <- as.data.frame(logCPM)
row.names(logCPM) <- logCPM[,1]
logCPM <- logCPM[,-1]

#create correlation martix
logCPM_t <- t(logCPM)
cor_matrix <- cor(logCPM_t, method = 'pearson')
diag(cor_matrix) <- 0

#create dissmiliraity matrix
dissim_mat <- as.dist(1-as.data.frame(cor_matrix))

#create gene tree with hierarchical clustering
gene_tree <- hclust(dissim_mat_smad, method = "average")

#detect gene modules by setting a threshold (h)
modules <- cutree(gene_tree, h=0.65)

#save file
write.csv(modules, file = results_path, row.names = TRUE)
