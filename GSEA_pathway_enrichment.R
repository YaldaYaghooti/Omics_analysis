##Pathway enrichment analysis with GSEA

#load list of genes sorted based on their fold change between two experimental groups. 
#The file consist of two columns: gene_names and fold_change in a descending order.
library(readxl)
sorted_list <- read_excel(FCsorted_path)

mean_vector <- sorted_list$fold_change`
names(mean_vector) <- sorted_list$gene_names

require(fgsea)
require(dplyr)

#load a file containing the gene sets in a text format. The file should contain name of the sets as rows and each gene contained in the set in a column (a gene set consisting of 5 genes should be represented in a single row and 5 columns).
#the file is also called a gmt file.
myGO = gmtPathways(gmt_path)

#Enrichment analysis
fgRes_mean <- fgsea(pathways = myGO, stats = mean_vector)
fgress <- apply(X =fgRes_mean,FUN = as.character, MARGIN = 2)

#save the results
#define the destination directory as results_path
write.csv(as.data.frame(fgress), results_path, row.names=FALSE)

