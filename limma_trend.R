
### limma-trend pipeline for differential expression analysis

# create the count dataframe and design matrix

# count dataframe consists of sample names in the columns, genes in the rows, and gene counts before normalization as values

library(readxl)
#define the directory to the counts file (counts_path)
counts <- read_excel(counts_path)
counts <- as.data.frame(counts)
row.names(counts) <- counts[,1]
counts <- counts[,-1]

# design matrix consists of two columns: sample names and the class of data which each sample belongs to (as 0 (e.g, healthy) or 1 (e.g, disease))
# shows the design of the experiment

#define the directory to the design file (design_path)
design_matrix <- read_excel(design_path)
design <- as.data.frame(design_matrix)
row.names(design) <- design[,1]
design <- design[,-1]

# create logCPM (logarithm of counts per million) from counts

library(edgeR)
dge <- DGEList(counts=counts)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log=TRUE, prior.count=3)

# Differential expression with the limma-trend pipeline (limma voom differs from limma trend from this point)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
#define n_top_genes
top_table <- topTable(fit, coef=ncol(design), number = n_top_genes)
# output: genes with most significant expression changes and their significance