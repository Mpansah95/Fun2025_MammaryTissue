# Loading the libraries
library(edgeR)
library(limma)
library(ggplot2)

setwd("/Users/princemensahansah/Library/CloudStorage/Box-Box/SpringClass/BIOL 6850/individualproject/")

# Loading the count matrix and sample metadata
# Counts: genes are in rows, samples are in columns
counts <- read.table("gene_count_matrix.csv", header = TRUE, row.names = 1, sep = ",")
meta <- read.table("phenotype.txt", header = TRUE)
View(meta)
colnames(counts)

# Checking to see differences in count and metadata
counts <- counts[, meta$ID]

keep <- filterByExpr(dge, group = meta$Size)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Creating a DGEList object
dge <- DGEList(counts = counts, group = meta$Size)

barplot(dge$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")


cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)

keep.exprs08 <- rowSums(cpm>1)>=3
x_filtered_cpm1_3 <- dge[keep.exprs08, keep.lib.sizes=FALSE]
dim(x_filtered_cpm1_3)


x2_2008<- x_filtered_cpm1_3

x2_2008norm <- calcNormFactors(x2_2008, method = "TMM")
x2_2008norm$samples$norm.factors

x2_2008norm

library(RColorBrewer)

lcpm <- cpm(x2_2008norm, log=TRUE)
par(mfrow=c(1,1))
col.group <- meta$Size
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")

col.group <- as.character(col.group)
plotMDS(lcpm, labels=meta$Size, col=col.group)
title(main="A. Treatments")


# Setting up the design matrix
meta$Size <- factor(meta$Size, levels = c("medium", "small"))
design <- model.matrix(~Size, data = meta)

# Applying the voom transformation

v <- voom(dge, design, plot = TRUE)

# Fitting the linear model
fit <- lmFit(v, design)

# Applying empirical Bayes moderation
fit <- eBayes(fit)

# Extracting results (medium vs small)
results <- topTable(fit, coef = 2, adjust = "BH", number = Inf)
dim(results)

DEGs <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
dim(DEGs)

library(ggplot2)

results$Significant <- ifelse(results$adj.P.Val < 0.05 & results$logFC > 1, "Upregulated",
                              ifelse(results$adj.P.Val < 0.05 & results$logFC < -1, "Downregulated", "Ns"))
logfc_threshold <- 1
hline <- log10(0.05)
#png("Volcano_Dr.no_initialfilter.png", width = 6, height = 6, units = "in", res = 300)
ggplot(results, aes(x=logFC, y=-log10(adj.P.Val), color=Significant)) +
  geom_point(alpha=0.8) +
  scale_color_manual(values=c("red", "grey","blue")) +
  theme_classic() +
  labs(x="Log2 Fold Change", y="-Log10 ( P-Value)") +
  geom_vline(xintercept=c(-logfc_threshold, logfc_threshold), linetype="dashed", color="blue")+
  geom_hline(yintercept=-hline,linetype="dashed", color="red")+
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()
# Writing the differentials into a csv files
write.csv(DEGs, "limma_voom_differential_expression_schwartz2.csv")

# Optional: basic volcano plot
volcanoplot(fit, coef = 2, highlight = 10, names = rownames(fit))




library(dplyr)
library(stringr)
library(UpSetR)
limma <- read.csv("~/Library/CloudStorage/Box-Box/SpringClass/BIOL 6850/Individualproject/limma_voom_differential_expression_schwartz2.csv")
deseq <- read.csv("~/Library/CloudStorage/Box-Box/SpringClass/BIOL 6850/Individualproject/DESeq2_significant_results.csv")
dim(limma)
dim(deseq)

Limma_upregulated <- limma %>%
  filter(logFC > 1, adj.P.Val < 0.05)

Limma_Downregulated <- limma %>%
  filter(logFC < -1, adj.P.Val < 0.05)

Deseq2_upregulated <- deseq %>%
  filter(log2FoldChange > 1, padj < 0.05)

Deseq2_Downregulated <- deseq %>%
  filter(log2FoldChange < -1, padj < 0.05)



limma$X <- sapply(strsplit(limma$X, "\\|"), `[`, 2)

# samples <- list(
#   'Limma' = toupper(as.character(limma$X)),
#   'Limma Downregulated' = toupper(as.character(Limma_Downregulated$X)),
#   'Limma Upregulated' = toupper(as.character(Limma_upregulated$X)),
#   'DESeq2' = toupper(as.character(deseq$X)),
#   'DESeq2 Upregulated' = toupper(as.character(Deseq2_upregulated$X)),
#   'DESeq2 Downregulated' = toupper(as.character(Deseq2_Downregulated$X))
#   
#   )

# Making a list for Upset plot

samples <- list(
  'Limma' = toupper(as.character(limma$X)),
  'DESeq2' = toupper(as.character(deseq$X))
)

binaryData <- fromList(samples)
binaryData <- as.data.frame(lapply(binaryData, as.numeric))

png("Limma vs DESeq2.png", width = 6, height = 6, units = "in", res = 300)
upset(binaryData, nsets = 11,matrix.color = "blue",sets.bar.color = "red",
      mainbar.y.label = "Intersection DEGs", 
      order.by ="freq" )

dev.off()
