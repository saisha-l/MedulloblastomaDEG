install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GEOquery")
library(GEOquery)
my_id <- "GSE50161"
gse <- getGEO(my_id)
length(gse)
gse <- gse[[1]]
gse
summary(exprs(gse))

exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)
library(dplyr)
sampleInfo <- pData(gse)
sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)
sampleInfo <- rename(sampleInfo, group=source_name_ch1, patient=characteristics_ch1.1)
sampleInfo


BiocManager::install("pheatmap")
library(pheatmap)
corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)
rownames(sampleInfo)
colnames(corMatrix)

rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo) 


library(readr)
full_output <- cbind(fData(gse),exprs(gse))
write_csv(full_output, "gse_full_output.csv")
features <- fData(gse)
View(features)
library(limma)
design <- model.matrix(~0+sampleInfo$patient)
design 
colnames(design) <- make.names(colnames(design))
summary(exprs(gse))
design

cutoff <- median(exprs(gse))

## is gene expressed in sample
is_expressed <- exprs(gse) > cutoff

## identify genes expressed in 2+ samples

keep <- rowSums(is_expressed) > 2

table(keep)

gse <- gse[keep,]

fit <- lmFit(exprs(gse), design)
head(fit$coefficients)
contrasts <- makeContrasts(sampleInfo.patienttissue..medulloblastoma - sampleInfo.patienttissue..medulla, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
topTable(fit2, coef=1)
anno <- fData(gse)
anno
print("hello everyone")
full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()


p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

