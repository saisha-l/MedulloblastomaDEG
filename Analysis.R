  install.packages("BiocManager")
  install.packages("forcats")
  install.packages("stringr")
  install.packages("readr")
  install.packages("tidyr")
  install.packages("survminer")
  BiocManager::install("limma")
  BiocManager::install("org.Hs.eg.db")
  BiocManager::install("GEOquery")
  BiocManager::install("pheatmap")

library(GEOquery)
library(dplyr)
library(limma)
library(pheatmap)
library(ggplot2)
library(readr)

my_id <- "GSE50161"
gse <- getGEO(my_id)[[1]]
exprs(gse) <- log2(exprs(gse))

dev.new()
boxplot(exprs(gse), outline=FALSE)

sampleInfo <- pData(gse) %>%
  select(source_name_ch1, characteristics_ch1.1) %>%
  rename(group=source_name_ch1, patient=characteristics_ch1.1)

# diff exprs 1
design <- model.matrix(~0+sampleInfo$patient)
colnames(design) <- make.names(colnames(design))

fit <- lmFit(exprs(gse), design)
contrasts <- makeContrasts(Diff = sampleInfo.patienttissue..medulloblastoma - sampleInfo.patienttissue..medulla, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)


full_results <- topTable(fit2, number=Inf) %>%
  tibble::rownames_to_column("ID")

write_csv(full_results, "gse_full_output.csv")

p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>%
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff) %>%
  ggplot(aes(x = logFC, y = B, col = Significant)) + geom_point()


full_results_csv <- read_csv("gse_full_output.csv")
View(full_results_csv)

# randomForest 


library(randomForest)
install.packages("caret")
library(caret)

data <- exprs(gse)
labels <- factor(sampleInfo$group) 

set.seed(123) 
trainIndex <- createDataPartition(labels, p = 0.8, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]
trainLabels <- labels[trainIndex]
testLabels <- labels[-trainIndex]

rf_model <- randomForest(x = trainData, y = trainLabels, ntree = 500, importance = TRUE)

rf_predictions <- predict(rf_model, testData)
confusionMatrix(rf_predictions, testLabels)

importance_scores <- importance(rf_model)
top_genes <- head(sort(importance_scores, decreasing = TRUE), 250) # Adjust the number as needed

View(top_genes)
print(top_genes)
importance_scores <- importance(rf_model)

top_genes_scores <- sort(importance_scores[, "MeanDecreaseAccuracy"], decreasing = TRUE)[1:100]

names(top_genes_scores) <- rownames(importance_scores)[order(importance_scores[, "MeanDecreaseAccuracy"], decreasing = TRUE)][1:100]

top_genes_df <- data.frame(Gene = names(top_genes_scores), Importance = top_genes_scores)

write.csv(top_genes_df, "top_genes.csv", row.names = FALSE)
top_genes_df <- read.csv("top_genes.csv")
head(top_genes_df)
tail(top_genes_df)
View(top_genes_df)
