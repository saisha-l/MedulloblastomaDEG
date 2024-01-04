# Identifying Differentially Expressed Genes for Medulloblastoma Cancer

## Background 
Medulloblastomas are a type of brain tumor predominantly found in the cerebellum of  pediatric patients. Given the vital role of the cerebellum in bodily functions, these tumors can cause neurological deterioration, and in some cases are life-threatening. Comparing the sequenced gene expression of those with Medulloblastomas and healthy patients, can help us identify specific biological pathways that may be affected by the disease. 

## Methods
### Data
I compiled my dataset by combining various experimental datasets from Gene Expression Omnibus, an open source functional genomics data repository. Multiple datasets were used in order to decrease experimental bias, and get more accurate results. 

### Analysis using R
Data was cleaned and preprocessed before using "limma" to perform differential expression analysis, by fitting a linear model. The results were adjusted and filtered by p-value.


### GEO2R
Differentially expressed genes were also found through using GEO2R, _____. These results were compared with the original R analysis genes, in order to find the top 250 genes (sorted by P-Value) to be differentially expressed


## Results
Full list of genes can be found in the significant_genes.csv file. 
Top 250 Genes were used to perform functional enrichment analysis with DAVID and Metascape. 

