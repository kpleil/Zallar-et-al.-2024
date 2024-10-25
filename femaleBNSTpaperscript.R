```{r}
#old stuff 
library(tidyverse)
library(Seurat)
library(calibrate)
library(cowplot)
library(Matrix)
library(readr)
```
```{r}
## FEMALES
# Read in `matrix.mtx`
mtx <- readMM('female/GSE126836_BNST_only_Fem_matrix.mtx.gz')

# Read in `genes.tsv`
genes <- read_tsv("female/GSE126836_BNST_only_Fem_genes.csv.gz", col_names = T)
gene_ids <- genes$gene

# Read in `barcodes.tsv`
cell_ids <- read_tsv("female/GSE126836_BNST_only_Fem_barcodes.csv.gz", col_names = F)$X1
cell_ids <- cell_ids[-1]

# Make the column names as the cell IDs and the row names as the gene IDs
rownames(mtx) <- gene_ids
colnames(mtx) <- cell_ids

# Turn count matrix into a Seurat object (output is a Seurat object)
seurat_female <- CreateSeuratObject(mtx)

saveRDS(seurat_female, "seurat_bnst_female.rds") #this file is female bnst only nuclei
```

```{r}
###new stuff
# load in seurat_bnst_female.rds
#seurat_female <- seurat_bnst_female
#seurat_bnst_female <- NULL

# female only analysis
seurat_female <- NormalizeData(object = seurat_female, normalization.method = "LogNormalize", scale.factor = 10000) 
seurat_female <- FindVariableFeatures(object = seurat_female, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
seurat_female <- ScaleData(object = seurat_female)
seurat_female <- RunPCA(object = seurat_female)
ElbowPlot(seurat_female) # 15 PC


# Determine the K-nearest neighbor graph
seurat_female <- FindNeighbors(object = seurat_female, 
                               dims = 1:15)

# Determine the clusters for various resolutions                                
seurat_female <- FindClusters(object = seurat_female, resolution = 0.3)#changed from .6 before
seurat_female <- RunTSNE(seurat_female, dims = 1:15)
seurat_female <- RunUMAP(seurat_female, dims = 1:15)

saveRDS(seurat_female, "seurat_bnst_female_analyzed.rds") #this file is female bnst only analyzed seurat 

DimPlot(seurat_female, reduction = 'tsne', label = TRUE, label.size = 6)
```
```{r}
FeaturePlot (seurat_female, c('Crh'), reduction = 'umap', label = TRUE, label.size = 5, pt.size = .2)
FeaturePlot (seurat_female, c('Gfap'), reduction = 'umap', label = TRUE, label.size = 5, pt.size = .2)
FeaturePlot (seurat_female, c('Slc17a6'), reduction = 'umap', label = TRUE, label.size = 5, pt.size = .2)
```


```{r}
FCRH_subset <- subset(seurat_female, Crh>0)
FCRH_subset <- RunUMAP(FCRH_subset, dims = 1:12)
```

```{r}
Fvglut_subset <- RunUMAP(Fvglut_subset, dims=1:15)
Fvglut_subset <- subset(seurat_female, Slc17a6>0)
FeaturePlot(Fvglut_subset, 
            c('Esr1' , 'Esr2'), cols = c("lightgrey", "#FF10F0"),
            reduction = 'umap',  , label.size = 5, label = T, order=T,
            pt.size = 1.7)

```{r}
FGFAP_subset <- subset(seurat_female, Gfap>0)
FGFAP_subset <- RunUMAP(FGFAP_subset, dims = 1:5)


FeaturePlot(FGFAP_subset, 
            c('Esr1' , 'Esr2'), cols = c("lightgrey", "#FF10F0"),
            reduction = 'umap',  , label.size = 5, label = T, order=T,
            pt.size = 1.7)

```