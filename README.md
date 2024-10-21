![alt text]([http://url/to/img.png](https://github.com/jonasns/scODIN/blob/main/scODIN_cover_image.jpeg)

# scODIN: Optimized Detection and Inference of Names in scRNA-seq data

`scODIN` is an R package designed for analyzing single-cell RNA sequencing (scRNA-seq) data using the scODIN scoring method. It provides tools for scoring cell types based on priority genes, clustering, simplifying labels, downsampling, and performing k-Nearest Neighbors (kNN) classification

## Installation

You can install the `scODIN` package from GitHub using the `devtools` package:

```r
devtools::install_github("jonasns/scODIN")
```

## Dependencies

please install the following dependencies separately according to their instructions:

`Seurat`

`dplyr`

`readxl`

`stringr`


## Tutorial and presentation of functions
### Load the package and dependencies

```r
library(scODIN)
library(Seurat)
library(readxl)
library(ggplot2)
```

### Load test data
First download the dataset from here: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

```r
# Load the PBMC dataset
pbmc_3k <- Read10X(data.dir = "~/Desktop/filtered_gene_bc_matrices/hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc_3k <- CreateSeuratObject(counts = pbmc_3k, project = "pbmc3k", min.cells = 3, min.features = 200)

#QC
pbmc_3k[["percent.mt"]] <- PercentageFeatureSet(pbmc_3k, pattern = "^MT-")
pbmc_3k <- subset(pbmc_3k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalize Seurat object
pbmc_3k <- NormalizeData(pbmc_3k, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable features
pbmc_3k <- FindVariableFeatures(pbmc_3k, selection.method = "vst", nfeatures = 2000)

#scale the Seurat data
all.genes <- rownames(pbmc_3k)
pbmc_3k <- ScaleData(pbmc_3k, features = all.genes)
```

### UMAP of the data
```r
pbmc_3k <- RunPCA(pbmc_3k, features = VariableFeatures(object = pbmc_3k))
pbmc_3k <- FindNeighbors(pbmc_3k, dims = 1:10)
pbmc_3k <- FindClusters(pbmc_3k, resolution = 0.5)
pbmc_3k <- RunUMAP(pbmc_3k, dims = 1:10)
DimPlot(pbmc_3k, reduction = "umap")
```

### Add the original cell type IDs to the Seurat object
How it was in the original paper: (https://www.nature.com/articles/ncomms14049)

```r
#Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc_3k)
pbmc_3k <- RenameIdents(pbmc_3k, new.cluster.ids)
```

```r
# plot
DimPlot(pbmc_3k, reduction = "umap", label = TRUE, repel = T, pt.size = 0.5, raster = F) + ggtitle("original annotation") +
  theme(plot.title = element_text(hjust = 0.5))
```

### scODIN scoring on the top level
**scODIN function #1 - odin_scoring()**  
This function calculates scODIN scores for cell types based on priority genes from a provided gene priority table and a Seurat object with gene expression values.  

@param `gene_priority_table` A tibble containing gene priorities for different cell types. Usually provided in an Excel sheet.  
@param `scData` A Seurat object containing normalized and scaled scRNA-seq gene expression data.  
@param `core_cell_cutoff` A numeric value indicating the cutoff for what is considered a core cell type. Default is 5, but should be adjusted lower if using less deeply sequenced dataset.  
@param `similarity_threshold` A numeric value defining the similarity threshold for determining double labels. Default is 1.  
@param `accepted_doubles_table` A tibble for accepted double labels. Usually provided in an Excel sheet.  
@param `cell_level` A string telling which level the analysis should be done on (e.g., "Top", "CD4_T").  

@return A Seurat object with scODIN scores and metadata added.  

```r
# Note the gene_priority_table has been provided to reviewers and will be made available upon paper acceptance.
gene_priority_table <- read_excel("~/Desktop/241016_TableS1_gene_priority_table.xlsx")
accepted_doubles_table <- read_excel("~/Desktop/241016_TableS1_gene_priority_table.xlsx", sheet = "accepted_doubles")

scData = odin_scoring(gene_priority_table = gene_priority_table,
                       scData = pbmc_3k,
                       cell_level = "Top",
                       core_cell_cutoff = 2,
                       similarity_threshold = 1,
                       accepted_doubles_table = accepted_doubles_table) 
```

### Combining cell number and octim_score per cluster for a cluster-labeling
**scODIN function #2 - odin_cluster_scoring**  
Combine cell number and scODIN score for labeling cells on the cluster level  
This function calculates a combined score based on cell number and scODIN score per cluster, and assigns the top scoring label to each cluster in the given Seurat object.
 
@param `scData` A Seurat object containing the clustering data and scODIN scores.  
@param `clustering_column` The name of the column that contains the clustering information (default is "seurat_clusters").  

@return The Seurat object with an updated `odin_classification` column, containing the highest scoring label for each cluster.

```{r}
scData <- odin_cluster_scoring(scData, "seurat_clusters")
```

```{r}
DimPlot(scData, reduction = "umap", label = T, repel = T, group.by = 'odin_classification')      + ggtitle("scODIN_results") +
  theme(plot.title = element_text(hjust = 0.5))   
```

### Now subset to CD4 T cells
```{r}
scCD4 = subset(scData, final_labels == "CD4_T")
```

### rerun scODIN scoring on CD4 T cell level
```{r}
scCD4 = odin_scoring(gene_priority_table = gene_priority_table,
                       scData = scCD4,
                       cell_level = "CD4_T",
                       core_cell_cutoff = 2,
                       similarity_threshold = 1,
                       accepted_doubles_table = accepted_doubles_table) 
```

### Simplify final double labels
**scODIN function #3 - simplify_double_labels()**  
Simplify final double labels in scRNA-seq data  
This function simplifies the final double labels in a Seurat object by replacing accepted double combinations of cell types with a simplified name provided by the user.  

@param `scData` A Seurat object that contains `final_labels` in its metadata.  
@param `accepted_doubles_table` A tibble for accepted double labels. Usually coming from an Excel sheet.  
@param `cell_level` Character string specifying the subset level (e.g., "CD4_T").  
 
@return The Seurat object with updated `final_labels` in the metadata.

```{r}
# labels before
table(scCD4$final_labels)

scCD4 = simplify_double_labels(scCD4, accepted_doubles_table, "CD4_T")

# labels after
table(scCD4$final_labels)
```

In the previous steps we have identified the core cells. In the next steps we will conduct kNN search for nearest neighbours. First we downsample and run PCA on the downsampled dataset.  

### Downsample core cells
**scODIN function #4 - downsample_core()**  
Downsample the identified core cells to a specified cell count per final label.  
This function downsamples a Seurat object to have a similar number of core cells per final label. It randomly selects cells up to a user-defined target count or the median cell count if specified.  

@param `scData` A Seurat object containing the single-cell data and the core cell types as `final_labels`.  
@param `target_count` An integer specifying the target count for downsampling. It is overruled by use_median = TRUE.  
@param `use_median` A logical indicating whether to use the median cell count (TRUE) or the target count (FALSE).  
@param `labels_to_exclude` A character vector specifying the labels to exclude (default is `c("unknown")`).  
@param `seed` An integer value for setting the seed to ensure reproducibility (default is 1984).  

@return A downsampled Seurat object with the specified cell level and a similar number of cells per final label.  

```{r}
# number of cells in each subtype before downsampling
table(scCD4$final_labels)

pbmc_core_ds = downsample_core(scCD4, target_count = NULL, use_median = T, labels_to_exclude = c("unknown"), seed = 123)

# number of cells in each subtype after downsampling
table(pbmc_core_ds$final_labels)
```

### PCA on downsampled cells
**scODIN function #5 - core_pca()**  
Wrapper for Seurat functions: FindVariableFeatures, ScaleData, and RunPCA.  
This function sequentially applies FindVariableFeatures, ScaleData, and RunPCA on a given Seurat object, with user-defined input for parameters such as the number of variable features, number of principal components, etc.  

@param `seurat_obj` A Seurat object on which the functions will be applied.  
@param `assay` Assay to use (default is `"RNA"`).  
@param `selection_method` Method for selecting variable features in FindVariableFeatures (default is `"vst"`).  
@param `nfeatures` Number of variable features to select in FindVariableFeatures (default is 2000).  
@param `verbose` Logical indicating whether to print progress (default is TRUE).  
@param `npcs` Number of principal components to compute in RunPCA (default is 50).  
@param `reduction_name` Name of the dimensional reduction set to store (default is `"core_PCA"`).  

@return A Seurat object with variable features found, data scaled, and PCA run.  

```{r}
pbmc_core_ds <- core_pca(pbmc_core_ds, 
                         assay = "RNA", 
                         selection_method = "vst", 
                         nfeatures = 2000, 
                         npcs = 50, 
                         verbose = TRUE, 
                         reduction_name = "core_PCA")
```

### kNN search for identities of unknown cell types
**scODIN function #6 - odin_knn()**  
k-Nearest Neighbors (kNN) search using core cells as reference and unknown cells as query.  
This function performs a k-Nearest Neighbors (kNN) search on a Seurat object to classify undefined cells using a reference dataset of core cell types. It wraps around Seurat's `FindTransferAnchors` and `TransferData` functions with additional parameters.  

@param `scData_sub` A Seurat object subset containing the cells to analyze. For instance subsetted to CD4 T cells.  
@param `reference_obj` A Seurat object containing the reference data (e.g., downsampled core cell types).  
@param `knn_method` The kNN method to use, either "annoy" or "rann" (default is "annoy").  
@param `npcs` Number of principal components to use (default is 50).  
@param `dims` Dimensions to use for kNN search (default is `1:50`).  
@param `normalization_method` The normalization method for the Seurat object (default is `"LogNormalize"`).  
@param `verbose` Logical flag to indicate if progress should be printed (default is TRUE).  
@param `n_trees` Number of trees to use for the kNN method (default is 50 for annoy).  
@param `k_weight` The number of neighbors to use when predicting labels (default is 50).  
@param `k_filter` Whether to filter anchors by distance (default is `NA`).  
@param `mapping_score_k` Logical, whether to compute the mapping score (default is TRUE).  

@return A Seurat object with transferred labels and prediction scores stored in the "predictions" assay.

```{r}
pbmc_undefined_transf <- odin_knn(scData_sub = scCD4, 
                            reference_obj = pbmc_core_ds, 
                            knn_method = "annoy", 
                            npcs = 50, 
                            dims = 1:50
                            )
```

### Filter out low scoring predicted cells
**scODIN function #7 - apply_score_filter()**  
Apply a score filter to kNN predictions.  
This function applies a score filter to the predictions in a specified assay of a Seurat object. It assigns "Unassigned" to any cell where the maximum prediction score is below the specified threshold and returns the name of the cell type with the highest score for those above the threshold.  

@param `scData` A Seurat object containing the predictions.  
@param `assay` A character string specifying the assay from which to retrieve the predictions. Defaults to "predictions".  
@param `slot` A character string specifying the slot from which to retrieve the data within the assay. Defaults to "data".  
@param `score.filter` A numeric value indicating the threshold below which cells are labeled as "Unassigned". Defaults to 0.5. This should be adjusted empirically, as few cell types in the core requires a higher cutoff than a core with many cell types.  

@return A character vector of predicted cell types, with "Unassigned" for cells that do not meet the score threshold.  

```{r}
# number of cells before score filter. Every cell is given an ID based on which cell has the highest score
table(pbmc_undefined_transf$predicted.celltype.l1)

pbmc_undefined_transf$predicted_id_knn = apply_score_filter(
  pbmc_undefined_transf,
  assay = "predictions",
  slot = "data",
  score.filter = 0.2
)

# number of cells after score filter.
table(pbmc_undefined_transf$predicted_id_knn)
```

### Combining core and predicted cells to one Seurat object
**scODIN function #8 - odin_merge()**  
Merge core cell type Seurat object and kNN predicted Seurat object.  
This function merges a core Seurat object with a kNN Seurat object, adds a source column to indicate the origin of each cell, and corrects the final labels by replacing dashes with underscores based on the original labels.  

@param `core_object` A Seurat object containing core cell types with original labels.  
@param `knn_object` A Seurat object containing kNN results with predicted labels.  
 
@return A merged Seurat object with updated final labels and a source column indicating whether the cell originated from the core or the kNN object.

```{r}
merged_seurat <- odin_merge(core_object = scCD4, knn_object = pbmc_undefined_transf)

# check that dimensions of original Seurat object fits with the recombined one
dim(scCD4)
dim(merged_seurat)
```

```{r}
# compare labels before and after kNN
table(scCD4$final_labels)
table(merged_seurat$final_labels)
```


## License

This project is licensed under the MIT License - see the LICENSE.md file for details.  
