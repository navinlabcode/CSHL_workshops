# Introduction to VisiumHD
?

## the difference of VisiumHD and Visium
?

2um resolution is subcellular

<img width="371" alt="image" src="https://github.com/user-attachments/assets/17f3f711-3a47-4cf1-a8f4-e5e62ea659bb" />

## how to choose the right bin size

?
![image](https://github.com/user-attachments/assets/29d73e75-14d7-4105-ac67-7042e2aa2d96)


## Data Downloads

Before starting this tutorial, please making sure to downloaded/find VisiumHD output.

In this tutorial, we will use one of 10X public datasets, **Human Breast Cancer (Fresh Frozen)**. 
More different tissue type data can be found at [here](https://www.10xgenomics.com/datasets?configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&query=HD)

``` bash
# standard spaceranger output for visiumHD
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.2/Visium_HD_FF_Human_Breast_Cancer/Visium_HD_FF_Human_Breast_Cancer_web_summary.html
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.2/Visium_HD_FF_Human_Breast_Cancer/Visium_HD_FF_Human_Breast_Cancer_cloupe_008um.cloupe
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.2/Visium_HD_FF_Human_Breast_Cancer/Visium_HD_FF_Human_Breast_Cancer_feature_slice.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.2/Visium_HD_FF_Human_Breast_Cancer/Visium_HD_FF_Human_Breast_Cancer_metrics_summary.csv
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.2/Visium_HD_FF_Human_Breast_Cancer/Visium_HD_FF_Human_Breast_Cancer_molecule_info.h5
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.1.2/Visium_HD_FF_Human_Breast_Cancer/Visium_HD_FF_Human_Breast_Cancer_spatial.tar.gz
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/3.1.2/Visium_HD_FF_Human_Breast_Cancer/Visium_HD_FF_Human_Breast_Cancer_binned_outputs.tar.gz

## unzip
tar -xvf Visium_HD_FF_Human_Breast_Cancer_spatial.tar.gz
tar -xvf Visium_HD_FF_Human_Breast_Cancer_binned_outputs.tar.gz
```
Just a quick chekcing without coding, check **cloupe** with [Loupe Browser 8.1.2](https://www.10xgenomics.com/support/software/loupe-browser/downloads)

## Check 10X web_summary

Always the easiest and first step to understand your sample with XXXX_web_summary.html

![image](https://github.com/user-attachments/assets/2c8b6554-0667-4011-88dc-b7680edbd902)


## That's get started!

### loading R packages

<details open>
  <summary>Click to expand R pacakges list</summary>
  
```r
{
  library(matrixStats)  ## remotes::install_version("matrixStats", version="1.1.0") ## packageVersion('matrixStats')
  library(scCustomize)
  
  library(Seurat) # packageVersion('Seurat'), 5.2.0
  
  library(fields) 
  library(grid)
  library(ggplot2)
  library(schard)
  library(cowplot)
  library(paletteer)
  library(SeuratExtend)
  library(reshape2)
  library(patchwork) 
  library(plotly)
  library(RColorBrewer)
  library(arrow)
  library(magick) 
  library(ggpubr)
  
  library(kableExtra)
  library(rvest)
  
  library(ggrepel) 
  
  ## visiumHD annotation
  library(spacexr)
  library(Rfast)
  
  ## copykat
  library(copykat)

  set.seed(1234) ## for Reproducible
}
```
</details>

Note, this tutorial was prepared using  ***Seurat 5***

you might need to install the packages if you are running your analysis on a local computer

Most of R packages can be installed using:
```r
## from CRAN repositories
install.packages('XXX')

##alternative ways 
BiocManager::install('XXX')
remotes::install_github('XXX')

## Install multiple packages
install.packages(c("P1","P2","P3",dependencies = TRUE)

## Install packages with a specific version using: 
remotes::install_version('XXX',version=1.2-3)
```

### Functions in R

Here is some pre-defined Functions will be used in this tutorial, which help to organize, reuse and simplify your code.

<details open>
  <summary>Click to expand R function code</summary>
  
```r

Fitler_cluster <- function(Object_8um){
## Extract metadata with cluster info
before <- Object_8um@meta.data %>%
  count(cluster = default_10X_cluster, name = "n_before")

after <- Object_8um_filtered@meta.data %>%
  count(cluster = default_10X_cluster, name = "n_after")

## Join and calculate filtered-out percentage
filtered_stats_8bin <- left_join(before, after, by = "cluster") %>%
  mutate(
    n_after = ifelse(is.na(n_after), 0, n_after),
    n_filtered = n_before - n_after,
    percent_filtered = round(100 * n_filtered / n_before, 1),
    percent_Keep = 100-percent_filtered
  )

if( any(is.na(filtered_stats_8bin$cluster)) ){
filtered_stats_8bin <- filtered_stats_8bin[-which(is.na(filtered_stats_8bin$cluster)),]
}
 
p2_8u_BF <- ggplot(filtered_stats_8bin, aes(x = factor(cluster), y = percent_Keep)) +
  geom_bar(stat = "identity", fill = "#3182bd", width = 0.7) +
  geom_text(
    aes(label = n_after),
    vjust = -0.6,
    size = 4.5,
    fontface = "bold"
  ) +
  labs(
    title = "Qualified bins per Cluster",
    x = "Cluster",
    y = "Percentage of bins Pass QC"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text( size = 16, hjust = 0.5),
    axis.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

return(p2_8u_BF)
}

PieChart_Spot <- function(objectPie){
    # objectPie <-Object_8um_filteredQ
    spot_counts <- as.data.frame(table(objectPie$spot_class))  # Count cell types
    colnames(spot_counts) <- c("spot_class","count")
    spot_counts$fraction = round(spot_counts$count / sum(spot_counts$count) * 100,1)
    spot_counts$ymax = cumsum(spot_counts$fraction)
    spot_counts$ymin = c(0, head(spot_counts$ymax, n=-1))
    spot_counts$labelPosition <- (spot_counts$ymax + spot_counts$ymin) / 2
    spot_counts$label <- paste0(spot_counts$spot_class, "\n ", 
                            spot_counts$count," (",spot_counts$fraction,"%)")
  
  ## label only on the legend
  set1_colors <- brewer.pal(n = length(unique(spot_counts$spot_class)), name = "Set1")
  spot_counts$label_legend <- paste0(spot_counts$spot_class," (",spot_counts$fraction,"%)")
      
  Bin8_RCTD_spot_Pie <-plot_ly(
  data = spot_counts,
  labels = ~label_legend,
  values = ~count,
  type = 'pie',
  textinfo = 'none',  # hide labels on chart
  #hoverinfo = 'label+value+percent',  # keep tooltip
  hoverinfo = 'label+value',  # keep tooltip
  marker = list(colors = set1_colors)
  ) %>%
  layout(
    showlegend = TRUE,  # show the legend
    legend = list(orientation = "v", x = 1.05, y = 1)  # position legend to right
  )
    
    return(Bin8_RCTD_spot_Pie)
}
  
## pie chart for cell type
PieChart_CellType <- function(objectPie){
    #objectPie<-Object_8um_filteredQ
    spot_counts <- as.data.frame(table(objectPie$RCTD_CellType))  # Count cell types
    colnames(spot_counts) <- c("RCTD_CellType","count")
    spot_counts$fraction = round(spot_counts$count / sum(spot_counts$count) * 100,1)
    spot_counts$ymax = cumsum(spot_counts$fraction)
    spot_counts$ymin = c(0, head(spot_counts$ymax, n=-1))
    spot_counts$labelPosition <- (spot_counts$ymax + spot_counts$ymin) / 2
    spot_counts$label <- paste0(spot_counts$RCTD_CellType, "\n ", 
                            spot_counts$count," (",spot_counts$fraction,"%)")
    
    spot_counts$label_legend <- paste0(spot_counts$RCTD_CellType," (",spot_counts$fraction,"%)")
 
   ## label only on the legend
   Bin8_RCTD_CellType_Pie <-plot_ly(
  data = spot_counts,
  labels = ~label_legend,
  values = ~count,
  type = 'pie',
  textinfo = 'none',  # hide labels on chart
  #hoverinfo = 'label+value+percent',  # keep tooltip
  hoverinfo = 'label+value',  # keep tooltip
  marker = list(colors = DNA_col[seq_len(nrow(spot_counts))])
  ) %>%
  layout(
    showlegend = TRUE,  # show the legend
    legend = list(orientation = "v", x = 1.05, y = 1)  # position legend to right
  )
   
    return(Bin8_RCTD_CellType_Pie)
}

## clone color 
DNA_col <- c("#ffe5cd","#F70373","#822E1CFF","#AA0DFEFF","#1CBE4FFF",
             "#ff89cc","#1CFFCEFF","#ed968c","#f35d36","#bbb1ff","#D3FC5E",
             "#7c4b73","#CD69C9","#31c7ba","#5f5c0b","#4d1da9","#f4d451","#d0ffb7","#239f6e","#1b80ad","#2F4F4F",
             "#acb1b4","#080000","#2ab2e5","#b97e45","#f03c2b","#daa9d9","#63ff55","#ebc2b9",
             "#fae70f","#c9ceda","#564c6c","#4539dd","#dd0cc5","#c6662f","#105c13","#dd7d6d","#b1d8ff","#FEAF16FF",
             "#ffd000","#6596cd","#b90303","#aabf88","#534e46","#974949","#828282","#bd8399","#5373a7")
 ```
</details>

### Create Seurat spatial Objects

```r
Object <- Load10X_Spatial(data.dir = params$data_path, bin.size = c(8, 16))
```

![image](https://github.com/user-attachments/assets/5323e02d-9b3d-4d89-92a8-cc95979872a9)


### load 10X default graph cluster

we only load and domenstrate 8um here, feel free to try other bin size.

```r
## split into 8 um and load 10X umap and cluster
DefaultAssay(Object) <- "Spatial.008um"
Object_8um <- DietSeurat(Object,
                       assays = "Spatial.008um")
projection_8um <- read.csv(paste0(params$data_path,"/binned_outputs/square_008um/analysis/clustering/gene_expression_graphclust/clusters.csv",sep=""), row.names = 1,stringsAsFactors=F,check.names = FALSE)
cluster_8um <- read.csv(paste0(params$data_path,"/binned_outputs/square_008um/analysis/umap/gene_expression_2_components/projection.csv"), row.names = 1,stringsAsFactors=F,check.names = FALSE)
projection_cluster_8um <- cbind(projection_8um,cluster_8um)
#head(projection_cluster_8um)

Object_8um[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(projection_cluster_8um[, c("UMAP-1", "UMAP-2")]), 
                                             key = "UMAP-", assay = DefaultAssay(Object_8um))
Object_8um$default_10X_cluster <- projection_cluster_8um$Cluster[match(colnames(Object_8um), rownames(projection_cluster_8um))]

p1_8u <- DimPlot(Object_8um, reduction = "umap", 
        group.by = c("default_10X_cluster"),
        ncol = 1,label=T,raster = T)

p2_8u <-SpatialDimPlot(Object_8um, 
               label = F, repel = F,  
               pt.size.factor = 1.5, 
               image.alpha = 0.6,
               stroke = 0.01,
               group.by = c("default_10X_cluster"))


plot_grid(p1_8u,p2_8u,ncol=1)
```

![image](https://github.com/user-attachments/assets/61aa92f6-382a-415b-a417-a282c7a72884)


### QC metrics Visualization in 8 um Bin

```r
Object_8um[["percent.mt"]] <- PercentageFeatureSet(object = Object_8um, pattern = "^MT-")


bin8_bin_ncount <- ggplot(Object_8um@meta.data, aes(x = nCount_Spatial.008um)) +
  xlim(0,500) + 
  #ylim(0,100000) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "black") +
  labs(x = "nCount_Spatial.008um",
       y = "Cell Count") +
  geom_vline(xintercept = 25, color = "red", linetype = "dashed", size = 0.8) +
  theme_minimal()

bin8_bin_nFeature <- ggplot(Object_8um@meta.data, aes(x = nFeature_Spatial.008um)) +
  xlim(0,500) + 
  #ylim(0,100000) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "black") +
  labs(x = "nFeature_Spatial.008um",
       y = "Cell Count") +
  geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 0.8) +
  theme_minimal()

Bin8_cellViability4 <-ggplot(Object_8um@meta.data,aes(x=percent.mt))+stat_ecdf(aes(colour=orig.ident))+ scale_x_continuous(name = "percent.Mitochondrial",breaks=seq(0,100,10),limits=c(0,100))+theme_bw()+ylab("The percentage of bins")+ theme(legend.position="none")

Bin8_QC <-plot_grid(bin8_bin_nFeature,bin8_bin_ncount,Bin8_cellViability4,ncol=3)

Bin8_Feature <- FeaturePlot(Object_8um,features = c("nFeature_Spatial.008um","nCount_Spatial.008um")) + theme(legend.position = "right")
Bin8_S_Feature <- SpatialFeaturePlot(Object_8um, features = c("nFeature_Spatial.008um","nCount_Spatial.008um")) + theme(legend.position = "top")

plot_grid(Bin8_QC,Bin8_Feature,Bin8_S_Feature,ncol=1)
```
![image](https://github.com/user-attachments/assets/7c67e86e-054b-4962-94d6-8b28bcde67f6)


### Filtering criterion for 8um Bin

1. Bins detected in fewer than 3 spots were excluded  
2. Bins with fewer than 25 total detected genes were excluded

```r
Object_8um_filtered <- subset(
  Object_8um,
  subset = nFeature_Spatial.008um > 3 &
           nCount_Spatial.008um > 25) #& percent.mt < 20)

Bin8_pass <- (round(ncol(Object_8um_filtered) / ncol(Object_8um)*100,2))

table_8um <-data.frame(`8um Bin`=ncol(Object_8um_filtered),
                       Bin8_pass,
                       round(median(Object_8um_filtered@meta.data$nFeature_Spatial.008um)) )

colnames(table_8um)<-c("Number of Bin","Filter Ratio(%)","Median Genes per Bin")
table_8um
```
![image](https://github.com/user-attachments/assets/46e68af1-dd5c-4e27-94a0-d7baeb427977)

### checking the Filtering result

```r

## check the filter cell in each cluster 
p2_8u_BF <- Fitler_cluster(Object_8um)

Keep_cells <- colnames(Object_8um) %in% colnames(Object_8um_filtered)
Object_8um$QC_passed <- Keep_cells

p1_8u_F <- DimPlot(Object_8um, reduction = "umap", 
        group.by = c("QC_passed"),
        ncol = 1,label=T,raster = T)

p2_8u_F <-SpatialDimPlot(Object_8um, 
               label = F, repel = F,  
               pt.size.factor = 1.5, 
               image.alpha = 0.6,
               stroke = 0.01,
               group.by = c("QC_passed"))

plot_grid(p2_8u_BF, plot_grid(p1_8u_F,p2_8u_F,ncol=2),
          rel_heights=c(0.4,0.6),ncol=1)
```
![image](https://github.com/user-attachments/assets/7ce95c98-85dd-4384-90fe-b5738e18cf43)


### UMAP Cluster and Marker Genes for 8um Bin (10k bins subset)

```r

## subset for 10k cell
sampled_cells <- sample(Cells(Object_8um_filtered), size = 10000, replace = FALSE)
Object_8um_filteredQ <- subset(Object_8um_filtered, cells = sampled_cells)}

## standard seurat pre-processing like scRNA
{
Object_8um_filteredQ <- NormalizeData(Object_8um_filteredQ)
Object_8um_filteredQ <- FindVariableFeatures(Object_8um_filteredQ,nfeatures = 2000)
Object_8um_filteredQ <- ScaleData(Object_8um_filteredQ)

Object_8um_filteredQ <- RunPCA(Object_8um_filteredQ, verbose = FALSE)
Object_8um_filteredQ <- RunUMAP(Object_8um_filteredQ, dims = 1:30, verbose = FALSE)
Object_8um_filteredQ <- FindNeighbors(Object_8um_filteredQ, dims = 1:30, verbose = FALSE)
Object_8um_filteredQ <- FindClusters(Object_8um_filteredQ, verbose = FALSE,resolution = 0.1)
}

bin8_F_d <- DimPlot(Object_8um_filteredQ, label = TRUE)

bin8_F_sd <-SpatialDimPlot(Object_8um_filteredQ, 
               label = F, repel = F,  
               pt.size.factor = 10, 
              image.alpha = 0.5,
              stroke = 0.1)

Bin8_filter_VP <- VlnPlot(Object_8um_filteredQ, features = c("nCount_Spatial.008um"), pt.size = 0, log = TRUE)

plot_grid(bin8_F_d,bin8_F_sd,ncol = 2)

```
![image](https://github.com/user-attachments/assets/6a6de858-cde0-4b8c-ae16-246f151d45ae)


```r
## cluster markers
bin8_markers <- FindAllMarkers(Object_8um_filteredQ, only.pos = TRUE)
bin8_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> bin8_top10

bin8_F_G <- DoHeatmap(Object_8um_filteredQ, features = bin8_top10$gene, size = 2.5) + theme(axis.text = element_text(size = 8)) + NoLegend()
bin8_F_G
```

![image](https://github.com/user-attachments/assets/e12da33f-c402-42ec-aa20-5820ff0509c6)

Due to the probe-based chemistry (similar to FLEX) used in Visium HD, and the platform’s inherent characteristics (such as targeted transcript capture, limited probe design coverage, and reduced sensitivity to low-abundance transcripts), not all canonical scRNA-seq markers are reliably captured in Visium HD data.

```r
## canonical scRNA markers

Major_CellType_Markers  <- list(
  Endo_genes = c("CDH5","VWF","PECAM1","ESAM","THBD","CD36", "CD34", "HEG1"),
  Immun_genes = c("PTPRC","SRGN","TYROBP"),
  myeloid_genes = c("CD14", "CD33", "TREM1"),
  Epi_genes = c("EPCAM","KRT5","KRT17","KRT6","KRT18","KRT19", "KRT8", "KRT14"),
  Adipo_genes = c("APMAP","ADIPOQ","ADIPOQ-AS1","TPRA1"),
  Fibro_genes = c("COL1A1","LUM","FBLN1", "DCN", "FBN1", "COL1A2", "PCOLCE")
)

Dotplot_Group <- function(objectQ_subset,CellType_Markers) {
  marker_genes <- unlist(CellType_Markers)
  gene_annotation <- stack(marker_genes)
  colnames(gene_annotation) <- c("Gene", "Category")
  gene_annotation$Category <- gsub("_genes[0-9]", "", gene_annotation$Category)
  head(gene_annotation)
  
  dotplot <- DotPlot(objectQ_subset, features = as.character(marker_genes)) +
    coord_flip() +  # Rotate the plot
    theme_minimal()
  
  dotplot$data <- merge(dotplot$data, gene_annotation, by.x = "features.plot", by.y = "Gene", all.x = TRUE)
  
  # Re-plot with facets for annotation
  P <- dotplot + facet_grid(rows = vars(Category), scales = "free_y", space = "free_y")
  
  return(P)
}

bin8_F_C <- Dotplot_Group(Object_8um_filteredQ,Major_CellType_Markers)
bin8_F_C
```
![image](https://github.com/user-attachments/assets/ca96f2f6-b44e-44db-9585-fe6f2036e137)

## Spot Decomposition with RCTD

[Robust Cell Type Decomposition (RCTD)](https://www.nature.com/articles/s41587-021-00830-w) is a probabilistic method used to infer cell type composition in spatial transcriptomics data by mapping expression profiles to a single-cell RNA-seq reference. Each spatial bin is classified based on how confidently the model matches the expression to one or more cell types.

For the current human breast tissue VisiumHD data, we will utilize a single-cell reference dataset sourced from: [Kumar, T., Nee, K., Wei, R. et al. A spatially resolved single-cell genomic atlas of the adult human breast. Nature 620, 181–191 (2023)](https://doi.org/10.1038/s41586-023-06252-9)

```r
## load the RCTD reference made from Navin lab HBCA scRNA data
reference_path <- paste("THE_PATH_TO_Navin_HBCA_rds")
reference <- readRDS(reference_path)

## prepare the query 
counts_hd <- GetAssayData(Object_8um_filteredQ, slot = "counts")
cortex_cells_hd <- colnames(Object_8um_filteredQ)
coords <- GetTissueCoordinates(Object_8um_filteredQ)[cortex_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
  
## run RCTD (it take around 4 hours to run on my server)
Bin8_RCTD <- create.RCTD(query, reference, max_cores = 1, UMI_min = 10) ## the number of core need to set as 1 on core server, but it take all cores anyway
Bin8_RCTD <- run.RCTD(Bin8_RCTD, doublet_mode = "doublet")

## RCTD Result
Object_8um_filteredQ <- AddMetaData(Object_8um_filteredQ, metadata = Bin8_RCTD@results$results_df)
Object_8um_filteredQ$RCTD_CellType <- as.character(Object_8um_filteredQ$first_type)
Object_8um_filteredQ$RCTD_CellType[is.na(Object_8um_filteredQ$RCTD_CellType)] <- "Unknown"

## Skip waiting time by load pre-run rds object
Object_8um_filteredQ <- readRDS( paste(wd, params$SampleName,"_8um_Filtered_Subset10k.rds",sep=""))
```

### Spot Decomposition with RCTD: bin classifications

RCTD fits a probabilistic Bayesian model at each spatial location (spot or bin) to determine whether the observed gene expression is best explained. 
It calculates posterior probabilities for each candidate cell type (or pair of types).

Classification Summary:
- **Singlet**: The bin is confidently assigned to a single cell type when one cell type has high posterior probability (e.g., >0.8 or 0.9).
- **Doublet_certain**: Two distinct cell types are confidently predicted with balanced contributions.  
- **Doublet_uncertain**: Two cell types are detected, but the second has a lower or uncertain contribution.
- **Reject**: No reliable cell type prediction due to low confidence or ambiguous expression profiles.

```r
## spot class
Idents(Object_8um_filteredQ) <- "spot_class"

Bin8_RCTD_spot_VP <- VlnPlot(Object_8um_filteredQ, features = c("nCount_Spatial.008um"), pt.size = 0, log = TRUE) +  scale_fill_brewer(palette="Set1")

Bin8_RCTD_dimplot_spot <- DimPlot(Object_8um_filteredQ, reduction = "umap",pt.size = 1 ,
        group.by = c("spot_class"), ncol = 1) + ggtitle("UMAP of RCTD clustering") + scale_fill_brewer(palette="Set1")


Bin8_RCTD_spot_DP <- SpatialDimPlot(Object_8um_filteredQ, label = F, repel = F,  
                      pt.size.factor = 10, 
                      image.alpha = 0.5,
                      stroke = 0.1) +  scale_fill_brewer(palette="Set1")

Bin8_RCTD_spot_Pie <- PieChart_Spot(Object_8um_filteredQ)

Bin8_RCTD_spot_ALL <- plot_grid(Bin8_RCTD_spot_Pie,
                                Bin8_RCTD_spot_VP,
                                Bin8_RCTD_dimplot_spot,
                                Bin8_RCTD_spot_DP, 
                                rel_heights = c(0.2,0.2,0.3,0.3), ncol=1)
Bin8_RCTD_spot_ALL
```
![image](https://github.com/user-attachments/assets/3e2a65b6-d7bf-418f-80ec-42a5cdb828a2)

![image](https://github.com/user-attachments/assets/03c7632f-a136-4864-9e72-7b3cbfec1d92)

### Spot Decomposition with RCTD: Cell Type

```r
## cell type
Idents(Object_8um_filteredQ) <- "RCTD_CellType"

## fix the cell type order
celltypes_8 <- as.character(Idents(Object_8um_filteredQ))
current_levels_8 <- unique(sort(celltypes_8))
unknown_levels_8 <- sort(current_levels_8[grepl("(?i)^unknown$", current_levels_8, perl = TRUE)])
other_levels_8 <- setdiff(current_levels_8, unknown_levels_8)
new_levels_8 <- c(other_levels_8, unknown_levels_8)
Object_8um_filteredQ$RCTD_CellType <- factor(celltypes_8, levels = new_levels_8)

Bin8_RCTD_celltype_VP <- VlnPlot(Object_8um_filteredQ, features = c("nCount_Spatial.008um"), group.by = "RCTD_CellType",pt.size = 0, log = TRUE, cols = DNA_col) + NoLegend()

Bin8_RCTD_dimplot_celltype <- DimPlot(Object_8um_filteredQ, reduction = "umap",pt.size = 1 , group.by = c("RCTD_CellType"), cols = DNA_col, ncol = 1) + ggtitle("UMAP of RCTD clustering")

## cluster markers
bin8_celltye_markers <- FindAllMarkers(Object_8um_filteredQ, only.pos = TRUE)
bin8_celltye_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> bin8_ct_top5

bin8_CT_Heat <- DoHeatmap(Object_8um_filteredQ, features = bin8_ct_top5$gene, size = 2.5) + theme(axis.text = element_text(size = 8)) + NoLegend()

Bin8_RCTD_CellType_DP <- SpatialDimPlot(Object_8um_filteredQ, label = F, repel = F,  
                                        group.by = "RCTD_CellType",
                                        pt.size.factor = 10, 
                                        image.alpha = 0.5,
                                        stroke = 0.1)  + 
                                        scale_fill_manual(values = DNA_col)

Bin8_RCTD_CellType_Pie <- PieChart_CellType(Object_8um_filteredQ)

Bin8_RCTD_CellType_ALL <- plot_grid(Bin8_RCTD_CellType_Pie,
                                    Bin8_RCTD_celltype_VP,
                                    bin8_CT_Heat,
                                    Bin8_RCTD_dimplot_celltype,
                                    Bin8_RCTD_CellType_DP, 
                                    rel_heights = c(0.2,0.2,0.3,0.3,0.3), ncol=1)
Bin8_RCTD_CellType_ALL
```
![image](https://github.com/user-attachments/assets/d5b3a6b5-b9ad-49be-86ba-6737a953f907)

![image](https://github.com/user-attachments/assets/f287cd39-84d0-4648-b5d5-bdda50fe7049)

## copykat Infering CNV

Copykat is the Navin lab tool that infers genome-wide copy number variations (CNVs) from single-cell RNA-seq data, distinguishing aneuploid tumor cells from diploid normal cells.
Please refer to our published paper [Delineating copy number and clonal substructure in human tumors from single cell transcriptomes](https://www.nature.com/articles/s41587-020-00795-2)

```r

SampleName <- "BreastCancer_10X_FF_copykat_8um"

## run copykat
objectQS_cpk <- copykat(rawmat=Object_8um_filteredQ@assays$Spatial.008um$counts,
        ngene.chr=5,
        win.size=25,
        KS.cut=0.2, ##  segmentation parameter 0.05-0.15, Increasing KS.cut decreases sensitivity, i.e. less segments/breakpoints
        min.gene.per.cell = 100, ## mini
        sam.name=SampleName, distance="euclidean", n.cores=50)

knitr::include_graphics(paste0(SampleName,"_copykat_heatmap.jpeg"))
```
![image](https://github.com/user-attachments/assets/bdda191f-0654-4b8b-a31e-1b7703de9242)

### Cell Ploidy Mapping

```r
ck <- read.table(paste(SampleName,"_copykat_prediction.txt",sep=""),header = T)

Object_8um_filteredQ@meta.data$CopyKat<-"Undefined"
diploid<-ck[which((ck$copykat.pred=="diploid")|(ck$copykat.pred=="c1:low.confidence")),]
aneuploid<-ck[which((ck$copykat.pred=="aneuploid")|(ck$copykat.pred=="c2:low.confidence")),]

Object_8um_filteredQ@meta.data[which(row.names(Object_8um_filteredQ@meta.data) %in% aneuploid$cell.names),"CopyKat"]<-"Aneuploid"
Object_8um_filteredQ@meta.data[which(row.names(Object_8um_filteredQ@meta.data) %in% diploid$cell.names),"CopyKat"]<-"Diploid"

saveRDS(Object_8um_filteredQ, paste(params$SampleName,"_8um_Filtered_Subset10k.rds",sep=""))

Bin8_dimplot_copykat <- DimPlot(Object_8um_filteredQ, reduction = "umap",pt.size = 1 , group.by = c("CopyKat"), ncol = 1) + ggtitle("UMAP of RCTD clustering")

Bin8_Sdimplot_copykat <- SpatialDimPlot(Object_8um_filteredQ, label = F, repel = F,  
                        pt.size.factor = 5, 
                        image.alpha = 0.5,
                        stroke = 0.1,
                        group.by = c("CopyKat") )

plot_grid(Bin8_dimplot_copykat, Bin8_Sdimplot_copykat, ncol=2)
```
![image](https://github.com/user-attachments/assets/a6beff9e-4903-48da-9c58-6b8226f4621b)

## Spatial ToolBox

| Software  | Language | Unique Key Feature | Link |
|-----------|----------|---------------------|------|
| **Seurat** | R | Flexible integration and analysis of single-cell and spatial transcriptomics data | [Link](https://satijalab.org/seurat/) |
| **semla** | R | Efficient spatial transcriptomics visualization and analysis tailored for Seurat objects | [Link](https://github.com/ludvigla/semla) |
| **Squidpy** | Python | Graph-based spatial analysis with image integration and neighborhood statistics | [Link](https://squidpy.readthedocs.io/en/stable/index.html) |
| **stLearn** | Python | Joint modeling of spatial distance, tissue morphology, and gene expression patterns | [Link](https://stlearn.readthedocs.io/en/latest/index.html) |
| **StarDist** | Python | Deep learning-based cell segmentation from microscopy images using star-convex polygons | [Link](https://github.com/stardist/stardist) |
| **bin2cell** | Python | Assigns VisiumHD gene expression bins to segmented cells (expansion with on StarDist result) for single-cell spatial resolution | [Link](https://github.com/yizhouh/cell2location/tree/master/bin2cell) |





