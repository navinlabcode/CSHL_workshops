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



