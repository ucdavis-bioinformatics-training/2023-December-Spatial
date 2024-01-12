---
title: "Advanced Topics in Single Cell RNA-Seq: Spatial Transcriptomics"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---



# Advanced Topics in Single Cell RNA-Seq: Spatial Transcriptomics

This documentation was produced following the workshop in December 2023 to remedy issues with the documentation used at the time. An accompanying recording of a documentation walk-through will be made available to December 2023 registrants.

## Workspace set-up

A number of libraries are required for proper function of this documentation including Seurat, Rfast2, and spacexr for the analysis, ggplot2, viridis, and patchwork for visualizations, and knitr for formatting.


```r
if (!any(rownames(installed.packages()) == "knitr")){
  BiocManager::install("knitr")
}
if (!any(rownames(installed.packages()) == "Seurat")){
  BiocManager::install("Seurat")
}
if (!any(rownames(installed.packages()) == "Rfast2")){
  BiocManager::install("Rfast2")
}
```

<div class='r_output'> 
 The downloaded binary packages are in
 	/var/folders/f6/xwh6hmj94v3bpnmr3yhs9vm80000gn/T//RtmpaxyOi2/downloaded_packages
</div>
```r
if (!any(rownames(installed.packages()) == "hdf5r")){
  BiocManager::install("hdf5r")
}
if (!any(rownames(installed.packages()) == "Matrix")){
  BiocManager::install("Matrix")
}
if (!any(rownames(installed.packages()) == "ggplot2")){
  BiocManager::install("ggplot2")
}
if (!any(rownames(installed.packages()) == "viridis")){
  BiocManager::install("viridis")
}
if (!any(rownames(installed.packages()) == "patchwork")){
  BiocManager::install("patchwork")
}
if(!any(rownames(installed.packages()) == "devtools")) {
  install.packages("devtools")
}
if(!any(rownames(installed.packages()) == "spacexr")) {
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
}
library(knitr)
library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)
library(spacexr)
```

Download the Space Ranger output from tadpole to your working directory with the code below. I suggest copying the commands one at a time into a working terminal (e.g. the "Terminal" tab of RStudio) instead of running the chunk. Don't forget to replace "username" with your username!


```bash
scp username@tadpole.genomecenter.ucdavis.edu:/share/workshop/Spatial_Workshop/Downloads/data.tar.gz .
tar -xzvf data.tar.gz
<div class='r_output'>
The data directory downloaded by this command contains one directory for each sample with the necessary output files from Space Ranger. To replicate this structure, you could simply rename the "outs" directory with the sample name. The design.tsv file used to run Space Ranger has been copied to the data directory as well, along with a reference expression data set we will make use of later in the analysis.


```r
set.seed(1234)
project.name <- "Mouse Brain Sagittal Sections"
dataset.loc <- "data/"
design <- read.delim(paste0(dataset.loc, "design.tsv"),
                     header = FALSE,
                     col.names = c("sample", "image", "slide", "capture.area"))

experiment.slices <- lapply(c(1:dim(design)[1]), function(i){
  # read in image file
  image = Read10X_Image(paste0(dataset.loc, design$sample[i], "/outs/spatial"))
  # create Seurat object with "Spatial" assay
  sce = Load10X_Spatial(paste0(dataset.loc, design$sample[i], "/outs"),
                        filename = "filtered_feature_bc_matrix.h5",
                        assay = "Spatial",
                        slice = design$sample[i],
                        filter.matrix = TRUE,
                        to.upper = FALSE,
                        image = image)
  sce$orig.ident = design$sample[i]
  # add slide serial number to metadata
  sce = AddMetaData(sce, metadata = design$slide[i], col.name = "slide")
  # add capture area to metadata
  sce = AddMetaData(sce, metadata = design$capture.area[i], col.name = "capture.area")
  # create shorter sample label for visualization
  sce = AddMetaData(sce, metadata = sub("Section_", "", sub("V1_Mouse_Brain_Sagittal_", "", design$sample[i])), col.name = "sample.ident")
   rename image to match
  names(sce@images) = sub("Section_", "", sub("V1_Mouse_Brain_Sagittal_", "", design$sample[i]))
  # append sample name to cell barcode
  RenameCells(sce, new.names = paste(sapply(strsplit(Cells(sce), split = "-"), "[[", 1), design$sample[i], sep = "-"))
  })

# merge matrices into a single object for QC, normalization, scaling, etc
experiment.merged <- merge(experiment.slices[[1]], experiment.slices[2:4])
experiment.merged <- JoinLayers(experiment.merged)

#tidy workspace
rm(experiment.slices)

experiment.merged
```

</div>## An object of class Seurat 
## 32285 features across 12162 samples within 1 assay 
## Active assay: Spatial (32285 features, 0 variable features)
##  1 layer present: counts
##  4 images present: Anterior_1, Posterior_1, Anterior_2, Posterior_2
<div class='r_output'>
 QA/QC

# Custom metrics

In many tissues, a high proportion of UMIs corresponding to mitochondrial transcripts is an indicator of poor viability. However, in an energy intensive tissue, high expression may be typical. In this case, a high percent mitochondrial reads is typical of the cells, making percent mitochondrial an inappropriate filtering criterion for this data.


```r
experiment.merged$percent.mito <- PercentageFeatureSet(experiment.merged,
                                                       pattern = "^mt-")
```

# Spatial plots

Indeed, overlaying the per-spot values for UMIs, gene count, and percent mitochondrial onto the slice images reveals a clear relationship between these features and the tissue structures.


```r
SpatialFeaturePlot(experiment.merged, "nCount_Spatial")
```

![](spatial_analysis_files/figure-html/per_spot_vis-1.png)<!-- -->

```r
SpatialFeaturePlot(experiment.merged, "nFeature_Spatial")
```

![](spatial_analysis_files/figure-html/per_spot_vis-2.png)<!-- -->

```r
SpatialFeaturePlot(experiment.merged, "percent.mito")
```

![](spatial_analysis_files/figure-html/per_spot_vis-3.png)<!-- -->

The legends and titles of figures can become compressed as the number of images Seurat is tiling increases. To view a single slice, use the "images" argument of the SpatialFeaturePlot function.


```r
SpatialFeaturePlot(experiment.merged, "percent.mito", images = c("Posterior_1", "Posterior_2"))
```

![](spatial_analysis_files/figure-html/per_spot_vis_single-1.png)<!-- -->

# Ridge plots


```r
RidgePlot(experiment.merged,
          features = "nCount_Spatial",
          group.by = "sample.ident") +
  guides(fill = "none") +
  scale_fill_viridis_d()
```

![](spatial_analysis_files/figure-html/RidgePlot-1.png)<!-- -->

```r
RidgePlot(experiment.merged,
          features = "nFeature_Spatial",
          group.by = "sample.ident") +
  guides(fill = "none") +
  scale_fill_viridis_d()
```

![](spatial_analysis_files/figure-html/RidgePlot-2.png)<!-- -->

```r
RidgePlot(experiment.merged,
          features = "percent.mito",
          group.by = "sample.ident") +
  guides(fill = "none") +
  scale_fill_viridis_d()
```

![](spatial_analysis_files/figure-html/RidgePlot-3.png)<!-- -->

# Scatter plots


```r
FeatureScatter(experiment.merged,
               feature1 = "nCount_Spatial",
               feature2 = "nFeature_Spatial",
               group.by = "sample.ident",
               shuffle = TRUE) +
  scale_color_viridis_d()
```

![](spatial_analysis_files/figure-html/FeatureScatter-1.png)<!-- -->

```r
FeatureScatter(experiment.merged,
               feature1 = "nCount_Spatial",
               feature2 = "percent.mito",
               group.by = "sample.ident",
               shuffle = TRUE) +
  scale_color_viridis_d()
```

![](spatial_analysis_files/figure-html/FeatureScatter-2.png)<!-- -->

# Filtering

This workflow uses the filtered_feature_bc_matrix.h5 file, which contains UMI counts per gene from each spot identified as underneath tissue by Space Ranger. The automatic tissue detection removes the majority of empty spots, however, visual artifacts or unsuccessful permeablization can result in spots with no (or very low) UMIs. The inclusion of zero-count spots interferes with normalization.

Unlike in single cell experiments, the presence of multiplets is not a concern; Visium experiments generally have multiple cells per spot, though the number varies from tissue to tissue.

In order to ensure that no empty spots are incorporated, we will remove any spots with fewer than 200 genes detected.


```r
experiment.merged
```

</div>## An object of class Seurat 
## 32285 features across 12162 samples within 1 assay 
## Active assay: Spatial (32285 features, 0 variable features)
##  1 layer present: counts
##  4 images present: Anterior_1, Posterior_1, Anterior_2, Posterior_2
<div class='r_output'>
```r
experiment.merged <- subset(experiment.merged, nFeature_Spatial >= 200)
experiment.merged
```

</div>## An object of class Seurat 
## 32285 features across 12146 samples within 1 assay 
## Active assay: Spatial (32285 features, 0 variable features)
##  1 layer present: counts
##  4 images present: Anterior_1, Posterior_1, Anterior_2, Posterior_2
<div class='r_output'>As only a small number of spots are removed, we will not reproduce the QC plots using the filtered data.

 Analysis

# Transform

```r
experiment.merged <- SCTransform(experiment.merged,
                                 assay = "Spatial",
                                 verbose = FALSE)
```

# PCA

```r
experiment.merged <- RunPCA(experiment.merged,
                            assay = "SCT",
                            npcs = 50,
                            verbose = FALSE)
```

# UMAP

```r
experiment.merged <- RunUMAP(experiment.merged,
                             reduction = "pca",
                             dims = 1:30,
                             verbose = FALSE)
DimPlot(experiment.merged,
        reduction = "umap",
        group.by = "sample.ident",
        shuffle = TRUE) +
  scale_color_viridis_d()
```

![](spatial_analysis_files/figure-html/UMAP-1.png)<!-- -->

```r
saveRDS(experiment.merged, file = "experiment.merged_UMAP.rds")
```

# Integrate

In this experiment, all four slices were located on the same slide, eliminating the most prominent source of batch variation. In larger Visium experiments, this may not be possible. When samples are spread across two or more slides, the effect of slide to slide variation can be mitigated with integration.

This code is provided as an example. For the remainder of the workshop, the un-integrated data will be used.


```r
experiment.harmony <- experiment.merged
DefaultAssay(experiment.harmony) <- "Spatial"
experiment.harmony@assays$SCT <- NULL
experiment.harmony@assays$Spatial <- split(experiment.harmony@assays$Spatial,
                                           f = experiment.harmony$sample.ident) # for real batch correction, use "slide"
experiment.harmony <- NormalizeData(experiment.harmony)
use.features <- rownames(experiment.merged)[rowSums(experiment.merged@assays$Spatial@layers$counts) >= 1]
experiment.harmony <- ScaleData(experiment.harmony, features = use.features)
experiment.harmony <- RunPCA(experiment.harmony, features = use.features)
experiment.harmony <- IntegrateLayers(object = experiment.harmony,
                                      method = "HarmonyIntegration",
                                      orig.reduction = "pca",
                                      new.reduction = "harmony",
                                      features = use.features,
                                      verbose = FALSE)
experiment.harmony <- RunUMAP(experiment.harmony,
                              reduction = "harmony",
                              dims = 1:50,
                              verbose = FALSE)
p1 <- DimPlot(experiment.harmony,
              reduction = "umap",
              group.by = "sample.ident",
              shuffle = TRUE) +
  scale_color_viridis_d() +
  ggtitle("Integrated")
p2 <- DimPlot(experiment.merged,
              reduction = "umap",
              group.by = "sample.ident",
              shuffle = TRUE) +
  scale_color_viridis_d() +
  ggtitle("Uncorrected")
p1 + p2
```

![](spatial_analysis_files/figure-html/harmony-1.png)<!-- -->

```r
rm(use.features, experiment.harmony, p1, p2)
```

Clustering and visualizations may be performed on the integrated data, but for differential expression analysis, the layers should be re-joined, and the uncorrected data used.

Running Harmony does not appear to have altered the UMAP much. Though it has changed, the overall relationship between points is quite similar before and after batch correction. For an example of Harmony results demonstrating a true batch correction, see [this 10x-provided vignette](https://www.10xgenomics.com/resources/analysis-guides/correcting-batch-effects-in-visium-data).

# Cluster


```r
experiment.merged <- FindNeighbors(experiment.merged,
                                   reduction = "pca",
                                   verbose = FALSE)
experiment.merged <- FindClusters(experiment.merged,
                                  resolution = seq(from = 0.15, to = 0.6, by = 0.15),
                                  verbose = FALSE)
lapply(grep("snn", colnames(experiment.merged@meta.data), value = TRUE),
       function(res){
         DimPlot(experiment.merged,
                 reduction = "umap",
                 group.by = res,
                 shuffle = TRUE) +
           scale_color_viridis_d(option = "turbo")
       })
```

</div>## [[1]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/cluster-1.png)<!-- -->

</div>## 
## [[2]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/cluster-2.png)<!-- -->

</div>## 
## [[3]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/cluster-3.png)<!-- -->

</div>## 
## [[4]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/cluster-4.png)<!-- -->

```r
# merged object visualization is too small, plus color scale change applies only to final image
SpatialDimPlot(experiment.merged,
               group.by = "SCT_snn_res.0.3") +
  scale_fill_viridis_d(option = "turbo")
```

![](spatial_analysis_files/figure-html/cluster-5.png)<!-- -->

```r
# create palette
cluster.palette <- viridis(length(levels(experiment.merged$SCT_snn_res.0.3)), option = "turbo")
names(cluster.palette) <- levels(experiment.merged$SCT_snn_res.0.3)
# lapply to plot each slice individually
p <- lapply(names(experiment.merged@images), function(image){
  SpatialDimPlot(experiment.merged,
                 group.by = "SCT_snn_res.0.3",
                 images = c(image)) +
    scale_fill_manual(values = cluster.palette) +
    ggtitle(paste(image))
})
p[[1]] + p[[3]]
```

![](spatial_analysis_files/figure-html/cluster-6.png)<!-- -->

```r
p[[2]] + p[[4]]
```

![](spatial_analysis_files/figure-html/cluster-7.png)<!-- -->

```r
rm(p)
# save object
saveRDS(experiment.merged, "experiment.merged_clusters.rds")
```

# Find spatially variable genes

Spatially variable features are genes for which the spatial coordinates (i.e. location within the tissue) of spots explain expression level. In a layered tissue like the mouse brain samples we are using, spatial coordinates correspond closely to tissue structures and cell types. In other tissues, this may not be the case.

FindSpatiallyVariableFeatures should only be run on spots belonging to the same slice, as calculating spatial variability across discontinuous spots will produce artifacts and errors. Expression profiles vary widely across the mouse brain sagittal slices used in this experiment. To identify variation of interest (and run the algorithm in a timely manner), we can select one or more groups of spots corresponding to anatomical structures to interrogate for spatially variable features.

Here we select cluster 3, which appears to correspond roughly to the striatum (in blue, below).

![Allen Mouse Brain Atlas (Sagittal), image 13 of 21 zoom: 12.5%](figures/Sagittal_striatum_image13_zoom12.5.jpg)


```r
# subset object to slice and cluster of interest
anterior1.cluster3 <- subset(experiment.merged, sample.ident == "Anterior_1" & SCT_snn_res.0.3 == "3")
# remove extra images
anterior1.cluster3@images[2:4] <- NULL
# find spatially variable features
anterior1.cluster3 <- FindSpatiallyVariableFeatures(anterior1.cluster3,
                                                    assay = "SCT",
                                                    slot = "scale.data",
                                                    features = VariableFeatures(anterior1.cluster3),
                                                    selection.method = "markvariogram")
```

The SVFInfo and SpatiallyVariableFeatures functions, designed to access the slot modified by FindSpatiallyVariableFeatures, are producing errors at the time of writing this documentation, but addressing the slot directly produces results. There are [a number of open issues](https://github.com/satijalab/seurat/issues?q=is%3Aissue+is%3Aopen+findspatiallyvariablefeatures) in the Seurat 5 GitHub repository pertaining to FindSpatiallyVariableFeatures, so I expect that future updates will address the problem.


```r
svf.rank <- anterior1.cluster3@assays$SCT@meta.features
svf.rank <- svf.rank[order(svf.rank$markvariogram.spatially.variable.rank),]
svf.rank <- svf.rank[1:length(which(svf.rank$markvariogram.spatially.variable)),]
svf.rank[1:10,]
```

</div>##         r.metric.5 markvariogram.spatially.variable
## Ttr     0.03686007                             TRUE
## Enpp2   0.13260687                             TRUE
## Tnnt1   0.21006631                             TRUE
## Ecrg4   0.29415873                             TRUE
## Nrgn    0.31975442                             TRUE
## Slc17a7 0.33873826                             TRUE
## Kl      0.35503928                             TRUE
## Rasgrf2 0.36069930                             TRUE
## Necab3  0.36878837                             TRUE
## Igfbp2  0.37085088                             TRUE
##         markvariogram.spatially.variable.rank
## Ttr                                         1
## Enpp2                                       2
## Tnnt1                                       3
## Ecrg4                                       4
## Nrgn                                        5
## Slc17a7                                     6
## Kl                                          7
## Rasgrf2                                     8
## Necab3                                      9
## Igfbp2                                     10
<div class='r_output'>
Visualizing the expression with SpatialFeaturePlot can provide a dramatic illustration of variation across the selected region.

```r
lapply(rownames(svf.rank)[c(1,3,5,8)], function(feature){
    SpatialFeaturePlot(anterior1.cluster3, features = feature, crop = FALSE)
})
```

</div>## [[1]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/vis_top_spatial_cluster3-1.png)<!-- -->

</div>## 
## [[2]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/vis_top_spatial_cluster3-2.png)<!-- -->

</div>## 
## [[3]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/vis_top_spatial_cluster3-3.png)<!-- -->

</div>## 
## [[4]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/vis_top_spatial_cluster3-4.png)<!-- -->

However, expanding the view to include the entire tissue section reveals that some of the variation appears to be primarily influenced by expression in neighboring regions.

```r
lapply(rownames(svf.rank)[c(1,3,5,8)], function(feature){
  SpatialFeaturePlot(experiment.merged,
                     features = feature,
                     images = c("Anterior_1", "Anterior_2"))
})
```

</div>## [[1]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/vis_top_spatial_all-1.png)<!-- -->

</div>## 
## [[2]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/vis_top_spatial_all-2.png)<!-- -->

</div>## 
## [[3]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/vis_top_spatial_all-3.png)<!-- -->

</div>## 
## [[4]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/vis_top_spatial_all-4.png)<!-- -->

```r
rm(anterior1.cluster3)
```

# Cell type prediction and spot "decomposition"

Each spot on a Visium slide may incorporate one or more cell types, depending on the placement. Under the previous version of Seurat (4.4.0), the recommended method for cell type detection was integration. Instead of attempting to assign a single cell type to each spot, the Seurat 4 integration method calculated a prediction score for each cell type based on the expression profiles from a single cell data set. This methodology is still available (in version 5), but the authors now recommend using [Robust Cell Type Decomposition](https://www.nature.com/articles/s41587-021-00830-w)(RCTD).

The spacexr library installed at the beginning of this session supports RCTD. Here we use the spacexr function SpatialRNA to create an a list of objects (one for each slice) containing counts data and centroid coordinates for each spot.


```r
experiment.slices <- experiment.merged
DefaultAssay(experiment.slices) <- "Spatial"
experiment.slices <- split(experiment.slices, f = experiment.slices$sample.ident)
queries <- lapply(seq_along(Layers(experiment.slices)), function(i){
  coords = GetTissueCoordinates(experiment.slices, image = names(experiment.slices@images)[i], which = "centroids")
  counts = GetAssayData(experiment.slices, assay = "Spatial", layer = Layers(experiment.slices)[i])
  SpatialRNA(coords = coords, counts = counts, nUMI = colSums(counts))
})
```

The Allen Brain Atlas has a comprehensive collection of [publicly available single cell data sets](https://portal.brain-map.org/atlases-and-data/rnaseq). For this experiment, we have selected a small number of cells corresponding to cell types present in the slices. This sub-setting was necessary in order to create a reference set that will run quickly with limited resources. The Seurat object containing the atlas was created for a previous workshop under Seurat 4 and must be updated to Seurat 5 prior to continuing with the workflow. The spacexr reference is then created from the atlas object.


```r
atlas <- readRDS(paste(dataset.loc, "allen_subset_atlas.rds", sep = "/"))
atlas <- UpdateSeuratObject(atlas)
# RCTD requires a minimum of 25 cells per type. drop types with fewer cells.
atlas <- atlas[,!(atlas$cell_type_accession_label %in% names(which(table(atlas$cell_type_accession_label) < 25)))]
# alias labels (most readable) contain some special characters. create conversion table
atlas.labels <- atlas@meta.data[!duplicated(atlas$cell_type_alias_label),c("cell_type_accession_label", "cell_type_alias_label")]
rownames(atlas.labels) <- NULL
reference <- Reference(counts = GetAssayData(atlas, assay = "RNA", layer = "counts"),
                       cell_types = as.factor(atlas$cell_type_accession_label),
                       nUMI = atlas$nCount_RNA)
```

Once both the query (Visium) and reference (atlas) objects have been created, we can run the decomposition itself. This is more computationally intensive than earlier steps and can be run on more than one core to improve speed. The code below


```r
RCTD.list <- lapply(queries, function(query){
  RCTD = create.RCTD(query, reference, max_cores = 8)
  RCTD = run.RCTD(RCTD, doublet_mode = "doublet")
})
```

</div>## 
## CS202106160_114 CS202106160_115  CS202106160_12 CS202106160_167 CS202106160_189 
##              36              32             163             119              86 
## CS202106160_200 CS202106160_204 CS202106160_228 CS202106160_238 CS202106160_287 
##              64             105             137              90             231 
## CS202106160_288 CS202106160_293 CS202106160_348 CS202106160_363 CS202106160_375 
##             245              82              39             376              28 
## CS202106160_376  CS202106160_82  CS202106160_98  CS202106160_99 
##              28              70              50              26 
## [1] "gather_results: finished 1000"
## [1] "gather_results: finished 2000"
## 
## CS202106160_114 CS202106160_115  CS202106160_12 CS202106160_167 CS202106160_189 
##              36              32             163             119              86 
## CS202106160_200 CS202106160_204 CS202106160_228 CS202106160_238 CS202106160_287 
##              64             105             137              90             231 
## CS202106160_288 CS202106160_293 CS202106160_348 CS202106160_363 CS202106160_375 
##             245              82              39             376              28 
## CS202106160_376  CS202106160_82  CS202106160_98  CS202106160_99 
##              28              70              50              26 
## [1] "gather_results: finished 1000"
## [1] "gather_results: finished 2000"
## [1] "gather_results: finished 3000"
## 
## CS202106160_114 CS202106160_115  CS202106160_12 CS202106160_167 CS202106160_189 
##              36              32             163             119              86 
## CS202106160_200 CS202106160_204 CS202106160_228 CS202106160_238 CS202106160_287 
##              64             105             137              90             231 
## CS202106160_288 CS202106160_293 CS202106160_348 CS202106160_363 CS202106160_375 
##             245              82              39             376              28 
## CS202106160_376  CS202106160_82  CS202106160_98  CS202106160_99 
##              28              70              50              26 
## [1] "gather_results: finished 1000"
## [1] "gather_results: finished 2000"
## 
## CS202106160_114 CS202106160_115  CS202106160_12 CS202106160_167 CS202106160_189 
##              36              32             163             119              86 
## CS202106160_200 CS202106160_204 CS202106160_228 CS202106160_238 CS202106160_287 
##              64             105             137              90             231 
## CS202106160_288 CS202106160_293 CS202106160_348 CS202106160_363 CS202106160_375 
##             245              82              39             376              28 
## CS202106160_376  CS202106160_82  CS202106160_98  CS202106160_99 
##              28              70              50              26 
## [1] "gather_results: finished 1000"
## [1] "gather_results: finished 2000"
## [1] "gather_results: finished 3000"
<div class='r_output'>
The complex RCTD object created contains cell type predictions within the results slot, in a data.frame called "results_df." To add the predictions to the experiment.merged object, we collect the relevant column from the data.frame and collapse it into a named vector before using the new vector as input to AddMetaData.


```r
# add results from each slice to metadata table
annotations.list <- lapply(RCTD.list, function(slice){
  annotations = atlas.labels$cell_type_alias_label[match(slice@results$results_df$first_type, atlas.labels$cell_type_accession_label)]
  names(annotations) = rownames(slice@results$results_df)
  annotations
})
experiment.merged <- AddMetaData(experiment.merged,
                                 unlist(annotations.list),
                                 col.name = "predicted.celltype")

length(is.na(experiment.merged$predicted.celltype))
```

</div>## [1] 12146
<div class='r_output'>
```r
table(experiment.merged$predicted.celltype)
```

</div>## 
##       114_Pvalb       115_Pvalb        12_Lamp5 167_L2/3 IT CTX 189_L4/5 IT CTX 
##            1357            1142              29            1277             361 
##   200_L5 IT CTX 204_L5/6 IT CTX   228_L6 IT CTX        238_Car3   287_L6 CT CTX 
##            1624             569             947              36              29 
##   288_L6 CT CTX   293_L6 CT CTX      348_CA1-do          363_DG       375_Oligo 
##             113              31             131             131            2149 
##       376_Astro          82_Sst          98_Sst          99_Sst 
##            1444             264              69             443
<div class='r_output'>
```r
# establish color palette
celltype.palette <- viridis(dim(atlas.labels)[1], option = "mako")
names(celltype.palette) <- atlas.labels$cell_type_alias_label
# Spatial plot
p <- lapply(names(experiment.merged@images), function(slice){
  SpatialDimPlot(experiment.merged,
                 group.by = "predicted.celltype",
                 images = c(slice)) +
    scale_fill_manual(values = celltype.palette)
})
p[[1]] + p[[3]]
```

![](spatial_analysis_files/figure-html/predictions-1.png)<!-- -->

```r
p[[2]] + p[[4]]
```

![](spatial_analysis_files/figure-html/predictions-2.png)<!-- -->

```r
# UMAP
DimPlot(experiment.merged,
        group.by = "predicted.celltype",
        shuffle = TRUE) +
  scale_color_manual(values = celltype.palette)
```

![](spatial_analysis_files/figure-html/predictions-3.png)<!-- -->

Let's take a few minutes to explore the results of the cell type prediction. The [Cell Type Knowledge Explorer](https://knowledge.brain-map.org/celltypes/CCN202002013) offers a look at the types described by the alias label. Click on a wedge in the plot to go to information about the cell type(s) highlighted.

![Cell Type Knowledge Explorer](figures/Cell_Type_Knowledge_Explorer.png)

Using ggplot2, we can generate a custom plot to illustrate the cell type composition of each cluster.


```r
colnames(experiment.merged@meta.data)
```

</div>##  [1] "orig.ident"         "nCount_Spatial"     "nFeature_Spatial"  
##  [4] "slide"              "capture.area"       "sample.ident"      
##  [7] "percent.mito"       "nCount_SCT"         "nFeature_SCT"      
## [10] "SCT_snn_res.0.15"   "SCT_snn_res.0.3"    "SCT_snn_res.0.45"  
## [13] "SCT_snn_res.0.6"    "seurat_clusters"    "predicted.celltype"
<div class='r_output'>
```r
ggplot(data = experiment.merged@meta.data, mapping = aes(x = SCT_snn_res.0.3, fill = predicted.celltype)) + geom_bar() + scale_fill_manual(values = celltype.palette) + theme_classic()
```

![](spatial_analysis_files/figure-html/cluster_composition-1.png)<!-- -->

Cluster 3, which we selected for the spatially variable features analysis, appears to be enriched in cell types 167_L2/3 IT CTX and 228_L6 IT CTX. Use the table function to get exact counts.

Does the expression of the top spatially variable genes from cluster three differ between the cell types?


```r
lapply(rownames(svf.rank)[c(1,3,5,8)], function(feature){
  VlnPlot(experiment.merged,
          group.by = "predicted.celltype",
          idents = "3",
          features = feature) +
    scale_fill_manual(values = celltype.palette)
})
```

</div>## [[1]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/tnnt1_vln_cluster3-1.png)<!-- -->

</div>## 
## [[2]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/tnnt1_vln_cluster3-2.png)<!-- -->

</div>## 
## [[3]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/tnnt1_vln_cluster3-3.png)<!-- -->

</div>## 
## [[4]]
<div class='r_output'>
![](spatial_analysis_files/figure-html/tnnt1_vln_cluster3-4.png)<!-- -->

Does the relationship between predicted cell type and Tnnt1 expression hold up outside of cluster 3?


```r
VlnPlot(experiment.merged,
        group.by = "predicted.celltype",
        features = "Tnnt1") +
  scale_fill_manual(values = celltype.palette)
```

![](spatial_analysis_files/figure-html/tnnt1_vln_all-1.png)<!-- -->


```r
rm(annotations, annotations.list, atlas, atlas.labels, experiment.slices, queries, reference, RCTD, RCTD.list)
```


# Build "niches"

Seurat 5 introduces the concept of "niches," which organize cells by both cell type composition and spatial adjacency. While the clustering performed earlier relies on a network constructed on transcriptional similarity (via the PCA dimensionality reduction), niche-building uses a network describing the cell type identities of each spot's *spatial* neighbors. Then, k-means clustering identifies spots that have similar neighbors as members of the same "spatial niche." Like other functions that use spatial coordinates, BuildNicheAssay() produces errors if performed on multiple slices. At the moment it is simplest to subset the object to a single slice before creating the field of view object. See [Issue #8126](https://github.com/satijalab/seurat/issues/8126) for details.


```r
anterior1 <- subset(experiment.merged, sample.ident == "Anterior_1")
anterior1@images[2:4] <- NULL
anterior1[["fov"]] <- CreateFOV(coords = GetTissueCoordinates(anterior1, which = "centroids"),
                                type = "centroids")
anterior1 <- BuildNicheAssay(object = anterior1, fov = "fov", group.by = "predicted.celltype")
celltype.plot <- ImageDimPlot(object = anterior1, group.by = "predicted.celltype", size = 1.5, cols = celltype.palette, dark.background = FALSE)
niche.plot <- ImageDimPlot(object = anterior1, group.by = "niches", size = 1.5, cols = viridis(length(unique(anterior1$niches))), dark.background = FALSE)
celltype.plot + niche.plot
```

![](spatial_analysis_files/figure-html/BuildNicheAssay-1.png)<!-- -->
As with clustering, the resolution of niche building can be tuned to suit the biology of the tissue, in this case using the neighbors.k and niches.k.

The niche classification, which incorporates both transcriptional information (through cell type assignment) and spatial information (adjacency), may do a better job of representing the functional and biological divisions present in a tissue slice that the clustering performed earlier.

 Further exploratory visualizations

In addition to the gene expression data, we can plot relationships between continuous metadata and cluster, niche, or cell type assignment.


```r
VlnPlot(experiment.merged,
        group.by = "predicted.celltype",
        features = "percent.mito") +
  scale_fill_manual(values = celltype.palette)
```

![](spatial_analysis_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
VlnPlot(experiment.merged,
        group.by = "SCT_snn_res.0.3",
        features = "percent.mito") +
  scale_fill_manual(values = cluster.palette)
```

![](spatial_analysis_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```r
ggplot(data = anterior1@meta.data, mapping = aes(x = niches, fill = predicted.celltype)) + geom_bar() + scale_fill_manual(values = celltype.palette) + theme_classic()
```

![](spatial_analysis_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

 Session information


```r
sessionInfo()
```

</div>## R version 4.3.1 (2023-06-16)
## Platform: aarch64-apple-darwin20 (64-bit)
## Running under: macOS Monterey 12.4
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Los_Angeles
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] spacexr_2.2.1      patchwork_1.1.3    viridis_0.6.4      viridisLite_0.4.2 
## [5] ggplot2_3.4.4      Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-2          
## [9] knitr_1.45        
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3     rstudioapi_0.15.0      jsonlite_1.8.8        
##   [4] magrittr_2.0.3         spatstat.utils_3.0-4   farver_2.1.1          
##   [7] rmarkdown_2.25         vctrs_0.6.5            ROCR_1.0-11           
##  [10] spatstat.explore_3.2-5 htmltools_0.5.7        sass_0.4.8            
##  [13] sctransform_0.4.1      parallelly_1.36.0      KernSmooth_2.23-22    
##  [16] bslib_0.6.1            htmlwidgets_1.6.4      ica_1.0-3             
##  [19] plyr_1.8.9             plotly_4.10.3          zoo_1.8-12            
##  [22] cachem_1.0.8           igraph_1.6.0           mime_0.12             
##  [25] lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3       
##  [28] Matrix_1.6-4           R6_2.5.1               fastmap_1.1.1         
##  [31] fitdistrplus_1.1-11    future_1.33.1          shiny_1.8.0           
##  [34] digest_0.6.33          colorspace_2.1-0       tensor_1.5            
##  [37] RSpectra_0.16-1        irlba_2.3.5.1          labeling_0.4.3        
##  [40] progressr_0.14.0       fansi_1.0.6            spatstat.sparse_3.0-3 
##  [43] httr_1.4.7             polyclip_1.10-6        abind_1.4-5           
##  [46] compiler_4.3.1         bit64_4.0.5            withr_2.5.2           
##  [49] doParallel_1.0.17      fastDummies_1.7.3      highr_0.10            
##  [52] MASS_7.3-60            tools_4.3.1            lmtest_0.9-40         
##  [55] httpuv_1.6.13          future.apply_1.11.1    goftest_1.2-3         
##  [58] quadprog_1.5-8         glue_1.6.2             nlme_3.1-164          
##  [61] promises_1.2.1         grid_4.3.1             Rtsne_0.17            
##  [64] cluster_2.1.6          reshape2_1.4.4         generics_0.1.3        
##  [67] hdf5r_1.3.8            gtable_0.3.4           spatstat.data_3.0-3   
##  [70] tidyr_1.3.0            data.table_1.14.10     utf8_1.2.4            
##  [73] spatstat.geom_3.2-7    RcppAnnoy_0.0.21       ggrepel_0.9.4         
##  [76] RANN_2.6.1             foreach_1.5.2          pillar_1.9.0          
##  [79] stringr_1.5.1          spam_2.10-0            RcppHNSW_0.5.0        
##  [82] later_1.3.2            splines_4.3.1          dplyr_1.1.4           
##  [85] lattice_0.22-5         bit_4.0.5              survival_3.5-7        
##  [88] deldir_2.0-2           tidyselect_1.2.0       miniUI_0.1.1.1        
##  [91] pbapply_1.7-2          gridExtra_2.3          scattermore_1.2       
##  [94] RhpcBLASctl_0.23-42    xfun_0.41              matrixStats_1.2.0     
##  [97] stringi_1.8.3          lazyeval_0.2.2         yaml_2.3.8            
## [100] evaluate_0.23          codetools_0.2-19       tibble_3.2.1          
## [103] BiocManager_1.30.22    cli_3.6.2              uwot_0.1.16           
## [106] xtable_1.8-4           reticulate_1.34.0      munsell_0.5.0         
## [109] jquerylib_0.1.4        harmony_1.2.0          Rcpp_1.0.11           
## [112] globals_0.16.2         spatstat.random_3.2-2  png_0.1-8             
## [115] parallel_4.3.1         ellipsis_0.3.2         dotCall64_1.1-1       
## [118] listenv_0.9.0          scales_1.3.0           ggridges_0.5.5        
## [121] leiden_0.4.3.1         purrr_1.0.2            rlang_1.1.2           
## [124] cowplot_1.1.2
<div class='r_output'>