suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(Seurat)
  library(patchwork)
})

### Utility Functions ###

## typeline
# take user input
# msg: message to display (prompt)
# output: text entered by user
typeline <- function(msg) {
  if (interactive() ) {
    txt <- readline(msg)
  } else {
    cat(msg);
    txt <- readLines("stdin",n=1);
  }
  return(txt)
}

## userEnterList
# get list of entries from the user
# msg: messsage to display (prompt)
# output: text entered by user
userEnterList <- function(msg) {
  if (interactive() ) {
    EXP <- as.integer(strsplit(readline(msg), " ")[[1]])
  } else {
    cat(msg)
    EXP <- as.integer(strsplit(readLines("stdin",n=1), " ")[[1]])
  }
}

### identifyNormalRNAEP
# custom function for Epithelioid Sarcoma to enable the user to identify the normal cell populations in a sample based on expression of SmarcB1, Rgs5, and Pdgfra
# cells are labeled as tumor or normal and this label is stored in "tumNormLabels" in the seurat object
# suob: seurat object containing epithelioid sarcoma expression levels
# output: seurat object with labels
identifyNormalRNAEP <- function(suob) {
  X11()
  p <- FeaturePlot(suob, features = c("Pdgfra","Rgs5","Smarcb1"))
  print(p)
  normclusters <- userEnterList("enter the non-tumor clusters (space-separated list) \n")
  suob$tumNormLabels <- ifelse(suob$seurat_clusters %in% normclusters, "normal", "tumor")
  print("Identified normal cell clusters: ", normclusters)
  return(suob)
}

### Major Functionality ###

## ReadRNA
# wrapper for reading 10X genomics library into seurat
# path.to.levels: path to output of 10X cellranger containing matrix.mtx, barcodes.tsv, features.tsv, and genes.tsv (usually located in filtered_gene_bc_matrices/genome/)
# output: seurat object
ReadRNA <- function(path.to.levels) {
  suob.rna <- Read10X(data.dir = path.to.levels)
  return(suob.rna)
}

### TODO: depreciate runSeuratRNA, move all preprocessing steps into filterScaleRNA, make new function that performs only the PCA and clustering (or one that only does PCA, and only does clustering, better options when replacing reduction)

## runSeuratRNA
# runs the standard seurat RNA filtering/clustering pipeline, generating relevant plots along the way
# suob: seurat object
# verbose: if TRUE, prints full set of figures with more details
# standard_filters: if FALSE, prompts user for filter values and clustering hyperparameters, otherwise uses standard values selected based on experience (if not sure, leave as false)
# usevariablefeatures: set to true causes the scaling of data using only variable features, false will use all genes
# outputs: seurat object which has been filtered, dimensionally reduced, and clustered
runSeuratRNA <- function(suob, verbose=FALSE, standard_filters=FALSE, usevariablefeatures = TRUE) {
  suob[["percent.mt"]] <- PercentageFeatureSet(suob, pattern = "^mt-")
  X11()
  p = VlnPlot(suob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  print(p)
  X11()
  plot1 <- FeatureScatter(suob, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(suob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p = plot1 + plot2
  print(p)
  if(standard_filters) {
    mtcutoff <- 10
    nfeatmin <- 3000
    nfeatmax <- 8000
    ncountmin <- 8000
  } else {
    mtcutoff <- typeline("Filtering: Enter value of mt-percent cutoff (default = 5): ")
    print(strtoi(mtcutoff, base = 0L))
    nfeatmin <- typeline("Filtering: Enter minimum features (default = 200): ")
    print(strtoi(nfeatmin, base = 0L))
    nfeatmax <- typeline("Filtering: Enter maximum features: ")
    print(strtoi(nfeatmax, base = 0L))
    ncountmin <- typeline("Filtering: Enter minimum count: ")
    print(strtoi(ncountmin, base = 0L))
  }
  old.cell.count <- ncol(suob[["RNA"]]@counts)
  suob <- subset(suob, subset = nFeature_RNA > strtoi(nfeatmin, base = 0L) & nFeature_RNA < strtoi(nfeatmax, base = 0L) & percent.mt < strtoi(mtcutoff, base = 0L)& nCount_RNA > strtoi(ncountmin, base = 0L))
  filtered.cell.count <- ncol(suob[["RNA"]]@counts)
  cat("Cells before filtering: ", old.cell.count)
  cat("Cells after filtering: ", filtered.cell.count)
  X11()
  p = VlnPlot(suob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  print(p)
  suob <- NormalizeData(suob)
  suob <- FindVariableFeatures(suob, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(suob), 10)
  if(verbose) {
    X11()
    par(mfrow=c(1,2))
    plot1 = VariableFeaturePlot(suob)
    p = LabelPoints(plot = plot1, points = top10, repel = TRUE)
    print(p)
  }
  if(usevariablefeatures) {
    genes <- VariableFeatures(suob)
  } else {
    genes <- rownames(suob)
  }
  suob <- ScaleData(suob, features = genes)
  #num_dims <- typeline("PCA: Enter number PCA dimensions (default = 50): ")
  num_dims <- 50
  suob <- RunPCA(suob, features = VariableFeatures(object = suob), npcs = strtoi(num_dims, base = 0L))
  print(suob[["pca"]], dims = 1:5, nfeatures = 5)
  if(verbose) {
    X11()
    p = VizDimLoadings(suob, dims = 1:2, reduction = "pca")
    print(p)
    X11()
    p = DimPlot(suob, reduction = "pca")
    print(p)
    X11()
    p = DimHeatmap(suob, dims = 1:15, cells = 500, balanced = TRUE, fast=FALSE)
    print(p)
  }
  X11()
  p = ElbowPlot(suob, ndims = 50)
  print(p)
  if(standard_filters) {
    clu_dims <- 50
  } else {
    clu_dims <- typeline("Clustering: Enter number of PCA dims (consider elbow plot): ")
    print(strtoi(clu_dims))
  }
  #res <- typeline("Clustering: Enter resolution (between 0.4-1.2, higher for larger datasets): ")
  #as.double(res)
  res <- 0.8
  suob <- FindNeighbors(suob, dims = 1:strtoi(clu_dims, base = 0L))
  suob <- FindClusters(suob, resolution = as.double(res))
  head(Idents(suob), 5)
  #map <- typeline("What map to use? (enter umap or tsne): ")
  map <- "umap"
  if (map == "umap") {
    suob <- RunUMAP(suob, dims = 1:strtoi(clu_dims))
    X11()
    p = DimPlot(suob, reduction="umap")
    print(p)
  } else if (map == "tsne") {
    suob <- RunTSNE(suob, dims.use = 1:strtoi(clu_dims))
    X11()
    p = DimPlot(suob, reduction="tsne")
    print(p)
  }
  return(suob)
}

## reclusterSeuratRNA
# runs seurat findclusters with louvain to regenerate clusters (and corresponding Idents of seurat object), meant to enable re-clustering at a finer or more granular resolution to control the number of clusters
# suob: seurat object
# res: resolution, required, default for runSeuratRNA is 0.8
# plot: generate a umap dimplot with the new clusters
# outputs: seurat object which has been reclustered
reclusterSeuratRNA <- function(suob, res, plot=TRUE) {
  suob <- FindClusters(object = suob, resolution = as.double(res))
  if (plot) {
    X11()
    p = DimPlot(suob, reduction="umap")
    print(p)
  }
  return(suob)
}

## subsetSeuratRNA
# subset a seurat object with a specific list of genes and then re-cluster explicitly for those genes
# suob: seurat object
# genes: list of genes to keep
# res: resolution, required, default for runSeuratRNA is 0.8
# plot: generate a umap dimplot with the new clusters
# verbose: generate elbowplot for new PCA
subsetSeuratRNA <- function(suob, genes, res, plot=TRUE, verbose=FALSE) {
  suob.subset <- subset(suob, features = genes)
  suob.subset <- NormalizeData(suob.subset)
  all.genes <- rownames(suob.subset)
  suob.subset <- ScaleData(suob.subset, features = all.genes)
  suob.subset <- RunPCA(suob.subset, features = all.genes)
  if (verbose) {
    X11()
    p = ElbowPlot(suob.subset)
    print(p)
  }
  suob.subset <- FindNeighbors(suob.subset, dims = 1:length(suob.subset@reductions$pca))
  suob.subset <- FindClusters(suob.subset, resolution = res)
  suob.subset <- RunUMAP(suob.subset, dims = 1:length(suob.subset@reductions$pca))
  if (plot) {
    X11()
    p = DimPlot(suob.subset, reduction="umap")
    print(p)
  }
  return(suob.subset)
}


### TODO: add description
filterScaleRNA <- function(suob, verbose=FALSE, standard_filters=FALSE, usevariablefeatures = TRUE, nvarfeatures = 2000) {
  suob[["percent.mt"]] <- PercentageFeatureSet(suob, pattern = "^mt-")
  X11()
  p = VlnPlot(suob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  print(p)
  X11()
  plot1 <- FeatureScatter(suob, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(suob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p = plot1 + plot2
  print(p)
  if(standard_filters) {
    mtcutoff <- 10
    nfeatmin <- 3000
    nfeatmax <- 8000
    ncountmin <- 8000
  } else {
    mtcutoff <- typeline("Filtering: Enter value of mt-percent cutoff (default = 5): ")
    print(strtoi(mtcutoff, base = 0L))
    nfeatmin <- typeline("Filtering: Enter minimum features (default = 200): ")
    print(strtoi(nfeatmin, base = 0L))
    nfeatmax <- typeline("Filtering: Enter maximum features: ")
    print(strtoi(nfeatmax, base = 0L))
    ncountmin <- typeline("Filtering: Enter minimum count: ")
    print(strtoi(ncountmin, base = 0L))
  }
  old.cell.count <- ncol(suob[["RNA"]]@counts)
  suob <- subset(suob, subset = nFeature_RNA > strtoi(nfeatmin, base = 0L) & nFeature_RNA < strtoi(nfeatmax, base = 0L) & percent.mt < strtoi(mtcutoff, base = 0L)& nCount_RNA > strtoi(ncountmin, base = 0L))
  filtered.cell.count <- ncol(suob[["RNA"]]@counts)
  cat("Cells before filtering: ", old.cell.count)
  cat("Cells after filtering: ", filtered.cell.count)
  X11()
  p = VlnPlot(suob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  print(p)
  suob <- NormalizeData(suob)
  ### How much difference does nfeatures make on the outcome of clustering (?!)
  suob <- FindVariableFeatures(suob, selection.method = "vst", nfeatures = nvarfeatures)
  top10 <- head(VariableFeatures(suob), 10)
  if(verbose) {
    X11()
    par(mfrow=c(1,2))
    plot1 = VariableFeaturePlot(suob)
    p = LabelPoints(plot = plot1, points = top10, repel = TRUE)
    print(p)
  }
  if(usevariablefeatures) {
    genes <- VariableFeatures(suob)
  } else {
    genes <- rownames(suob)
  }
  suob <- ScaleData(suob, features = genes)
  return(suob)
}