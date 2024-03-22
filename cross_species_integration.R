library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)

library(reticulate)
path_to_python <- ".../anaconda3/envs/scvi-env/bin/python"
use_python(path_to_python)

macro_mouse <- readRDS(file = "macro_mouse.rds")
DefaultAssay(macro_mouse)<-"RNA"
macro_mouse@meta.data$origin <- "mouse"

macro_axl <- readRDS(file = "macro_axl.rds")
DefaultAssay(macro_axl)<-"RNA"
macro_axl@meta.data$origin <- "axl"


macro_axl[['SCT']] <- NULL
macro_mouse[['SCT']] <- NULL
macro_all <- merge(macro_mouse, y = macro_axl, add.cell.ids = c("mouse", "axl"))

table(macro_all@meta.data$origin)

macro_all <- NormalizeData(macro_all)
macro_all <- FindVariableFeatures(macro_all)
macro_all <- ScaleData(macro_all)
macro_all <- RunPCA(macro_all)

############# without integration

macro_all <- FindNeighbors(macro_all, dims = 1:30, reduction = "pca")
macro_all <- FindClusters(macro_all, resolution = 2,
                          cluster.name = "unintegrated_clusters")
macro_all <- RunUMAP(macro_all, dims = 1:30, reduction = "pca",
                     reduction.name = "umap.unintegrated")

# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(macro_all, reduction = "umap.unintegrated", group.by = c("origin", "ident"))


# ############## CCAIntegration
# 
# macro_all <- IntegrateLayers(
#   object = macro_all, method = CCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = FALSE
# )
# 
# macro_all <- FindNeighbors(macro_all, reduction = "integrated.cca", dims = 1:30)
# macro_all <- FindClusters(macro_all, resolution = 2, cluster.name = "cca_clusters")
# macro_all <- RunUMAP(macro_all, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
# 
# p1 <- DimPlot(
#   macro_all,
#   reduction = "umap.cca",
#   group.by = c("origin",  "cca_clusters"),
#   combine = FALSE, label.size = 2
# )
# 
# 
# 
# ############## HarmonyIntegration
# 
# macro_all <- IntegrateLayers(
#   object = macro_all, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "harmony",
#   verbose = FALSE
# )
# 
# macro_all <- FindNeighbors(macro_all, reduction = "harmony", dims = 1:30)
# macro_all <- FindClusters(macro_all, resolution = 2, cluster.name = "harmony_clusters")
# 
# macro_all <- RunUMAP(macro_all, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
# p2 <- DimPlot(
#   macro_all,
#   reduction = "umap.harmony",
#   group.by = c("origin",  "harmony_clusters"),
#   combine = FALSE, label.size = 2
# )
# 
# wrap_plots(c(p1, p2), ncol = 2, byrow = F)
# 
# ############## RPCAIntegration
# 
# macro_all <- IntegrateLayers(
#   object = macro_all, method = RPCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.rpca",
#   k.weight = 70,
#   verbose = FALSE
# )
# 
# macro_all <- FindNeighbors(macro_all, reduction = "integrated.rpca", dims = 1:30)
# macro_all <- FindClusters(macro_all, resolution = 2, cluster.name = "rcpa_clusters")
# macro_all <- RunUMAP(macro_all, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rcpa")
# 
# p3 <- DimPlot(
#   macro_all,
#   reduction = "umap.rcpa",
#   group.by = c("origin",  "rcpa_clusters"),
#   combine = FALSE, label.size = 2
# )

############## FastMNNIntegration

macro_all <- IntegrateLayers(
  object = macro_all, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

macro_all <- FindNeighbors(macro_all, reduction = "integrated.mnn", dims = 1:30)
macro_all <- FindClusters(macro_all, resolution = 0.5,
                          cluster.name = "mnn_clusters")
macro_all <- RunUMAP(macro_all, reduction = "integrated.mnn",
                     dims = 1:30, reduction.name = "umap.mnn")

p4 <- DimPlot(
  macro_all,
  reduction = "umap.mnn",
  group.by = c("origin",  "mnn_clusters","Wound.day"),
  combine = TRUE, label.size = 2
)

p4

#wrap_plots(c(p3, p4), ncol = 2, byrow = F)

#wrap_plots(c(p1, p4), ncol = 2, byrow = F)

#Once integrative analysis is complete, you can rejoin the layers - which 
#collapses the individual datasets together and recreates the original counts 
#and data layers. You will need to do this before performing any 
#differential expression analysis. However, you can always resplit the layers
#in case you would like to reperform integrative analysis.

macro_all
macro_all <- JoinLayers(macro_all)
macro_all <- ScaleData(macro_all, features = rownames(macro_all))

macro_all

DimPlot(
  macro_all[,macro_all$origin == "mouse"],
  reduction = "umap.mnn",
  group.by = c("origin",  "mnn_clusters"),
  combine = FALSE, label.size = 2
)

DimPlot(
  macro_all,
  reduction = "umap.mnn",
  group.by = c(  "mnn_clusters"),
  split.by = "origin",
  label.size = 2
)


macro_all@meta.data %>%
  group_by(mnn_clusters,origin) %>%
  count() %>%
  group_by(mnn_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=mnn_clusters,y=percent, fill=origin)) +
  geom_col() +
  ggtitle("Percentage of species per cluster")


FindMarkers(macro_all ,ident.1 = 0, min.pct = 0.25)

VlnPlot(macro_all ,features=c("SELENOP","C1QB","RPS29","SLC34A2"  ) )

# This loop just runs the FindMarkers function on all of the clusters
lapply(
  levels(macro_all[["mnn_clusters"]][[1]]),
  function(x)FindMarkers(macro_all,ident.1 = x,min.pct = 0.25)
) -> cluster.markers


# This simply adds the cluster number to the results of FindMarkers
sapply(0:(length(cluster.markers)-1),function(x) {
  cluster.markers[[x+1]]$gene <<- rownames(cluster.markers[[x+1]])
  cluster.markers[[x+1]]$cluster <<- x
})


mouse <- subset(x = macro_all , subset = origin == "mouse")
mouse <- ScaleData(mouse, features = rownames(mouse))

axl <- subset(x = macro_all , subset = origin == "axl")
axl <- ScaleData(axl, features = rownames(axl))

mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE)
mouse.markers %>%
  group_by(cluster) 

mouse.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(mouse, features = top10$gene) + NoLegend()


axl.markers <- FindAllMarkers(axl, only.pos = TRUE)
axl.markers %>%
  group_by(cluster) 

axl.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(axl, features = top10$gene) + NoLegend()


##### Markers for 4 different paths
library(readxl)

# read_excel reads both xls and xlsx files
path_markers <- read_excel("mac_path_markers.xlsx")

phagocytic_markers_late <- unname(unlist(path_markers[path_markers$`Target stage`== "Late.P1" ,c( "Gene"   )]))
phagocytic_markers_late <- c("CSF1R", "CD74" , "CD83"  )


##### Phagocytic
RidgePlot(macro_all, features = phagocytic_markers_late , ncol = 2)
VlnPlot(macro_all, features = phagocytic_markers_late, split.by = "origin" )
FeaturePlot(macro_all, features = phagocytic_markers_late, reduction = "umap.mnn" )
DotPlot(macro_all, features = phagocytic_markers_late, split.by = "origin") + RotatedAxis()
DoHeatmap(subset(macro_all, downsample = 100) , features = phagocytic_markers_late, size = 3, group.by = c("ident", "origin"))
DoHeatmap( macro_all , features = phagocytic_markers_late , size = 3)+ labs(title = "")

FeatureScatter(macro_all, feature1 = "CD63", feature2 ="FCGR3" )


macro_mouse <- NormalizeData(macro_mouse)
macro_mouse <- FindVariableFeatures(macro_mouse)
macro_mouse <- ScaleData(macro_mouse)
macro_mouse <- RunPCA(macro_mouse)

############# without integration

macro_mouse <- FindNeighbors(macro_mouse, dims = 1:30, reduction = "pca")
macro_mouse <- FindClusters(macro_mouse, resolution = 2, cluster.name = "unintegrated_clusters")
macro_mouse <- RunUMAP(macro_mouse, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(macro_mouse, reduction = "umap.unintegrated", group.by = c("origin", "ident"))

FeaturePlot(macro_mouse, features = inflammatory_markers   )
