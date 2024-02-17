library(Seurat)
library(SeuratData)
library(dplyr)

pbmc.data <- Read10X(data.dir = "~/milarge_files/Klatzmann_inner/SC_motiffs/filtered_feature_bc_matrix/")
pbmc.data <- Read10X_h5(filename = "10k_PBMC_5pv2_nextgem_Chromium_X_Multiplex_count_raw_feature_bc_matrix.h5")
pbmc.data <- Read10X_h5(filename = "~/milarge_files/Klatzmann_inner/SC_motiffs/20k_PBMC_5pv2_HT_nextgem_Chromium_X_Multiplex_count_raw_feature_bc_matrix.h5")

pbmc.data <- pbmc.data$`Gene Expression`

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "Klazmann_comp", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
levels(pbmc)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:7)
pbmc <- FindNeighbors(pbmc, dims = 1:10)

#pbmc <- FindNeighbors(pbmc, dims = 1:8)
pbmc <- FindClusters(pbmc, resolution = 0.6)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
#pbmc <- RunUMAP(pbmc, dims = 1:8)

DimPlot(pbmc, reduction = "umap")

cell_names <- read.csv('~/milarge_files/Klatzmann_inner/SC_motiffs/clustered_barcodes_alpha_3.csv')
cell_names <- cell_names$X0
cur_names <- pbmc@assays$RNA$counts@Dimnames[[2]]
cell_names<-intersect(cell_names, cur_names)

DimPlot(pbmc, reduction = "umap", cells.highlight = cell_names)

test.markers <- FindMarkers(pbmc,
                            ident.1 = cell_names,
                            min.pct = 0.25, test.use = "MAST")

test.markers <- FindMarkers(pbmc,
                            ident.1 = 1,
                            min.pct = 0.25, test.use = "MAST")


head(test.markers, n = 20)

top_genes <- rownames(head(test.markers, n = 16))


FeaturePlot(pbmc, features = c("CD8A", "CD8B", "LINC02446", 'CD3E', 'RPS14', 'RPSA', 'LEF1', 'CD27', 'RPL14', 'CD4', 'PTPRC', 'NCAM1', 'CD16'))

FeaturePlot(pbmc, features = c("CD8A", "CD8B", "CD4", 'CD3E', 'CD3G', 'CD3D'))


FeaturePlot(pbmc, features = c('CD4', 'CD8A', "CD5", "CD69", "CD24", 'CD62L', 'RAG1', 'CXCR4', 'CCR7', 'GATA3', 'ZBTB7B', 'RUNX3', 'PRNP'))

FeaturePlot(pbmc, features = c('CD4', 'CD8A', "CD5", "CD69", "CD24", 'CD62L', 'RAG1', 'CXCR4', 'CCR7', 'GATA3', 'ZBTB7B', 'RUNX3', 'PRNP'))

FeaturePlot(pbmc, features = c('CACNA1E', 'IL7R', "MPP4", "NKG7", "RFLNB", 'ITGAE', 'IKZF2', 'RAPSN', 'SLIT3', 'PLAUR', 'S1PR1', 'GPR68', 'F13A1', 'SLPI'))


FeaturePlot(pbmc, features = c('CACNA1E', 'IL7R', "MPP4", "NKG7", "RFLNB", 'ITGAE', 'IKZF2', 'RAPSN', 'SLIT3', 'PLAUR', 'S1PR1', 'GPR68', 'F13A1', 'SLPI'))

FeaturePlot(pbmc, features = c('TRAV10'))


FeaturePlot(pbmc, features = top_genes)

write.csv(x = test.markers, file = 'mast_res_cd8_alpha_all_3.csv')
test.markers

write.csv(pbmc$seurat_clusters, file = 'SC_clusters.csv')
