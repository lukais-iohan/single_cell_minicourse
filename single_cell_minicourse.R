library('Seurat')
library('dplyr')
library('slingshot')
library('DelayedMatrixStats')
library('RColorBrewer')

### Importar matriz de expressao

pbmc_matriz <- readRDS('dados/pbmc_matriz.rds')


### Como é representada a matriz de expressao?

pbmc_matriz[c("CD3D", "TCL1A", "MS4A1"), 1:10]

### Criar objeto Seurat

pbmc_obj <- CreateSeuratObject(counts = pbmc_matriz,project = 'pbmc3k',
            min.cells = 3, min.features = 200)

pbmc_obj

### QC e filtrar células

## Mesmo após o QC depois do sequenciamento, ainda podemos ter células de baixa qualidade
## i.e, empty droplets, doublets, ou ainda células em estado de apoptose (alta quantidade MT-genes)

## Calcular porcentagem de genes mitocondriais

pbmc_obj[["percent.mt"]] <- PercentageFeatureSet(pbmc_obj,pattern = 'MT-')

## Onde as métricas são armazenadas: metadata

### Visualizar métricas

# Visualize QC metrics as a violin plot
VlnPlot(pbmc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

plot1 <- FeatureScatter(pbmc_obj, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")

plot2 <- FeatureScatter(pbmc_obj, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2

## Subset

pbmc_obj <- subset(pbmc_obj, 
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc_obj

## Normalizacao Normalizacao - Scale Data - FindVariable Features

pbmc_obj <- SCTransform(pbmc_obj, vars.to.regress = 'percent.mt')

## Redução de dimensionalidade - PCA

pbmc_obj <- RunPCA(pbmc_obj,features = VariableFeatures(object = pbmc_obj))

DimPlot(pbmc_obj, reduction = "pca") + NoLegend()
ElbowPlot(pbmc_obj,ndims = 30)

### Clusterização das células

pbmc_obj <- FindNeighbors(pbmc_obj, dims = 1:10)
pbmc_obj <- FindClusters(pbmc_obj, resolution = c(0.4,0.6,0.8,1,1.2))

### Run Umap

pbmc_obj <- RunUMAP(pbmc_obj, dims = 1:10)

Idents(pbmc_obj) <- 'SCT_snn_res.0.4'

DimPlot(pbmc_obj, reduction = "umap",
        label = TRUE, pt.size = 0.5) + NoLegend()

### Expressao Diferencial

pbmc.markers <- FindAllMarkers(pbmc_obj, only.pos = TRUE)
pbmc.markers <- pbmc.markers %>% group_by(cluster) %>%
  filter(avg_log2FC > 1,p_val_adj < 0.01) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

top5 <- top5[!duplicated(top5$gene),]

DoHeatmap(pbmc_obj,features = top5$gene)


# 4 == CD8+ T
# 2 == CD14+ Mono
# 1 == Naive CD4+ T
# 0 == Memory CD4+
# 3 == B
# 4 == CD8+ T
# 6 == FCGR3A+ Mono
# 7 == DC
# 5 == NK
# 8 = Platelet

## O - Naive CD4+T
## 

p1 <- FeaturePlot(pbmc_obj, features = c('IL7R','CCR7'))
p2 <- DimPlot(pbmc_obj, reduction = "umap", 
              label = TRUE, pt.size = 0.5) + NoLegend()

p1+p2

new_cluster_ids <- c('Memory CD4+','Naive CD4+ T','CD14+ Mono',
                     'B','CD8+ T', 'NK','FCGR3A+ Mono','DC','Platelet')
names(new_cluster_ids) <- levels(pbmc_obj)
pbmc_obj <- RenameIdents(pbmc_obj, new_cluster_ids)
DimPlot(pbmc_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


pbmc_obj_sce <- as.SingleCellExperiment(pbmc_obj)

pbmc_obj_sce_sling <- slingshot(reducedDim(pbmc_obj_sce, "UMAP")[,1:2], 
                            clusterLabels = pbmc_obj$SCT_snn_res.0.4, 
                            start.clus = 0)

plot(reducedDims(pbmc_obj_sce)$UMAP, 
     col = brewer.pal(9,'Set1')[pbmc_obj_sce$SCT_snn_res.0.4], 
     pch=16, asp = 1)
lines(SlingshotDataSet(pbmc_obj_sce_sling), 
      lwd=2, type = 'lineages', col = 'black')

#### 

mice <- readRDS('dados/mice_neuro.rds')
