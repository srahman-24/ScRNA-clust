library(dplyr)
library(Seurat)
library(patchwork)
library(CellBench)
library(scran)
library(ggplot2)
library(scater)
library(flexclust)
library(LICORS)
library(SingleCellExperiment)
library(RJcluster)


###### Some historic perspective on Single-cell experiment data 
## The `SingleCellExperiment` class is a lightweight Bioconductor container for storing and manipulating single-cell 
## genomics data. Rows ahould represent genes or genomic regions. And column should represent cells. 
## It provides methods for storing dimensionality reduction results and data for alternative feature sets ( synhetc spike-in transcripts) 
## It is the central data structure for Bioconductor single-cell packages like `scater` and `scran`. 

#### preparing single-cell  sce OBJECT ######

counts = matrix(rpois(100, lambda = 10), ncol = 10, nrow = 10)
sce    = SingleCellExperiment::SingleCellExperiment(list(counts = counts))
sce
pretend.cell.labels  <- sample(letters, ncol(counts), replace = TRUE)
pretend.gene.lengths <- sample(10000, nrow(counts))
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
                            colData = DataFrame(label = pretend.cell.labels),
                            rowData = DataFrame(length = pretend.gene.lengths),
                            metadata = list(study = "GSE111111")
)
sce

libsizes       <- colSums(counts)
size.factors   <- libsizes/mean(libsizes)
logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)
sce

pca_data <- prcomp(t(logcounts(sce)), rank = 10)
library(Rtsne)
set.seed(5252)
tsne_data <- Rtsne::Rtsne(pca_data$x[,1:10], pca = FALSE, perplexity = 1)

reducedDims(sce) <- list(PCA = pca_data$x, TSNE = tsne_data$Y)
sce


##### assays : users can assign arbitrary names to entries of assays. standard names are used for assays to 
##### have parity between packages. 
##### counts     =  Raw counts data representing the number of reads or transcripts for a particular genes. 
##### normcounts =  Normalized values on the same scale as the original counts : counts/cell size 
##### logcounts  = Log-transformed counts : log2 transform on norm-counts. 
##### cpm        = counts-per-million. 


guo.data <- read.csv("/Users/srahman/Downloads/Guo_O-A-1.csv", header = TRUE, row.names = 1)
#pbmc     <- CreateSeuratObject(pbmc.data, min.cells = 3, min.features = 200)
dat     <- data.matrix(guo.data)
#colnames(guo.data)
#rownames(guo.data)
#imnames(dat) <- NULL 
#class(dat)
#dim(dat)
#sce    = SingleCellExperiment::SingleCellExperiment(list(counts = dat))
sce    = SingleCellExperiment::SingleCellExperiment(list(counts = dat),
                     colData = data.frame(label = colnames(guo.data)),
                     rowData = data.frame(length = rownames(guo.data)),
                     metadata = list(study = "Guo-O-A-1")
)
sce
libsizes       <- colSums(dat)
size.factors   <- libsizes/mean(libsizes)
SingleCellExperiment::logcounts(sce) <- log2(t(t(dat)/size.factors) + 1)

save(sce, file = "sce_O_A_1.RData")

#pca_data <- prcomp(t(logcounts(sce)), rank = 50)
#library(Rtsne)
#set.seed(5252)
#tsne_data <- Rtsne::Rtsne(pca_data$x[,1:50], pca = FALSE, perplexity = 10)

#reducedDims(sce) <- list(PCA = pca_data$x, TSNE = tsne_data$Y)
#sce

 
scran_high_var = function(sce,topn = 1000){
  var.fit <- trendVar(sce, method = "loess", use.spikes = FALSE)
  var.out <- decomposeVar(sce, var.fit)
  hvg.out <- var.out[order(var.out$bio, decreasing = TRUE)[1:topn], ]
  rowData(sce)$hi_var = FALSE
  rowData(sce)$hi_var[rownames(rowData(sce)) %in% rownames(hvg.out)] = TRUE
  return(sce)
}

sce1      = scran_high_var(sce)
var.genes = rownames(sce1)[rowData(sce1)$hi_var]
tmp       = logcounts(sce1)[var.genes,]
data      = t(as.matrix(tmp))
sdata      = scale(data, center = apply(data, 2, mean), scale =  T)
boxplot(sdata[1:100,1:100])
round(colMeans(sdata),3)
apply(sdata, 2, sd)

library(Rcpp)

library(RJcluster)
#data  = matrix(rnorm(100000, 0, 1), nrow = 1000, ncol = 100)
RJres = RJcluster::RJclust(sdata, penalty = "hockey_stick", scaleRJ = TRUE, C_max = 20, n_bins = 200)
RJres$K
save(RJres, file = 'RJres_0_1.RData')

RJres_bic = RJcluster::RJclust(sdata, penalty = "bic", scaleRJ = TRUE, C_max = 20, n_bins = 200)
RJres_bic$K

library(Rtsne)
library(M3C)
umap(tmp, labels = as.factor(RJres$class), controlscale = TRUE  , scale = 3, legendtitle = 'clusters')


##### Conversion from SingleCellExperiment to Seurat objects 

sce_seurat = as.Seurat(sce, counts = "counts", data = "logcounts")

data      <- NormalizeData(sce_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
data      <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(data)
data      <- ScaleData(data, features = all.genes)
data      <- RunPCA(data, features = VariableFeatures(object = data))

data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.03)
data <- RunUMAP(data, dims = 1:10)
data <- RunTSNE(data, seed.use = 5567, tsne.method = "Rtsne", dims = 1:10)

##### Visualization 
DimPlot(data, reduction = "umap")
DimPlot(data, reduction = "tsne")

table(data$seurat_clusters, RJres$class)

data <- AddMetaData(object = data, metadata = RJres$class, col.name = "RJcluster")
DimPlot(data, reduction = "tsne", group.by = "RJcluster", label = T)
save(data, file = "Guo_0_1_seurat003.RData")

##### Identify the 10 most highly variable genes 

top10 <- head(VariableFeatures(data), 10)
top10

plot(RJres$penalty, ylab = "Hockey Stick", pch = 19, main = "Guo_0_1")
abline(v = 7, lwd = 3, col = "red")



### Find Biomarkers for each clusters 
cluster1.markers <- FindMarkers(data, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_logFC)

#### Visualization 

FeaturePlot(data, features = c("ENSMUSG00000035042", "ENSMUSG00000004612", "ENSMUSG00000003379", 
                               "ENSMUSG00000020154", "ENSMUSG00000029304", "ENSMUSG00000069516", 
                               "ENSMUSG00000056071", "ENSMUSG00000057465", "ENSMUSG00000064373")) 


par(mfrow = c(4,4))
VlnPlot(data, features = c("ENSMUSG00000035042"), log = TRUE, pt.size = 0.001) #0
VlnPlot(data, features = c("ENSMUSG00000004612"), log = TRUE, pt.size = 0.001) #0
VlnPlot(data, features = c("ENSMUSG00000003379"), log = TRUE, pt.size = 0.001) #1
VlnPlot(data, features = c("ENSMUSG00000024610"), log = TRUE, pt.size = 0.001) #1
VlnPlot(data, features = c("ENSMUSG00000020154"), log = TRUE, pt.size = 0.001) #2
VlnPlot(data, features = c("ENSMUSG00000064373"), log = TRUE, pt.size = 0.001) #2
VlnPlot(data, features = c("ENSMUSG00000029304"), log = TRUE, pt.size = 0.001) #3
VlnPlot(data, features = c("ENSMUSG00000022037"), log = TRUE, pt.size = 0.001) #3
VlnPlot(data, features = c("ENSMUSG00000069516"), log = TRUE, pt.size = 0.001) #4
VlnPlot(data, features = c("ENSMUSG00000040809"), log = TRUE, pt.size = 0.001) #4
VlnPlot(data, features = c("ENSMUSG00000056071"), log = TRUE, pt.size = 0.001) #5
VlnPlot(data, features = c("ENSMUSG00000056054"), log = TRUE, pt.size = 0.001) #5
VlnPlot(data, features = c("ENSMUSG00000057465"), log = TRUE, pt.size = 0.001) #6
VlnPlot(data, features = c("ENSMUSG00000078672"), log = TRUE, pt.size = 0.001) #6
