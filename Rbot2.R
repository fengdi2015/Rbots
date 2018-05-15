library(akima);library(MASS);library(tsne)
colors = rainbow(length(unique(iris$Species)))
names(colors) = unique(iris$Species)
ecb = function(x,y) {plot(x,t='n');text(x,labels=iris$Species, col=colors[iris$Species]) }
tsne_iris = tsne(iris[,1:4], epoch_callback = ecb, perplexity=50)
data=tsne_iris;
###
library(Seurat)
library(Matrix)
pbmc33k.data <- Read10X(data.dir = "10data/")
pbmc33k  <- CreateSeuratObject(raw.data = pbmc33k.data, min.cells = 3, project = "10X_PBMC33K", names.field = 2, names.delim = "\\-")
# Calculate percentage of mitochondrial genes for filtering and add into object@meta.data,
# a great place to stash QC stats
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc33k@data), value = TRUE)
percent.mito <- colSums(pbmc33k@data[mito.genes, ]) / colSums(pbmc33k@data)
pbmc33k <- AddMetaData(object = pbmc33k, metadata = percent.mito, col.name = "percent.mito")
Seurat::VlnPlot(object = pbmc33k, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# We filter out cells that have unique gene counts over 2,500 and under 500, and > pbmc33k5% mitochondrial percentage
pbmc33k <- FilterCells(object = pbmc33k, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(2500, 0.05))

# Normalize the data
pbmc33k<- NormalizeData(object = pbmc33k)

# Choose gene outliers on mean-variability plot
pbmc33k <- FindVariableGenes(object = pbmc33k, x.low.cutoff = 0, y.cutoff = 0.8)
length(x = pbmc33k@var.genes)

# Perform negative-binomial regression on the variable genes, this sets their value in pbmc33k@scale.data, which is used for PCA/clustering
# We only do this on the variable genes to save time, but you can do this genome-wide
# We treat mitochondrial percentage, batch, and nUMI as confounding variables,

pbmc33k <- ScaleData(object = pbmc33k, vars.to.regress = c("percent.mito", "orig.ident", "nUMI"), genes.use = pbmc33k@var.genes, model.use = "negbinom")

# you can save the object at any time to save results, and can restore it back in using load()
#save(pbmc33k, file = "~/Projects/datasets/pbmc33k/pbmc33k_tutorial.Robj")

# Run PCA with the IRLBA package (iteratively computes the top dimensions, dramatic increase in speed since we only use a fraction of the PCs anyways)
# if you see the warning "did not converge–results might be invalid!; try increasing maxit or fastpath=FALSE", try increasing maxit
pbmc33k <- RunPCA(object = pbmc33k, pc.genes = pbmc33k@var.genes, pcs.compute = 40, pcs.print = 1:30, maxit = 500, weight.by.var = FALSE)

PCElbowPlot(object = pbmc33k, num.pc = 40)
PrintPCA(object = pbmc33k, pcs.print = 1:36)
PCHeatmap(object = pbmc33k, pc.use = 1:4,100)
PCHeatmap(object = pbmc33k, pc.use = 5:8,100)
PCHeatmap(object = pbmc33k, pc.use = 9:12,100)
# select 25 PCs for downstream analysis
pbmc33k <- RunTSNE(object = pbmc33k, dims.use = 1:25, do.fast = TRUE)
pbmc33k <- FindClusters(object = pbmc33k, reduction.type = "pca", dims.use = 1:25, resolution = 4, save.SNN = TRUE)
# save.SNN means that you can easily re-run with different resolution values.
# Here we run with a few different res values
### color test
require( gplots )
barplot( rep(1,100), width = rep(2,100) , col=rich.colors(100), space = 0 , border=0, axes=FALSE)

require( RColorBrewer )
display.brewer.pal(11 , "Spectral" )

require( colorRamps )
image(matrix(1:400, 20), col = blue2green2red(400) , axes = FALSE)

mycolor=blue2green2red(400)
##################################################
data=pbmc33k@dr$tsne@cell.embeddings
result=kde2d(data[,1],data[,2])
pdf("test.pdf")
image(result,col=blue2green2red(400),add=FALSE,xlab="",ylab="")
dev.off()
## create a fake height
n=length(data[,1])
h=sample(1:20,n,replace=T)
hcutoff=15
h[h<hcutoff] <- 0
s=interp(data[,1],data[,2],h,nx=150,ny=150,linear=TRUE,extrap=TRUE,jitter.random = TRUE)
pdf("test2.pdf")
with(s,image(x,y,z,col=mycolor,xlab="",ylab=""),add=TRUE)
with(s,contour(result,add=T,col = "#a5a298",lwd=0.9,drawlabels = F))
dev.off()


