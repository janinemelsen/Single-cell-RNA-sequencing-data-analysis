## packages
library(Seurat)
library(scran)
library(scater)
library(harmony)
library(batchelor)
library(BiocParallel)
library(BiocNeighbors)
library(pheatmap)
library(reshape2)
library(ggthemes)


# For this analysis 4 scRNAseq datasets were analyzed that included circulating NK cells
# GSE119562 1 sample (GSM3776778)
# GSE130430 2 samples (GSM3738542 and GSM3738543)
# GSE184329 5 samples (GSM5584151, GSM5584155, GSM5584156_1, GSM5584156_2, GSM5584156_3)
# GSE197037 7 samples (GSM5907300, GSM5907301, GSM5907302, GSM5907303, GSM5907304, GSM5907305, GSM5907306)


# Quality control was performed by clustering every individual sample and remove clusters with high mitochondrial content, low ribosomal content, low gene count and high gene count (doublets)
# SingleR was used to predict which cells are NK cells (using BlueprintEncodeData). Non-NK cells were removed from further analysis.

# This script starts at importing the filtered seurat objects (every object is one individual sample) 


### IMPORT DATASET  ###
file_paths <- list.files(path='path/filtered objects', pattern = ".Rds", full.names=TRUE)
filenames <-  gsub(pattern = ".Rds", replacement = "", x = basename(file_paths))
# Read all your data into a list
seurat_list <- lapply(file_paths, readRDS)
# Assign file names to list elements
names(seurat_list) <- filenames 
# convert to sce objects
sce_list <- lapply(seurat_list, as.SingleCellExperiment)
# lib size factors
lib.factor <- lapply(sce_list, function(x) librarySizeFactors(x))
#compute log-normalized expression values for downstream use (tis step can be skipped if you want to correct for sequencing depth).
for (i in 1:length(sce_list)){
  sce_list[[i]] <- logNormCounts(sce_list[[i]], size_factors=lib.factor[[i]])
}

#remove GSM76 and GSM77 (activated cells)
sce_list[["GSM3377676_filtered_by_cluster"]] <- NULL
sce_list[["GSM3377677_filtered_by_cluster"]] <- NULL

### CORRECTION SEQUENCING DEPTH  ###
# what are the overlapping genes?
minrow <- Reduce(intersect, lapply(sce_list, rownames))
sce_list_universal <- lapply(sce_list, function(x) x[minrow])
# correction sequencing depth
sce_list_rescaled <- multiBatchNorm(sce_list_universal, BPPARAM=SnowParam(workers=18))

### HIGHLY VARIABLE GENES ###
dec<- lapply(sce_list_rescaled, function(x) modelGeneVar(x, BPPARAM=SnowParam(workers=19)))
dec_combined <- do.call(combineVar, dec)
hvg.sce.var <- getTopHVGs(dec_combined, n=5000)
set.seed(1000101001)

### COMBINE DATASET IN 1 SCE ###
for (i in 1:length(sce_list_rescaled)){
  rowData(sce_list_rescaled[[i]]) <- NULL 
  altExps(sce_list_rescaled[[i]]) <- NULL
  colData(sce_list_rescaled[[i]]) <- colData(sce_list_rescaled[[i]])[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo", 'seurat_clusters')]
}
# correct orig.ident
sce_list_rescaled[[3]]$orig.ident <- 'CMVpos1_5_libA'
sce_list_rescaled[[4]]$orig.ident <- 'CMVpos1_5_libB'
# merge datasets into one sce object
merged <- do.call(cbind, sce_list_rescaled)

### INTEGRATION ###
# PCA (set to scale=TRUE for removal of clusters driven by highly expressed genes)
set.seed(100)
reducedDim(merged, 'PCA') <- NULL
reducedDim(merged, 'UMAP') <- NULL
merged <- runPCA(merged, subset_row=hvg.sce.var, name="PCA", exprs_values='logcounts', scale=TRUE)

  
# add Dataset
merged$Dataset <- merged$orig.ident
merged$Dataset <- gsub('GSM3377678$', 'Dataset1', merged$Dataset)
merged$Dataset <- gsub('GSM3738542$', 'Dataset2', merged$Dataset)
merged$Dataset <- gsub('GSM3738543$', 'Dataset2', merged$Dataset)
merged$Dataset <- gsub('GSM5584154$', 'Dataset3', merged$Dataset)
merged$Dataset <- gsub('GSM5584155$', 'Dataset3', merged$Dataset)
merged$Dataset <- gsub('GSM5584156_1$', 'Dataset3', merged$Dataset)
merged$Dataset <- gsub('GSM5584156_2$', 'Dataset3', merged$Dataset)
merged$Dataset <- gsub('GSM5584156_3$', 'Dataset3', merged$Dataset)
merged$Dataset <- gsub('CMVneg1_donorA$', 'Dataset4', merged$Dataset)
merged$Dataset <- gsub('CMVneg2_donorB$', 'Dataset4', merged$Dataset)
merged$Dataset <- gsub('CMVpos1_5_libA$', 'Dataset4', merged$Dataset)
merged$Dataset <- gsub('CMVpos1_5_libB$', 'Dataset4', merged$Dataset)
merged$Dataset <- gsub('CMVpos2_donorC$', 'Dataset4', merged$Dataset)
merged$Dataset <- gsub('CMVpos3_donorD$', 'Dataset4', merged$Dataset)
merged$Dataset <- gsub('CMVpos4_donorE$', 'Dataset4', merged$Dataset)

# add Chemistry
merged$Chemistry <- merged$orig.ident
merged$Chemistry <- gsub('GSM3377678$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('GSM3738542$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('GSM3738543$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('GSM5584154$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('GSM5584155$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('GSM5584156_1$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('GSM5584156_2$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('GSM5584156_3$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('CMVneg1_donorA$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('CMVneg2_donorB$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('CMVpos1_5_libA$', 'V3', merged$Chemistry)
merged$Chemistry <- gsub('CMVpos1_5_libB$', 'V3', merged$Chemistry)
merged$Chemistry <- gsub('CMVpos2_donorC$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('CMVpos3_donorD$', 'V2', merged$Chemistry)
merged$Chemistry <- gsub('CMVpos4_donorE$', 'V2', merged$Chemistry)

# Harmony integration Batch
merged <- RunHarmony(merged, group.by.vars=c('Dataset', 'Chemistry'), 
                     theta = NULL, lambda = NULL, sigma = 0.1,
                     nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                     max.iter.cluster = 20, epsilon.cluster = 1e-05,
                     epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                     reference_values = NULL, reduction.save = "HARMONY")



### CLUSTERING ###
set.seed(100)
g <- buildSNNGraph(merged, k=10, use.dimred = 'HARMONY' ,type='jaccard')
clustering <- igraph::cluster_louvain(g)$membership
merged$cluster <- factor(clustering)


### UMAP ###
merged <- runUMAP(merged, dimred="HARMONY",
                  name='UMAP',
                  ncomponents=2,
                  external_neighbors=TRUE, 
                  n_neighbors=30,
                  BNPARAM=AnnoyParam(),
                  BPPARAM=SnowParam(workers=18),
                  n_threads=bpnworkers(SnowParam(workers=18)),
                  min_dist=0.4)

### tSNE
merged <- runTSNE(merged, dimred="HARMONY",
                  name='TSNE',
                  ncomponents=2,
                  external_neighbors=TRUE, 
                  BNPARAM=AnnoyParam(),
                  BPPARAM=SnowParam(workers=18),
                  n_threads=bpnworkers(SnowParam(workers=18)))


### VISUALIZATION ###
plotColData(merged, x='nFeature_RNA', y='nCount_RNA', colour_by='orig.ident')
plotColData(merged, x='nFeature_RNA', y='nCount_RNA', colour_by='cluster')+geom_density_2d(colour = "black")
plotReducedDim(merged, dimred="UMAP", colour_by ="cluster", text_by = 'cluster')+ theme(text = element_text(size = 15), axis.text=element_text(size=15))+geom_density_2d(colour = "black")
plotReducedDim(merged, dimred="TSNE", colour_by ="cluster", text_by = 'cluster')+ theme(text = element_text(size = 15), axis.text=element_text(size=15))+geom_density_2d(colour = "black")


# bargraph
distribution <- as.data.frame(table(merged$cluster, merged$orig.ident))
colnames(distribution) <- c('cluster', 'sample', 'frequency')
ggplot(distribution, aes(x=sample, y=frequency, fill=cluster))+geom_bar(position='fill', stat='identity', color='black')+scale_fill_tableau(palette='Classic 10')


### DIFFERENTIAL GENE EXPRESSION
markers<- findMarkers(merged, test="t", groups=merged$cluster, direction='up')
top5=list()
for (i in 1:length(markers)){
  top5[[i]] <- rownames(markers[[i]])[markers[[i]]$Top<=5]
}
p <- plotDots(merged, features=unique(unlist(top5)), group='cluster')
p <- p$data[,c(1,2,4)]
p2 <- dcast(p,Feature ~Group,  value.var='Average')
row.names(p2) <- p2$Feature
p2 <- p2[,-1]
pheatmap(p2, scale="row", cluster_cols=FALSE, cluster_rows = TRUE, cellwidth=5, cellheight = 5, fontsize_row=5, fontsize_col=7)



### zoom in on CD56bright NK cells
bright <- merged[,merged$cluster=='3']
#8211 cels

##normalize again
#HIGHLY VARIABLE GENES #
dec<- modelGeneVar(bright, block=bright$orig.ident, BPPARAM=SnowParam(workers=19))
hvg.sce.var <- getTopHVGs(dec, n=2000)

# PCA (set to scale=TRUE for removal of clusters driven by highly expressed genes)
set.seed(100)
reducedDim(bright, 'PCA') <- NULL
reducedDim(bright, 'UMAP') <- NULL

bright <- runPCA(bright, subset_row=hvg.sce.var, name="PCA", exprs_values='logcounts', scale=TRUE)
bright <- RunHarmony(bright, group.by.vars=c('Dataset', 'Chemistry'), 
                     theta = NULL, lambda = NULL, sigma = 0.1,
                     nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                     max.iter.cluster = 20, epsilon.cluster = 1e-05,
                     epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                     reference_values = NULL, reduction.save = "HARMONY")

g <- buildSNNGraph(bright, k=20, use.dimred = 'HARMONY', type='jaccard')
clustering <- igraph::cluster_louvain(g)$membership
bright$cluster <- factor(clustering)


bright <- runUMAP(bright, dimred="HARMONY",
                  name='UMAP',
                  ncomponents=2,
                  external_neighbors=TRUE, 
                  n_neighbors=30,
                  BNPARAM=AnnoyParam(),
                  BPPARAM=SnowParam(workers=18),
                  n_threads=bpnworkers(SnowParam(workers=18)),
                  min_dist=0.4)

bright <- runUMAP(bright, dimred="HARMONY",
                  name='UMAP_mindis1',
                  ncomponents=2,
                  external_neighbors=TRUE, 
                  n_neighbors=30,
                  BNPARAM=AnnoyParam(),
                  BPPARAM=SnowParam(workers=18),
                  n_threads=bpnworkers(SnowParam(workers=18)),
                  min_dist=1)

bright <- runTSNE(bright, dimred="HARMONY",
                  name='TSNE',
                  ncomponents=2,
                  external_neighbors=TRUE, 
                  BNPARAM=AnnoyParam(),
                  BPPARAM=SnowParam(workers=18),
                  n_threads=bpnworkers(SnowParam(workers=18)))

plotReducedDim(bright, dimred="UMAP", colour_by ="cluster", text_by = 'cluster')+ theme(text = element_text(size = 15), axis.text=element_text(size=15))+geom_density_2d(colour = "black")
plotReducedDim(bright, dimred="UMAP_mindis1", colour_by ="cluster", text_by = 'cluster')+ theme(text = element_text(size = 15), axis.text=element_text(size=15))+geom_density_2d(colour = "black")



distribution <- as.data.frame(table(bright$cluster, bright$orig.ident))
colnames(distribution) <- c('cluster', 'sample', 'frequency')
ggplot(distribution, aes(x=sample, y=frequency, fill=cluster))+geom_bar(position='fill', stat='identity', color='black')+scale_fill_tableau(palette='Classic 10')
ggplot(distribution, aes(x=cluster, y=frequency, fill=sample))+geom_bar(position='fill', stat='identity', color='black')+scale_fill_tableau(palette='Classic 20')


#Differential gene expression
markers<- findMarkers(bright, test="t", groups=bright$cluster, direction='up')
top5=list()
for (i in 1:length(markers)){
  top5[[i]] <- rownames(markers[[i]])[markers[[i]]$Top<=5]
}
p <- plotDots(bright, features=unique(unlist(top5)), group='cluster')
p <- p$data[,c(1,2,4)]
p2 <- dcast(p,Feature ~Group,  value.var='Average')
row.names(p2) <- p2$Feature
p2 <- p2[,-1]
pheatmap(p2, scale="row", cluster_cols=FALSE, cluster_rows = TRUE, cellwidth=5, cellheight = 4, fontsize_row=5, fontsize_col=7)



### zoom in on CD56dim NK cells
dim <- merged[,merged$cluster%in%c('2', '7', '5')]
#37406 cels

##normalize again
#HIGHLY VARIABLE GENES #
dec<- modelGeneVar(dim, block=dim$orig.ident, BPPARAM=SnowParam(workers=19))
hvg.sce.var <- getTopHVGs(dec, n=2000)

# PCA (set to scale=TRUE for removal of clusters driven by highly expressed genes)
set.seed(100)
reducedDim(dim, 'PCA') <- NULL
reducedDim(dim, 'UMAP') <- NULL

dim <- runPCA(dim, subset_row=hvg.sce.var, name="PCA", exprs_values='logcounts', scale=TRUE)
dim <- RunHarmony(dim, group.by.vars=c('Dataset', 'Chemistry'), 
                     theta = NULL, lambda = NULL, sigma = 0.1,
                     nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                     max.iter.cluster = 20, epsilon.cluster = 1e-05,
                     epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                     reference_values = NULL, reduction.save = "HARMONY")

g <- buildSNNGraph(dim, k=30, use.dimred = 'HARMONY', type='jaccard')
clustering <- igraph::cluster_louvain(g)$membership
dim$cluster <- factor(clustering)


dim <- runUMAP(dim, dimred="HARMONY",
                  name='UMAP',
                  ncomponents=2,
                  external_neighbors=TRUE, 
                  n_neighbors=30,
                  BNPARAM=AnnoyParam(),
                  BPPARAM=SnowParam(workers=18),
                  n_threads=bpnworkers(SnowParam(workers=18)),
                  min_dist=0.4)

dim <- runUMAP(dim, dimred="HARMONY",
                  name='UMAP_mindis1',
                  ncomponents=2,
                  external_neighbors=TRUE, 
                  n_neighbors=30,
                  BNPARAM=AnnoyParam(),
                  BPPARAM=SnowParam(workers=18),
                  n_threads=bpnworkers(SnowParam(workers=18)),
                  min_dist=1)

dim <- runTSNE(dim, dimred="HARMONY",
                  name='TSNE',
                  ncomponents=2,
                  external_neighbors=TRUE, 
                  BNPARAM=AnnoyParam(),
                  BPPARAM=SnowParam(workers=18),
                  n_threads=bpnworkers(SnowParam(workers=18)))

plotReducedDim(dim, dimred="UMAP", colour_by ="cluster", text_by = 'cluster')+ theme(text = element_text(size = 15), axis.text=element_text(size=15))+geom_density_2d(colour = "black")
plotReducedDim(dim, dimred="UMAP_mindis1", colour_by ="cluster", text_by = 'cluster')+ theme(text = element_text(size = 15), axis.text=element_text(size=15))+geom_density_2d(colour = "black")



distribution <- as.data.frame(table(dim$cluster, dim$orig.ident))
colnames(distribution) <- c('cluster', 'sample', 'frequency')
ggplot(distribution, aes(x=sample, y=frequency, fill=cluster))+geom_bar(position='fill', stat='identity', color='black')+scale_fill_tableau(palette='Classic 10')
ggplot(distribution, aes(x=cluster, y=frequency, fill=sample))+geom_bar(position='fill', stat='identity', color='black')+scale_fill_tableau(palette='Classic 20')


#Differential gene expression
markers<- findMarkers(dim, test="t", groups=dim$cluster, direction='up')
top5=list()
for (i in 1:length(markers)){
  top5[[i]] <- rownames(markers[[i]])[markers[[i]]$Top<=10]
}
p <- plotDots(dim, features=unique(unlist(top5)), group='cluster')
p <- p$data[,c(1,2,4)]
p2 <- dcast(p,Feature ~Group,  value.var='Average')
row.names(p2) <- p2$Feature
p2 <- p2[,-1]
pheatmap(p2, scale="row", cluster_cols=FALSE, cluster_rows = TRUE, cellwidth=5, cellheight = 4, fontsize_row=5, fontsize_col=7)


#singleR prediction
pred <- SingleR(test=merged, ref=sce, labels=sce$label, de.method = 'wilcox')
merged$reflabel <- pred$labels
pred.sub <- SingleR(test=merged, ref=sce, labels=sce$sublabel, de.method = 'wilcox')
merged$reflabel.sub <- pred.sub$labels

plotReducedDim(merged, "UMAP",  colour_by="reflabel")+
  scale_color_manual(values = c('1' = "#729ece", '2' = "#67bf5c", '4'='#ed665d', '6'='#ed97ca', '7'='#ff9e4a', '8'='#a2a2a2'),aesthetics = "fill")+ theme(text = element_text(size = 15), axis.text=element_text(size=15))
plotReducedDim(merged, "UMAP",  colour_by="reflabel.sub")+
  scale_color_manual(values = c('1' = "#729ece", '2' = "#67bf5c", '4'='#ed665d', '6'='#ed97ca', '7'='#ff9e4a', '8'='#a2a2a2'),aesthetics = "fill")+ theme(text = element_text(size = 15), axis.text=element_text(size=15))

### zoom in on CD56adaptive NK cells
adaptive <- merged[,merged$cluster%in%c('1', '6', '8')]
#26172 cells

##normalize again?
#HIGHLY VARIABLE GENES #
dec<- modelGeneVar(adaptive, block=adaptive$orig.ident, BPPARAM=SnowParam(workers=19))
hvg.sce.var <- getTopHVGs(dec, n=2000)

# PCA (set to scale=TRUE for removal of clusters driven by highly expressed genes)
set.seed(100)
reducedadaptive(adaptive, 'PCA') <- NULL
reducedadaptive(adaptive, 'UMAP') <- NULL

adaptive <- runPCA(adaptive, subset_row=hvg.sce.var, name="PCA", exprs_values='logcounts', scale=TRUE)
adaptive <- RunHarmony(adaptive, group.by.vars=c('Dataset', 'Chemistry'), 
                  theta = NULL, lambda = NULL, sigma = 0.1,
                  nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                  max.iter.cluster = 20, epsilon.cluster = 1e-05,
                  epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                  reference_values = NULL, reduction.save = "HARMONY")

g <- buildSNNGraph(adaptive, k=20, use.dimred = 'HARMONY', type='jaccard')
clustering <- igraph::cluster_louvain(g)$membership
adaptive$cluster_2 <- factor(clustering)


adaptive <- runUMAP(adaptive, dimred="HARMONY",
               name='UMAP',
               ncomponents=2,
               external_neighbors=TRUE, 
               n_neighbors=30,
               BNPARAM=AnnoyParam(),
               BPPARAM=SnowParam(workers=18),
               n_threads=bpnworkers(SnowParam(workers=18)),
               min_dist=0.4)

adaptive <- runUMAP(adaptive, dimred="HARMONY",
               name='UMAP_mindis1',
               ncomponents=2,
               external_neighbors=TRUE, 
               n_neighbors=30,
               BNPARAM=AnnoyParam(),
               BPPARAM=SnowParam(workers=18),
               n_threads=bpnworkers(SnowParam(workers=18)),
               min_dist=1)

adaptive <- runTSNE(adaptive, adaptivered="HARMONY",
               name='TSNE',
               ncomponents=2,
               external_neighbors=TRUE, 
               BNPARAM=AnnoyParam(),
               BPPARAM=SnowParam(workers=18),
               n_threads=bpnworkers(SnowParam(workers=18)))

plotReducedDim(adaptive, dimred="UMAP", colour_by ="cluster", text_by = 'cluster')+ theme(text = element_text(size = 15), axis.text=element_text(size=15))+geom_density_2d(colour = "black")
plotReducedDim(adaptive, dimred="UMAP_mindis1", colour_by ="cluster", text_by = 'cluster')+ theme(text = element_text(size = 15), axis.text=element_text(size=15))+geom_density_2d(colour = "black")



distribution <- as.data.frame(table(adaptive$cluster, adaptive$orig.ident))
colnames(distribution) <- c('cluster', 'sample', 'frequency')
ggplot(distribution, aes(x=sample, y=frequency, fill=cluster))+geom_bar(position='fill', stat='identity', color='black')+scale_fill_tableau(palette='Classic 10')
ggplot(distribution, aes(x=cluster, y=frequency, fill=sample))+geom_bar(position='fill', stat='identity', color='black')+scale_fill_tableau(palette='Classic 20')


#Differential gene expression
markers<- findMarkers(adaptive, test="t", groups=adaptive$cluster, direction='up')
top5=list()
for (i in 1:length(markers)){
  top5[[i]] <- rownames(markers[[i]])[markers[[i]]$Top<=5]
}
p <- plotDots(adaptive, features=unique(unlist(top5)), group='cluster')
p <- p$data[,c(1,2,4)]
p2 <- dcast(p,Feature ~Group,  value.var='Average')
row.names(p2) <- p2$Feature
p2 <- p2[,-1]
pheatmap(p2, scale="row", cluster_cols=FALSE, cluster_rows = TRUE, cellwidth=5, cellheight = 4, fontsize_row=5, fontsize_col=7)

##conclusion: adaptive cluster 6 can be subclustered into two populations




### NKG2C HASHTAG ###

##For GSE19703, cells were sorted into a NKG2C+ and NKG2C- fraction and labeled with hashtags 
# import csv with hashtag information and label cells accordingly
# Lib1-5 can be split into two donors using the hashtag information

CMVneg1_hashtags <- read.csv("path/hashtags/CMVneg1_hashtags.csv")
CMVneg2_hashtags <- read.csv("path/hashtags/CMVneg2_hashtags.csv")
CMVpos4_hashtags <- read.csv("path/hashtags/CMVpos4_hashtags.csv")
CMVpos3_hashtags <- read.csv("path/hashtags/CMVpos3_hashtags.csv")
CMVpos2_hashtags <- read.csv("path/hashtags/CMVpos2_hashtags.csv")
CMVpos1_5_LibB_hashtags <- read.csv("path/hashtags/CMVpos1_5_LibB_hashtags.csv")
CMVpos1_5_LibA_hashtags <- read.csv("path/hashtags/CMVpos1_5_LibA_hashtags.csv")

CMVneg1_hashtags$orig.ident <- 'CMVneg1_donorA'
CMVneg1_hashtags$Barcode <- paste(CMVneg1_hashtags$X, "-1", sep="")
CMVneg2_hashtags$orig.ident <- 'CMVneg2_donorB'
CMVneg2_hashtags$Barcode <- paste(CMVneg2_hashtags$X, "-1", sep="")
CMVpos1_5_LibA_hashtags$orig.ident <- 'CMVpos1_5_libA'
CMVpos1_5_LibA_hashtags$Barcode <- paste(CMVpos1_5_LibA_hashtags$X, "-1", sep="")
CMVpos1_5_LibB_hashtags$orig.ident <- 'CMVpos1_5_libB'
CMVpos1_5_LibB_hashtags$Barcode <- paste(CMVpos1_5_LibB_hashtags$X, "-1", sep="")
CMVpos2_hashtags$orig.ident <- 'CMVpos2_donorC'
CMVpos2_hashtags$Barcode <- paste(CMVpos2_hashtags$X, "-1", sep="")
CMVpos3_hashtags$orig.ident <- 'CMVpos3_donorD'
CMVpos3_hashtags$Barcode <- paste(CMVpos3_hashtags$X, "-1", sep="")
CMVpos4_hashtags$orig.ident <- 'CMVpos4_donorE'
CMVpos4_hashtags$Barcode <- paste(CMVpos4_hashtags$X, "-1", sep="")

hashtags <- do.call(rbind,(list(CMVneg1_hashtags, CMVneg2_hashtags, CMVpos1_5_LibA_hashtags, CMVpos1_5_LibB_hashtags, CMVpos2_hashtags, CMVpos3_hashtags, CMVpos4_hashtags)))
df<- cbind(colnames(merged), as.data.frame(merged$orig.ident))
df$dataset <- 'filtered'
colnames(df) <- c('Barcode', 'orig.ident', 'dataset')
match_hashtag <- merge(df, hashtags, all.x=TRUE, all.y=TRUE)
match_hashtag <- match_hashtag%>%filter(match_hashtag$dataset =='filtered')

#rownames merged should match barcodes in colnames sce object
match_hashtag$code <- paste(match_hashtag$Barcode, match_hashtag$orig.ident, sep="_")
merged$code <- paste(colnames(merged), merged$orig.ident, sep="_")
match_hashtag <- match_hashtag[match(merged$code, match_hashtag$code),]
merged$nkg2c <- match_hashtag$x
merged$nkg2c <-as.character(merged$nkg2c)
merged$nkg2c[is.na(merged$nkg2c)] <- 'unknown'

table(merged$orig.ident, merged$nkg2c)
merged$donor <- as.character(merged$orig.ident)
merged$donor[merged$nkg2c == 'CMVpos1-NKG2Cneg' | merged$nkg2c == 'CMVpos1-NKG2Cpos'  ] <- "CMVpos1"
merged$donor[merged$nkg2c == 'CMVpos5-NKG2Cneg' | merged$nkg2c == 'CMVpos5-NKG2Cpos'  ] <- "CMVpos5"
merged$nkg2c[merged$nkg2c == 'CMVpos1-NKG2Cneg' | merged$nkg2c == 'CMVpos5-NKG2Cneg'  ] <- "NKG2Cneg"
merged$nkg2c[merged$nkg2c == 'CMVpos1-NKG2Cpos' | merged$nkg2c == 'CMVpos5-NKG2Cpos'  ] <- "NKG2Cpos"
table(merged$donor, merged$nkg2c)



