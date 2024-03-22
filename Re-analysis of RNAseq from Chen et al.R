#load packages needed
library(dplyr)
library(Seurat)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(readxl)

#import file and clean
counts<-read_excel("41586_2022_4756_MOESM3_ESM.xlsx")

colnames(counts)
counts1<-counts[,c(2,10:15)]
colnames(counts1)
counts1<-as.data.frame(counts1)
head(counts1)
counts1<-counts1[!duplicated(data.frame(counts1[,1])),]
length(unique(counts1$genes))==nrow(counts1)

#sampleinfo
sampleinfo_tet3<-c('Ctl1_count',"Ctrl_2count","Ctrl3_count","HG1_count","HG2_count","HG3_count")
sampleinfo_tet3
sampleinfo_tet3<-as.data.frame(sampleinfo_tet3)
sampleinfo_tet3$treatment<-c('control','control','control','diabetes','diabetes','diabetes')
row.names(sampleinfo_tet3)<-sampleinfo_tet3$sampleinfo_tet3

#put gene names as row names
rnames.tet3<-counts1$genes
rownames(counts1)<-counts1$genes
counts1$genes<-NULL
counts1
row.names(sampleinfo_tet3)==colnames(counts1)

#Create Seurat object
seurat.obj<-CreateSeuratObject(counts=counts1, min.cells=3, min.features=200, project="tet3")
seurat.obj<-AddMetaData(object=seurat.obj, metadata=sampleinfo_tet3)

#Count mito genes
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
mito_genes <- rownames(seurat.obj)[grep("mt-",rownames(seurat.obj))]
head(mito_genes,10)

#Plot mito_genes
total_counts_per_cell = colSums(as.matrix( seurat.obj@assays$RNA@counts))
seurat.obj$percent_mito <- colSums( as.matrix( seurat.obj@assays$RNA@counts[mito_genes,])  ) / total_counts_per_cell
seurat.obj$total_counts_per_cell = total_counts_per_cell
VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, group.by = "treatment")

#Select samples, filtering
selected_c <- WhichCells(seurat.obj, expression = nFeature_RNA >  3000 &nCount_RNA >40000  & percent_mito<0.1 )
selected_f <- rownames(seurat.obj)[ Matrix::rowSums(seurat.obj) > 1]
data.filt <- subset(seurat.obj, features=selected_f, cells=selected_c)
dim(data.filt)
tet3.obj <- data.filt

#Normalize
tet3.obj <- NormalizeData(object = tet3.obj, normalization.method = "LogNormalize", scale.factor = 10000)
tet3.obj <- ScaleData(object = tet3.obj)
tet3.obj <- FindVariableFeatures(object = tet3.obj, mean.function = ExpMean, 
                                dispersion.function = LogVMR, x.low.cutoff = 0.0125, nfeatures = 2000,
                                x.high.cutoff = 5, y.cutoff = 0.5) 

# write table with counts
write.table(data.filt@assays$RNA@counts,"counts.txt", row.names = T, sep = "\t")
y = read.table("counts.txt", header = T, row.names = 1)
colnames(y)

#DESEQ2
dds<-DESeqDataSetFromMatrix(as.matrix(y), colData=sampleinfo_tet3, design= ~ treatment)
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)

head(res)
dim(res)

#enhancedvolcano
EnhancedVolcano(res,
                lab = rownames(res),
                selectLab=c("Tet3"),
                x = 'log2FoldChange',
                title="Chen et al. Reanalysis 717 DEGs",
                legendLabels=c('Not Significant', '', '', 'Significant'),
                xlab=substitute(paste(bold(
                  -Log[2] ~ "fold change"))),
                ylab = substitute(paste(bold(~-Log[10] ~ "P-value"))),
                col=c('black','black','black','red'),
                y = 'padj', xlim=c(-5,5), ylim=c(0,15), FCcutoff=1, pCutoff=0.05, pointSize=3, labSize=5, drawConnectors=TRUE, widthConnectors=2)



       