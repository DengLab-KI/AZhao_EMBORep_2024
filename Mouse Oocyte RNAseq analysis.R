#load packages needed

library(clusterProfiler)
library(GOSemSim)
library(magrittr)
library(AnnotationHub)
library(dplyr)
library(DOSE)
library(enrichplot)
library(ggnewscale)
library(Seurat)
library(DESeq2)
library(sva)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(MAST)

#load files for 
oocyte_counts <- readRDS("~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Year 1/STZ project/oocyte/scRNAseq/Smartseq3_Oocyte.dgecounts.rds")

### For Figure 2c ### 

#extract introns and exons
oocyte_counts_2<-oocyte_counts$umicount$inex$all

#import sample info
oocyte_sampleinfo<-read.table("RNA_oocyte_sampleinfo_2.txt", header=TRUE)

#import gene names
oocyte_genenames<-read.table("Smartseq3_Oocyte.gene_names.txt", header=TRUE)

#merge files
df.oocyte<-merge(oocyte_counts_2, oocyte_genenames, by.x=0, by.y='gene_id')

#remove ensmus
df.oocyte.2<-dplyr::select(df.oocyte, -c(Row.names))

#add all transcript variants
df.oocyte.3<-aggregate(df.oocyte.2[,1:272], list(df.oocyte.2$gene_name), sum)

#put gene names as row names
rnames.tmp <- df.oocyte.3$Group.1
df.oocyte.3$Group.1 <- NULL
rownames(df.oocyte.3) <- rnames.tmp

#make df to df.O1
df.O1<-df.oocyte.3

#order
df.O1<-df.O1[,order(colnames(df.O1))]

#fix sampleinfo
rownames(oocyte_sampleinfo)<-oocyte_sampleinfo$Barcodes
oocyte_sampleinfo<-oocyte_sampleinfo[order(row.names(oocyte_sampleinfo)),]

#check
rownames(oocyte_sampleinfo)==colnames(df.O1)

#Create seurat object

seurat.obj.oocyte<-CreateSeuratObject(counts=df.O1, min.cells=3, min.features=1, project="STZ_oocyte")
seurat.obj.oocyte<-AddMetaData(object=seurat.obj.oocyte, metadata=oocyte_sampleinfo)

#Select samples, filtering
selected_c <- WhichCells(seurat.obj.oocyte, expression = nFeature_RNA >  500 & percent_mito<0.15 )
selected_f <- rownames(seurat.obj.oocyte)[ Matrix::rowSums(seurat.obj.oocyte) > 1]
data.filt <- subset(seurat.obj.oocyte, features=selected_f, cells=selected_c)
dim(data.filt)
oocyte.obj<- data.filt

#Normalize
oocyte.obj <- NormalizeData(object = oocyte.obj, normalization.method = "LogNormalize", scale.factor = 10000)
oocyte.obj <- ScaleData(object = oocyte.obj)
oocyte.obj <- FindVariableFeatures(object = oocyte.obj, selection.method="vst", mean.function = ExpMean, 
                                   dispersion.function = LogVMR, x.low.cutoff = 0.0125, nfeatures = 2000,
                                   x.high.cutoff = 5, y.cutoff = 0.5) 

#Differentialexpression
Idents(oocyte.obj) <- "Condition"
condition.diffgenes <- FindMarkers(oocyte.obj, ident.1 = "STZ", ident.2="Control", min.pct=0.01, logfc.threshold=0.05, test.use="wilcox")

#res.1 - Ctrl vs STZ
res.1<-condition.diffgenes

#enhancedVolcano 1
EnhancedVolcano(res.1,
                lab = rownames(res.1),
                x = 'avg_log2FC',
                col=c('black', 'black', 'black', 'red3'),
                xlab=substitute(paste(bold(Log[2] ~ "fold change"))),
                ylab = substitute(paste(bold(~-Log[10] ~ "P-value"))),
                legendLabels=c('Not Significant', '', '', ''),
                title="Control vs STZ oocytes",
                y = 'p_val_adj', xlim=c(-5,5), ylim=c(0,15), FCcutoff=1, pCutoff=0.05, pointSize=3, labSize=6, drawConnectors=TRUE, widthConnectors=0.6)

### For Figure 2b, e, EV2

#extract introns and exons
oocyte_counts_2<-oocyte_counts$readcount$inex$all

#import sample info
oocyte_sampleinfo<-read.table("RNA_oocyte_sampleinfo_2.txt", header=TRUE)

#import gene names
oocyte_genenames<-read.table("Smartseq3_Oocyte.gene_names.txt", header=TRUE)

#merge files
df.oocyte<-merge(oocyte_counts_2, oocyte_genenames, by.x=0, by.y='gene_id')

#remove ensmus
df.oocyte.2<-dplyr::select(df.oocyte, -c(Row.names))

#add all transcript variants
df.oocyte.3<-aggregate(df.oocyte.2[,1:272], list(df.oocyte.2$gene_name), sum)

#put gene names as row names
rnames.tmp <- df.oocyte.3$Group.1
df.oocyte.3$Group.1 <- NULL
rownames(df.oocyte.3) <- rnames.tmp

#make df to df.O1
df.O1<-df.oocyte.3

#order
df.O1<-df.O1[,order(colnames(df.O1))]

#fix sampleinfo
rownames(oocyte_sampleinfo)<-oocyte_sampleinfo$Barcodes
oocyte_sampleinfo<-oocyte_sampleinfo[order(row.names(oocyte_sampleinfo)),]

#check
rownames(oocyte_sampleinfo)==colnames(df.O1)

#Create seurat object

seurat.obj.oocyte<-CreateSeuratObject(counts=df.O1, min.cells=3, min.features=1, project="STZ_oocyte")
seurat.obj.oocyte<-AddMetaData(object=seurat.obj.oocyte, metadata=oocyte_sampleinfo)

#Count mito genes
seurat.obj.oocyte[["percent.mt"]] <- PercentageFeatureSet(seurat.obj.oocyte, pattern = "^mt-")
mito_genes <- rownames(seurat.obj.oocyte)[grep("mt-",rownames(seurat.obj.oocyte))]
head(mito_genes, 10)

#Count ribo genes
seurat.obj.oocyte<- PercentageFeatureSet(seurat.obj.oocyte, "rp[sl]", col.name = "percent_ribo")
ribo_genes <- rownames(seurat.obj.oocyte)[grep("rp[sl]",rownames(seurat.obj.oocyte))]
head(ribo_genes,10)
seurat.obj.oocyte$percent_ribo <- colSums(seurat.obj.oocyte@assays$RNA@counts[ribo_genes,]  ) / total_counts_per_cell*100

#Plot
VlnPlot(seurat.obj.oocyte, features = "nFeature_RNA", group.by = "Condition", cols=c("white","blue3")) + geom_boxplot()
VlnPlot(seurat.obj.oocyte, features = "nCount_RNA", group.by = "Condition", cols=c("white","blue3")) + geom_boxplot()
VlnPlot(seurat.obj.oocyte, features = "percent.mt", group.by = "Condition", cols=c("white","blue3")) + geom_boxplot()
VlnPlot(seurat.obj.oocyte, features = "percent_ribo", group.by = "Condition", cols=c("white","blue3")) + geom_boxplot()

#Select samples, filtering
selected_c <- WhichCells(seurat.obj.oocyte, expression = nFeature_RNA >  10000 & percent.mt<15 )
selected_f <- rownames(seurat.obj.oocyte)[ Matrix::rowSums(seurat.obj.oocyte) > 1]
data.filt <- subset(seurat.obj.oocyte, features=selected_f, cells=selected_c)
dim(data.filt)
oocyte.obj<- data.filt

#Normalize
oocyte.obj <- NormalizeData(object = oocyte.obj, normalization.method = "LogNormalize", scale.factor = 10000)
oocyte.obj <- ScaleData(object = oocyte.obj)
oocyte.obj <- FindVariableFeatures(object = oocyte.obj, selection.method="vst", mean.function = ExpMean, 
                                   dispersion.function = LogVMR, x.low.cutoff = 0.0125, nfeatures = 2000,
                                   x.high.cutoff = 5, y.cutoff = 0.5) 

#Run PCA
oocyte.obj <- RunPCA(object = oocyte.obj, features=VariableFeatures(object=oocyte.obj))
x<-oocyte.obj@reductions$pca@cell.embeddings[,1]
y<-oocyte.obj@reductions$pca@cell.embeddings[,2]

#calculate PC1 and PC2 percentage
stdev<-oocyte.obj@reductions$pca@stdev
var<-stdev^2
percent<-var*100/sum(var)
percent

#Plot PCA
z<-DimPlot(object = oocyte.obj, 
           reduction = "pca", 
           group.by = 'Condition', 
           pt.size = 2, combine=TRUE, cols=c("red3","blue3"), shape.by="Mouse")+scale_shape_manual(values=1:14)
z+xlab('PC1 16.3%')+ylab('PC2 5.4%')


#VlnPlots
oocyte_markers<-c("Dppa3","Gdf9","Bmp15","Mos","Zp3","Kit","Nlrp5","Ddx4")
VlnPlot(oocyte.obj, features=oocyte_markers, group.by="Condition", cols=c("white","blue3")) + geom_boxplot()

### For Figure 2g and h ###

#Aggregate expression of our samples
pseudobulk<-Seurat:::PseudobulkExpression(object = oocyte.obj, pb.method = 'aggregate', slot = 'counts', group.by="Mouse")
head(pseudobulk)
str(pseudobulk)
pseudo1<-as.data.frame(pseudobulk)
head(pseudo1)
colnames(pseudo1)<-c('Control_1','Control_2','Control_3','Control_4','Control_5','STZ_1','STZ_2','STZ_3','STZ_4','STZ_5','STZ_6','STZ_7','STZ_8','STZ_9')
#make sampleinfo for pseudo
pseudo_sampleinfo<-c('Control_1','Control_2','Control_3','Control_4','Control_5','STZ_1','STZ_2','STZ_3','STZ_4','STZ_5','STZ_6','STZ_7','STZ_8','STZ_9')
pseudo_sampleinfo<-as.data.frame(pseudo_sampleinfo)
pseudo_sampleinfo$Condition<-factor(c(rep("Control",5),rep("STZ",9)))
head(pseudo_sampleinfo)
rownames(pseudo_sampleinfo)<-pseudo_sampleinfo$pseudo_sampleinfo

#Create seurat object pseudo

seurat.obj.pseudo<-CreateSeuratObject(counts=pseudo1, min.cells=3, min.features=200, project="STZ_oocyte_pseudo")
seurat.obj.pseudo<-AddMetaData(object=seurat.obj.pseudo, metadata=pseudo_sampleinfo)

#Count mito genes
seurat.obj.pseudo[["percent.mt"]] <- PercentageFeatureSet(seurat.obj.pseudo, pattern = "^mt-")
mito_genes <- rownames(seurat.obj.pseudo)[grep("mt-",rownames(seurat.obj.pseudo))]
head(mito_genes, 10)

#Plot mito_genes
total_counts_per_cell = colSums(as.matrix( seurat.obj.pseudo@assays$RNA@counts))
seurat.obj.pseudo$percent_mito <- colSums( as.matrix( seurat.obj.pseudo@assays$RNA@counts[mito_genes,])  ) / total_counts_per_cell
seurat.obj.pseudo$total_counts_per_cell = total_counts_per_cell

#Select samples, filtering away samples of lower counts
selected_c <- WhichCells(seurat.obj.pseudo, expression = nFeature_RNA >  19000 & percent_mito<0.15 )
selected_f <- rownames(seurat.obj.pseudo)[ Matrix::rowSums(seurat.obj.pseudo) > 1]
data.filt <- subset(seurat.obj.pseudo, features=selected_f, cells=selected_c)
dim(data.filt)
pseudo.obj<- data.filt

# write table with counts
write.table(data.filt@assays$RNA@counts,"oocytePseudobulk.txt", row.names = T, sep = "\t")

### import file and clean from TET3 paper for integration ###
counts<-read_excel("41586_2022_4756_MOESM3_ESM.xlsx")

colnames(counts)
counts1<-counts[,c(2,10:15)]
colnames(counts1)
counts1<-as.data.frame(counts1)
head(counts1)
counts1<-counts1[!duplicated(data.frame(counts1[,1])),]
length(unique(counts1$genes))==nrow(counts1)

#sampleinfo TET3 paper
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

#merge with my data
pseudobulk<-read.table("oocytePseudobulk.txt")
head(pseudobulk)
pseudobulk$geneID<-rownames(pseudobulk)
head(counts1)
counts1$genes<-rownames(counts1)
dim(counts1)
dim(pseudobulk)
merged_file<-merge(counts1, pseudobulk, by.x="genes",by.y="geneID")

dim(merged_file)

colnames(merged_file)
rownames(merged_file)
rownames(merged_file)<-merged_file$genes

merged_file2<-merged_file
merged_file2$genes<-NULL
merged_file2<-as.data.frame(apply(merged_file2, 2, as.numeric))
str(merged_file2)
rownames(merged_file2)<-merged_file$genes


sampleinfo_merged<-colnames(merged_file2)
sampleinfo_merged<-as.data.frame(sampleinfo_merged)
head(sampleinfo_merged)
sampleinfo_merged$treatment<-c(rep('Control',3), rep('STZ_Chen et al.',3), rep('Control',5), rep('STZ This paper.',7))

sampleinfo_merged$Paper<-c(rep('Chen et al.',6), rep('This paper.', 12))
rownames(sampleinfo_merged)<-sampleinfo_merged$sampleinfo_merged

#Create Seurat object with merged data
seurat.obj<-CreateSeuratObject(counts=merged_file2, min.cells=3, min.features=200, project="tet3")
seurat.obj<-AddMetaData(object=seurat.obj, metadata=sampleinfo_merged)

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
selected_c <- WhichCells(seurat.obj, expression = nFeature_RNA >  10000 &nCount_RNA >40000  & percent_mito<0.1 )
selected_f <- rownames(seurat.obj)[ Matrix::rowSums(seurat.obj) > 1]
data.filt <- subset(seurat.obj, features=selected_f, cells=selected_c)
dim(data.filt)
tet3.obj <- data.filt

#batch effect correlation using ComBat-seq
combat_matrix<-as.matrix(tet3.obj@assays$RNA@counts)
combat_dat<-ComBat_seq(combat_matrix, tet3.obj$Paper, group=tet3.obj$treatment)
tet3.obj@assays$RNA@counts<-combat_dat
str(tet3.obj@assays$RNA@counts)
head(tet3.obj@assays$RNA@counts)

#Normalize
tet3.obj <- NormalizeData(object = tet3.obj, normalization.method = "LogNormalize", scale.factor = 10000)
tet3.obj <- ScaleData(object = tet3.obj)
tet3.obj <- FindVariableFeatures(object = tet3.obj, mean.function = ExpMean, 
                                 dispersion.function = LogVMR, x.low.cutoff = 0.0125, nfeatures = 2000,
                                 x.high.cutoff = 5, y.cutoff = 0.5) 

#Run PCA
tet3.obj <- RunPCA(object = tet3.obj, pc.genes = tet3.obj@var.genes,npcs = 17,
                   do.print = TRUE, pcs.print = 1:5,genes.print = 3)
x<-tet3.obj@reductions$pca@cell.embeddings[,1]
y<-tet3.obj@reductions$pca@cell.embeddings[,2]

#calculate PC1 and PC2 percentage
stdev<-tet3.obj@reductions$pca@stdev
var<-stdev^2
percent<-var*100/sum(var)
percent

#Plot PCA
z<-DimPlot(object = tet3.obj, reduction = "pca", group.by = 'treatment', shape.by='Paper',pt.size = 2.5, combine=TRUE, cols=c("red3","orange3","blue3"))
z+xlab('PC1 13.7%')+ylab('PC2 12.4%')+geom_label(aes(x,y,label=tet3.obj$sampleinfo_merged))
z+xlab('PC1 29.7%')+ylab('PC2 18.7%')

#VlnPlots
VlnPlot(tet3.obj, features="Tet3", group.by='treatment', cols=c("red3","orange3","blue3"))+geom_boxplot()

plots <- VlnPlot(tet3.obj, c("Tet3","Gdf9","Zp3"), group.by='treatment', cols=c("red3","orange3","blue3"), ncol=3, fill='treatment', combine=F)
for (i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot() + theme(legend.position = "none")
}
CombinePlots(plots)

