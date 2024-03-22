#load packages needed

library(clusterProfiler)
library(GOSemSim)
library(magrittr)
library(org.Mm.eg.db)
library(AnnotationHub)
library(dplyr)
library(DOSE)
library(enrichplot)
library(ggnewscale)
library(Seurat)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(MAST)


#load files
hum_oocyte_counts <- readRDS("~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Year 2/RNA/human oocytes/Smartseq3_Oocyte.dgecounts.rds")

#import gene names
oocyte_genenames<-read.table("Smartseq3_Oocyte.gene_names.txt", header=TRUE)

### extract introns and exons ##
hum_oocyte_counts_2<-hum_oocyte_counts$readcount$inex$all

#merge files
df.oocyte_hum<-merge(hum_oocyte_counts_2, oocyte_genenames, by.x=0, by.y='gene_id')

#remove ensmus
df.oocyte_hum.2<-dplyr::select(df.oocyte_hum, -c(Row.names))

#add all transcript variants
df.oocyte_hum.3<-aggregate(df.oocyte_hum.2[,1:66], list(df.oocyte_hum.2$gene_name), sum)

#put gene names as row names
rnames.tmp <- df.oocyte_hum.3$Group.1
df.oocyte_hum.3$Group.1 <- NULL
rownames(df.oocyte_hum.3) <- rnames.tmp

#make df to df.hum
df.hum<-df.oocyte_hum.3

#import sample info
oocyte_sampleinfo_hum<-read.table("sampleinfo_humanoocyte.txt", header=TRUE)

#fix sampleinfo
rownames(oocyte_sampleinfo_hum)<-oocyte_sampleinfo_hum$Barcode
oocyte_sampleinfo_hum<-oocyte_sampleinfo_hum[order(row.names(oocyte_sampleinfo_hum)),]

#order
df.hum<-df.hum[,order(colnames(df.hum))]

#sampleinfo 2, minus empty well
oocyte_sampleinfo_hum2<-oocyte_sampleinfo_hum[-12,]

#check again
rownames(oocyte_sampleinfo_hum2)==colnames(df.hum)

#remove 10mM cells
oocyte_sampleinfo_hum3<-oocyte_sampleinfo_hum2[-c(2,5,12,16,17,23,33,36,38,40,45,47,51,52,53,56,58,60,61,62,65,66),]
df.hum2<-df.hum[,-c(2,5,12,16,17,23,33,36,38,40,45,47,51,52,53,56,58,60,61,62,65,66)]

rownames(oocyte_sampleinfo_hum3)==colnames(df.hum2)

#Create seurat object
seurat.obj.oocyte_hum<-CreateSeuratObject(counts=df.hum2, min.cells=3, min.features=500, project="Human oocyte")
seurat.obj.oocyte_hum<-AddMetaData(object=seurat.obj.oocyte_hum, metadata=oocyte_sampleinfo_hum3)

#Count mito genes
seurat.obj.oocyte_hum[["percent.mt"]] <- PercentageFeatureSet(seurat.obj.oocyte_hum, pattern = "^MT-")

seurat.obj.oocyte_hum$Condition<-factor(seurat.obj.oocyte_hum$Condition, levels=c("Control","5mM"))
seurat.obj.oocyte_hum$Condition

#Count ribo genes
seurat.obj.oocyte_hum<- PercentageFeatureSet(seurat.obj.oocyte_hum, "RP[SL]", col.name = "percent_ribo")
ribo_genes <- rownames(seurat.obj.oocyte_hum)[grep("RP[SL]",rownames(seurat.obj.oocyte_hum))]
head(ribo_genes,10)
total_counts_per_cell = colSums(as.matrix( seurat.obj.oocyte_hum@assays$RNA@counts))
seurat.obj.oocyte_hum$percent_ribo <- colSums(seurat.obj.oocyte_hum@assays$RNA@counts[ribo_genes,]  ) / total_counts_per_cell*100

#Plot 
VlnPlot(seurat.obj.oocyte_hum, features = "nFeature_RNA", group.by="Condition", cols=c("white","blue3"))+geom_boxplot()
VlnPlot(seurat.obj.oocyte_hum, features = "nCount_RNA", group.by="Condition",cols=c("white","blue3"))+geom_boxplot()
VlnPlot(seurat.obj.oocyte_hum, features = "percent.mt", group.by="Condition",cols=c("white","blue3"))+geom_boxplot()
VlnPlot(seurat.obj.oocyte_hum, features = "percent_ribo", group.by = "Condition", cols=c("white","blue3")) + geom_boxplot()


#Select samples, filtering
selected_c <- WhichCells(seurat.obj.oocyte_hum, expression = nFeature_RNA >  1000 & percent.mt<25)
selected_f <- rownames(seurat.obj.oocyte_hum)[ Matrix::rowSums(seurat.obj.oocyte_hum) > 1]
data.filt<- subset(seurat.obj.oocyte_hum, features=selected_f, cells=selected_c)
dim(data.filt)
oocyte.obj_hum<- data.filt



#Normalize
oocyte.obj_hum <- NormalizeData(object = oocyte.obj_hum, normalization.method = "LogNormalize", scale.factor = 10000)
oocyte.obj_hum <- ScaleData(object = oocyte.obj_hum)
oocyte.obj_hum <- FindVariableFeatures(object = oocyte.obj_hum, selection.method="vst", mean.function = ExpMean, 
                                       dispersion.function = LogVMR, x.low.cutoff = 0.0125, nfeatures = 2000,
                                       x.high.cutoff = 5, y.cutoff = 0.5) 

#Run PCA
oocyte.obj_hum <- RunPCA(object = oocyte.obj_hum, features=VariableFeatures(object=oocyte.obj_hum), npcs=32)
x<-oocyte.obj_hum@reductions$pca@cell.embeddings[,1]
y<-oocyte.obj_hum@reductions$pca@cell.embeddings[,2]

#calculate PC1 and PC2 percentage
stdev<-oocyte.obj_hum@reductions$pca@stdev
var<-stdev^2
percent<-var*100/sum(var)
percent

#Plot PCA
oocyte.obj_hum$Patient
oocyte.obj_hum$patient<-c("2","22","6","7","17","2","23","16","2","23","10","18","18",
                              "9","21","2","7","6","21",
                              "23","11","22","23","23","21","22","1","1","6","6","22","11","11")
z<-DimPlot(object = oocyte.obj_hum, reduction = "pca",
           group.by = 'Condition', pt.size = 3, combine=TRUE, 
           cols=c("red3","blue3"),
           shape.by="patient")+scale_shape_manual(values=1:23)
z+xlab('PC1 9.1%')+ylab('PC2 6.3%')

#TET3 violinplot
VlnPlot(oocyte.obj_hum, features="TET3",group.by="Condition",cols=c("red","blue"))+geom_boxplot()

#VlnPlots
oocyte_markers_hum<-c("DPPA3","GDF9","BMP15","MOS","ZP3","KIT","NLRP5","DDX4")
VlnPlot(oocyte.obj_hum, features=oocyte_markers_hum, group.by="Condition", cols=c("white","blue3")) + geom_boxplot()
dplot<-DotPlot(oocyte.obj_hum, features=oocyte_markers_hum, cols=c("blue","red3"), group.by="Condition",scale.by="size")
dplot+scale_color_gradientn(name = "Average Expression",
                         colours = c("lightgrey","red3"), 
                         values = rescale(c(-5,5)),
                         limits=c(-5,5))+scale_size_continuous(
                           breaks = c(0,25,50,100),
                           limits = c(0, 100))

#DE
Idents(oocyte.obj_hum) <- "Condition"
condition.diffgenes_Cvs5 <- FindMarkers(oocyte.obj_hum, ident.1 = "5mM", ident.2="Control", min.pct=0.01, logfc.threshold=0.05, test.use="MAST")
condition.diffgenes_Cvs5


#C vs 5#

EnhancedVolcano(condition.diffgenes_Cvs5, x='avg_log2FC', y='p_val_adj',
                xlim=c(-5,5),
                ylim=c(0,15),
                lab=rownames(condition.diffgenes_Cvs5),
                xlab=substitute(paste(bold(-Log[2] ~ "fold change"))),
                ylab = substitute(paste(bold(~-Log[10] ~ "P-value"))),
                FCcutoff=1,
                col=c("black","black","black","red3"),
                title="Human Oocytes Control vs 5mM",
                legendLabels = c("Not Significant","","Not Significant","Significant"),
                pCutoff=0.05,
                drawConnectors=T)