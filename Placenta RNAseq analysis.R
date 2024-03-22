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
library(biomaRt)
library(msigdbr)

#Load file
counts <- readRDS("~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Year 1/STZ project/RNA/RNAseq analysis exp2/AZ_STZ_Plac_22.dgecounts.rds")

#get all introns and exons
counts2<-counts$umicount$inex$all

#import gene names
genenames<-read.table("AZ_STZ_Plac_22.gene_names.txt", header=TRUE)

#merge files
df.plac2<-merge(counts2, genenames, by.x=0, by.y='gene_id')

#remove ensmus
df.place2 <- dplyr::select(df.plac2, -c(Row.names))

#add all transcript variants
df.place2.all<-aggregate(df.place2[,1:25], list(df.place2$gene_name), sum)

#put gene names as row names
rnames.place2<-df.place2.all$Group.1
df.place2.all$Group.1<-NULL
rownames(df.place2.all)<-rnames.place2

#make df to df1
df1<-df.place2.all

#import sample info
sampleinfo<-read.table("sample_info_2.txt",sep = "\t", header=TRUE)

#Add treatment to sampleinfo
treatment<-factor(c(rep("Control",7),rep("Diabetes",18)))
sampleinfo$treatment<-treatment
sampleinfo
sampleinfo$mother<-c('C1','C1','C2','C2','C2','C3','C3','D4','D5','D5','D5','D5','D5','D5','D6','D7','D7','D7','D7','D7','D7','D8','D8','D8','D8')

#make order
df1<-df1[,order(colnames(df1))]

#fix sampleinfo
rownames(sampleinfo)<-sampleinfo$barcode
sampleinfo$barcode<-NULL
sampleinfo <- sampleinfo[order(row.names(sampleinfo)), ]

#Create Seurat object
seurat.obj<-CreateSeuratObject(counts=df1, min.cells=3, min.features=200, project="STZ")
seurat.obj<-AddMetaData(object=seurat.obj, metadata=sampleinfo)

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
selected_c <- WhichCells(seurat.obj, expression = nCount_RNA >2000000 )
selected_f <- rownames(seurat.obj)[ Matrix::rowSums(seurat.obj) > 1]
data.filt <- subset(seurat.obj, features=selected_f, cells=selected_c)
dim(data.filt)
placenta.obj <- data.filt

#Normalize
placenta.obj <- NormalizeData(object = placenta.obj, normalization.method = "LogNormalize", scale.factor = 10000)
placenta.obj <- ScaleData(object = placenta.obj)
placenta.obj <- FindVariableFeatures(object = placenta.obj, mean.function = ExpMean, 
                                dispersion.function = LogVMR, x.low.cutoff = 0.0125, nfeatures = 1500,
                                x.high.cutoff = 5, y.cutoff = 0.5) 

#Run PCA
placenta.obj <- RunPCA(object = placenta.obj, pc.genes = placenta.obj@var.genes,npcs = 21,
                  do.print = TRUE, pcs.print = 1:21,genes.print = 3)
x<-placenta.obj@reductions$pca@cell.embeddings[,1]
y<-placenta.obj@reductions$pca@cell.embeddings[,2]

#calculate PC1 and PC2 percentage
stdev<-placenta.obj@reductions$pca@stdev
var<-stdev^2
percent<-var*100/sum(var)
percent

#Plot PCA
z<-DimPlot(object = placenta.obj, reduction = "pca", 
           group.by = 'treatment', 
           pt.size = 4, 
           combine=TRUE, 
           cols=c("red3","blue"),
           shape.by="mother")+scale_shape_manual(values=1:10)
z+xlab('PC1 31.0%')+ylab('PC2 15.4%')+geom_label_repel(aes(x,y,label=placenta.obj$mother))
z+xlab('PC1 31.0%')+ylab('PC2 15.4%')

# write table with counts
write.table(data.filt@assays$RNA@counts,"counts.txt", row.names = T, sep = "\t")
y = read.table("counts.txt", header = T, row.names = 1)
colnames(y)

# Violinplots for paper
colors = c("white", "blue3")
plots <- VlnPlot(placenta.obj, features=c("Vegfd","Lox","Srgn","Pecam1"), group.by='treatment', cols=colors, ncol=2, fill="Condition", combine=F)
for (i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot() + theme(legend.position = "none")
}
CombinePlots(plots)


#fix sampleinfo before deseq2
sampleinfo<-sampleinfo[-c(3,8,12),]

#DESEQ2
dds<-DESeqDataSetFromMatrix(as.matrix(y), colData=sampleinfo, design= ~ treatment)
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)

head(res)
dim(res)
summary(res)
par(mfrow=c(1,1))

#enhancedvolcano
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj', xlim=c(-4,4), 
                xlab=substitute(paste(bold(-Log[2] ~ "fold change"))),
                ylab = substitute(paste(bold(~-Log[10] ~ "P-value"))),
                ylim=c(0,4), FCcutoff=0.5, 
                title="STZ Placenta",
                selectLab=c("Vegfd","Lox", "Srgn","Pecam1"),
                col=c("black","black","black","red3"),
                legendLabels = c("Not Significant","","Not Significant","Significant"),
                pCutoff=0.05, pointSize=3, labSize=5, 
                drawConnectors=TRUE, widthConnectors=1)



#make results tables
write.table(res, "PlacentaResultsSTZ_2.txt", sep="\t", row.names=TRUE)

#Make df with ALL genes ordered in descending order

df_gsea<-read.table("PlacentaResultsSTZ_2.txt", sep="\t", header=T, row.names=NULL)
dim(df_gsea)
df_gsea<-df_gsea[,c(1,3)]
head(df_gsea)

df_gsea<-group_by(df_gsea, log2FoldChange)%>%
  ungroup()%>%
  arrange(desc(log2FoldChange))
colnames(df_gsea)<-c('gene_id', 'log2FoldChange')

#Convert genenames and log2FC to ENTREZiD
ensembl <- useEnsembl('ensembl', dataset = 'mmusculus_gene_ensembl')

annot_FC <- getBM(
  attributes = c(
    'external_gene_name',
    'ensembl_gene_id',
    'entrezgene_id'),
  filters = 'external_gene_name',
  values = df_gsea$gene_id,
  mart = ensembl)
head(annot_FC)
str(annot_FC)

colnames(df_gsea)<-c('gene_id', 'log2FoldChange')
head(df_gsea)
str(df_gsea)
ent_gsea<-merge(df_gsea, annot_FC[,c(1,3)], by.x=c("gene_id"), by.y=c("external_gene_name"))

ent_gsea<-ent_gsea[,c(3,2)]
head(ent_gsea)


#Once again, arrange in descending order
ent_gsea<-group_by(ent_gsea, log2FoldChange)%>%
  ungroup()%>%
  arrange(desc(log2FoldChange))
head(ent_gsea)

ent_gsea<-na.omit(ent_gsea)
head(ent_gsea)
str(ent_gsea)
tail(ent_gsea, 10)

my_list<-as.list(ent_gsea$log2FoldChange)
names(my_list)<-ent_gsea$entrezgene_id
head(my_list)
vec<-unlist(my_list)
tail(vec)

set.seed(1)

#hallmark GSEA
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

gseares<-GSEA(geneList=vec, seed=TRUE, TERM2GENE=m_t2g)
gseares@result$ID

gseaplot2(gseares, geneSetID=c(9,10), pvalue_table=T)

