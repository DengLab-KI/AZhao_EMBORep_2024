library(ArchR)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(presto)
library(hexbin)
library(ggplot2)
library(Cairo)

#SET SEED
set.seed(1)

#add threads
addArchRThreads(threads=2)

#add genome
addArchRGenome("mm10")

#create genome annotation

ArrowFiles1 <- createArrowFiles(
  inputFiles = "mouse_oocyte_frags.tsv.gz",
  sampleNames = "mouseoocyte_1",
  minTSS = 0.5,
  minFrags = 1000, 
  maxFrags = 1e+06,
  minFragSize=50,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#Arrowfiles
ArrowFiles1

#Create ArchR Project

proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles1, 
  outputDirectory = "Oocyte",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

#Explore
proj1

getAvailableMatrices(proj1)

#Fix condition in metadata
head(proj1$cellNames)

proj1$cellNames

sampleinfo_ATAC<-read.table("sampleinfo_2.txt", header=TRUE)
condition<-sampleinfo_ATAC$Condition
batch<-sampleinfo_ATAC$Batch
condition<-unlist(condition)
batch<-unlist(batch)
proj1$Condition<-condition
proj1$Batch<-batch


# TSSEnrichment ridge plot
p1 <- plotGroups(
  ArchRProj = proj1, 
  groupBy = "Condition", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

# TSSEnrichment violin plot
p2 <- plotGroups(
  ArchRProj = proj1, 
  groupBy = "Condition", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2

# unique fragments ridge plot
p3 <- plotGroups(
  ArchRProj = proj1, 
  groupBy = "Condition", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p3

# unique fragments violin plot
p4 <- plotGroups(
  ArchRProj = proj1, 
  groupBy = "Condition", 
  colorBy = "cellColData", 
  name = "nFrags",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4

#QC continued

p5 <- plotFragmentSizes(ArchRProj = proj1,
                        groupBy = "Condition")
p5

p6 <- plotTSSEnrichment(ArchRProj = proj1,
                        groupBy = "Condition")
p6

#Dimensionality reduction and clustering
proj1 <- addIterativeLSI(ArchRProj =proj1, 
                               useMatrix = "TileMatrix", 
                               name = "IterativeLSI",
                               force = TRUE,
                               LSIMethod = 2,
                               depthCol = "nFrags",
                               projectCellsPre = TRUE)

#Clustering
proj1 <- addClusters(input = proj1, reducedDims = "IterativeLSI",method = "Seurat",
                           name = "Clusters",
                           resolution = 0.5,
                           force = TRUE)

#UMAP
proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "IterativeLSI", name = "ScranClusters",force = TRUE,
                 nNeighbors = 40)

#Plot UMAPs for condition and clusters
umap <- plotEmbedding(ArchRProj =proj1, colorBy = "cellColData", name = "Condition", embedding = "ScranClusters", size=2)
umap+theme(text=element_text(size=16, family="Arial"))


#Find differences between groups
markersCond<-getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerListCond<-getMarkers(markersCond, cutOff="FDR<=0.1 & Log2FC>=0")
markerListCond$Control
markerListCond$STZ

## making pseudo-bulk ATAC replicates #by a selected cluster ##

proj_pseudo_bulk<- addGroupCoverages(ArchRProj = proj1, groupBy = "Condition" )
pathToMacs2 <- findMacs2()

#peak merging #Clusters2
pseudo_bulk_Peaks_Macs2 <- addReproduciblePeakSet(
  ArchRProj = proj_pseudo_bulk, 
  groupBy = "Condition", 
  pathToMacs2 = pathToMacs2)
getPeakSet(pseudo_bulk_Peaks_Macs2)

#Add Peak Matrix
proj_Peak_Matrix <- addPeakMatrix(pseudo_bulk_Peaks_Macs2)
getAvailableMatrices(proj_Peak_Matrix)

#Identifying Marker Peaks
table(proj_Peak_Matrix$Condition)

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_Peak_Matrix, 
  useMatrix = "PeakMatrix", 
  groupBy = "Condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

markerPeaksList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

markerPeaksList$Control
markerPeaksList$STZ

#Differential expression ######
markerTest <- getMarkerFeatures(
  ArchRProj = proj_Peak_Matrix, 
  useMatrix = "PeakMatrix",
  groupBy = "Condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "STZ",
  bgdGroups = "Control"
)
pma <- markerPlot(seMarker = markerTest, name = "STZ", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pma

pv <- markerPlot(seMarker = markerTest, name = "STZ", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv


