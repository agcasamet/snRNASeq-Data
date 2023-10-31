#Simplified verison of the analysis

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DESeq2)
library(circlize)
library(nichenetr)
library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(enrichR)
require(DOSE)
library(stringr)

#Load Data
cac.data <- Read10X(data.dir = "D:/snRNA-Seq/snRNA_01/LLC/cac_filtered_feature_bc_matrix")
control.data <- Read10X(data.dir = "D:/snRNA-Seq/snRNA_01/Control/con_filtered_feature_bc_matrix")

#Create Seurat Object
Control <- CreateSeuratObject(control.data, project = "Control")
Cachexia <- CreateSeuratObject(cac.data, project = "Cachexia")

#Remove data for emptying space
rm(control.data)
rm(cac.data)

#Calculate percent mt-RNA content
Control[["percent.mt"]] <- PercentageFeatureSet(Control, pattern = "^mt-")
Cachexia[["percent.mt"]] <- PercentageFeatureSet(Cachexia, pattern = "^mt-")

#Visualize Feature, Count, mt-RNA
VlnPlot(Control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Cachexia, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter out data
Control <- subset(Control, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 1)
Cachexia <- subset(Cachexia, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 1)

#Create combined list - prepare for integration
combined_list <- list()
combined_list[["Control"]] <- Control
combined_list[["Cachexia"]] <- Cachexia

for (i in 1:length(combined_list)) {
  combined_list[[i]] <- NormalizeData(combined_list[[i]], verbose = T)
  combined_list[[i]] <- FindVariableFeatures(combined_list[[i]], selection.method = "vst", nfeatures = 3000, verbose = T)
}

#Integration steps
combined_anchors <- FindIntegrationAnchors(object.list = combined_list, dims = 1:40)
combined_seurat <- IntegrateData(anchorset = combined_anchors, dims = 1:40)

rm(combined_list)
rm(combined_anchors)

DefaultAssay(combined_seurat) <- "RNA"
combined_seurat <- NormalizeData(combined_seurat, verbose = F)
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 3000, verbose = F)
combined_seurat <- ScaleData(combined_seurat, verbose = F)
combined_seurat <- RunPCA(combined_seurat, npcs = 40, verbose = F)
combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:40, verbose = F)
#NON-INTEGRATED
DimPlot(combined_seurat,reduction = "umap", split.by = "orig.ident", pt.size=1) + plot_annotation(title = "TA muscle of Mouse, Healty and Cachexia Inoculated - Before Integration")

DefaultAssay(combined_seurat) <- "integrated"
combined_seurat <- ScaleData(combined_seurat, verbose = F)
combined_seurat <- RunPCA(combined_seurat, npcs = 40, verbose = F)
combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:40, verbose = F)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:40, k.param = 10, verbose = F)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.8, verbose = F)
#Integrated
DimPlot(combined_seurat,label = T, label.size = 4.5, pt.size=0.7,  split.by = "orig.ident") + NoLegend()

#Find marker genes for clusters
nuclei.markers <- FindAllMarkers(combined_seurat, split.by="orig.ident", only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)

nuclei.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>%
  top_n(n=10,wt=avg_log2FC) -> top10markers

#Visualize marker genes (Dot Plot)
DefaultAssay(combined_seurat) <- "RNA"
marker_genes <- c("Ttn","Myh7","Myh2","Myh1","Myh4","Myh11","Pecam1", "Pdgfra", "Mkx","Mrc1","Pax7","Col22a1","Ache","Reln")
DotPlot(combined_seurat, features = marker_genes, cols=c("lightgrey","red3")) + RotatedAxis()

#Assign clusters
new.cluster.ids <- c("Type IIx","Type IIb","Type IIb","Type IIb","Type IIx","Type IIx/b", "Type IIb", 
                     "Type IIb", "FAPs", "Endothelial", "MuSC","Smooth M.","MTJ","FAPs","Type IIx/b",
                     "FAPs", "Type IIa", "Macrophage","Tenocyte","NMJ","Endothelial","FAPs", "Pericyte", "Smooth M.")
names(new.cluster.ids) <- levels(combined_seurat)
combined_seurat <- RenameIdents(combined_seurat, new.cluster.ids)
DimPlot(combined_seurat,label = T, label.size = 4.5, pt.size=0.7,  split.by = "orig.ident", group.by="celltype") + NoLegend()

#Find DEGs
output <- "PATH"
combined_seurat$celltype.condition <- paste(Idents(combined_seurat), combined_seurat$orig.ident, sep ="_")
combined_seurat$celltype <- Idents(combined_seurat)
Idents(combined_seurat) <- "celltype.condition"
cluster.ids <-  c("Type IIx","Type IIb","FAPs","Endothelial","MTJ","MuSC","Smooth M.","Type IIa","Macrophage","NMJ","Tenocyte","Pericyte")

for (i in cluster.ids ){
  try({
    ident1 <- paste0(i,"_Control")
    ident2 <- paste0(i,"_Cachexia")
    condition.diffgenes <- FindMarkers(combined_seurat, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25)
    write.csv(condition.diffgenes, file=paste0(output,i,".csv"))
  })
}


#FeaturePlots
Idents(combined_seurat) <- "celltype"
FeaturePlot(combined_seurat, features = marker_genes, order=TRUE, ncol=5,cols=c("lightgrey","yellow2","red3")) & NoLegend()

#ViolinPlots
VlnPlot(combined_seurat, features = marker_genes ,split.by = "orig.ident",ncol=2)

#Heatmap
DoHeatmap(combined_seurat, features = top10$gene , size=4) + scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")

#GSEA - KEGG Pathway Analysis ----

# reading in data from deseq2
df = read.csv("PATH/X.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$avg_log2FC

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order
gene_list = sort(gene_list, decreasing = TRUE)

organism = "org.Mm.eg.db" 

gse <- gseGO(geneList      = gene_list, 
              ont          = "BP", 
              keyType      = "SYMBOL", 
              nPerm        = 10000, 
              minGSSize    = 3, 
              maxGSSize    = 800, 
              pvalueCutoff = 0.05, 
              verbose      = TRUE, 
              OrgDb        = organism, )


dotplot(gse, showCategory=20)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ridgeplot(gse) + labs(x = "enrichment distribution")

# Convert gene IDs for gseKEGG function
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

# remove duplicate IDS
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$X %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$avg_log2FC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "mmu"
kk <- gseKEGG(geneList     = kegg_gene_list,
              organism     = kegg_organism,
              nPerm        = 10000,
              minGSSize    = 5,
              maxGSSize    = 800,
              pvalueCutoff = 0.05,
              keyType      = "ncbi-geneid")

dotplot(kk, showCategory = 15, title = "KEGG - Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dotplot(kk, showCategory = 30)