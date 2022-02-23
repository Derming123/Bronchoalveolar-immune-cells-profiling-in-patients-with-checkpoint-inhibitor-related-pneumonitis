setwd('D:\\cpf\\1_Cellranger_result\\matrix')
library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratWrappers)
batch.color <-c('#15966E','#696758','#DC9A15')
sample.color <- c('#A9A9C5','#CAC6B9','#B38A8E','#4D828E','#ED8415','#4A8269','#A4CD78','#A4CD78','#EF9196','#EAC326','#69BCB9','#3D599D','#BB4393',
                  '#643576','#E83F2F','#3D2F30','#305E7A')
C1.data <- Read10X('./C1/')
colnames(x = C1.data) <- paste('C1', colnames(x = C1.data), sep = '_')
C1<- CreateSeuratObject(counts  =C1.data, 
                        min.cells = 10, min.features = 300,project = 'C1') 
C2.data <- Read10X('./C2/')
colnames(x = C2.data) <- paste('C2', colnames(x = C2.data), sep = '_')
C2<- CreateSeuratObject(counts  =C2.data, 
                        min.cells = 10, min.features = 300,project = 'C2') 
C3.data <- Read10X('./C3/')
colnames(x = C3.data) <- paste('C3', colnames(x = C3.data), sep = '_')
C3<- CreateSeuratObject(counts  =C3.data,
                        min.cells = 10, min.features = 300,project = 'C3')
C4.data <- Read10X('./C4/')
colnames(x = C4.data) <- paste('C4', colnames(x = C4.data), sep = '_')
C4<- CreateSeuratObject(counts  =C4.data,
                        min.cells = 10, min.features = 300,project = 'C4')
C5.data <- Read10X('./C5/')
colnames(x = C5.data) <- paste('C5', colnames(x = C5.data), sep = '_')
C5<- CreateSeuratObject(counts  =C5.data,
                        min.cells = 10, min.features = 300,project = 'C5')
C6.data <- Read10X('./C6/')
colnames(x = C6.data) <- paste('C6', colnames(x = C6.data), sep = '_')
C6<- CreateSeuratObject(counts  =C6.data,
                        min.cells = 10, min.features = 300,project = 'C6')
P1.data <- Read10X('./P1/')
colnames(x = P1.data) <- paste('P1', colnames(x = P1.data), sep = '_')
P1<- CreateSeuratObject(counts  =P1.data,
                        min.cells = 10, min.features = 300,project = 'P1')
P2.data <- Read10X('./P2/')
colnames(x = P2.data) <- paste('P2', colnames(x = P2.data), sep = '_')
P2<- CreateSeuratObject(counts  =P2.data,
                        min.cells = 10, min.features = 300,project = 'P2')
P3.data <- Read10X('./P3/')
colnames(x = P3.data) <- paste('P3', colnames(x = P3.data), sep = '_')
P3<- CreateSeuratObject(counts  =P3.data,
                        min.cells = 10, min.features = 300,project = 'P3')
P4.data <- Read10X('./P4/')
colnames(x = P4.data) <- paste('P4', colnames(x = P4.data), sep = '_')
P4<- CreateSeuratObject(counts  =P4.data,
                        min.cells = 10, min.features = 300,project = 'P4')
P5.data <- Read10X('./P5/')
colnames(x = P5.data) <- paste('P5', colnames(x = P5.data), sep = '_')
P5<- CreateSeuratObject(counts  =P5.data,
                        min.cells = 10, min.features = 300,project = 'P5')
P6.data <- Read10X('./P6/')
colnames(x = P6.data) <- paste('P6', colnames(x = P6.data), sep = '_')
P6<- CreateSeuratObject(counts  =P6.data,
                        min.cells = 10, min.features = 300,project = 'P6')


pbmc <- merge(C1,c(C2,C3,C4,C5,C6,P1,P2,P3,P4,P5,P6))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
percent.CD <- pmax(pbmc@assays$RNA@data["CD3D", ],pbmc@assays$RNA@data["CD79A", ])/Matrix::colSums(pbmc@assays$RNA@data[c("CD3D","CD79A"), ])
percent.CD[is.na(percent.CD)] <- 1
percent.CD[(pbmc@assays$RNA@data["CD3D", ] < 2) & (pbmc@assays$RNA@data["CD79A", ] < 2)] <- 1
pbmc <- AddMetaData(object = pbmc, metadata = percent.CD, col.name = "percent.CD")

percent.BD <- pmax(pbmc@assays$RNA@data["CD3D", ],pbmc@assays$RNA@data["CD68", ])/Matrix::colSums(pbmc@assays$RNA@data[c("CD3D","CD68"), ])
percent.BD[is.na(percent.BD)] <- 1
percent.BD[(pbmc@assays$RNA@data["CD3D", ] < 2) & (pbmc@assays$RNA@data["CD68", ] < 2)] <- 1
pbmc <- AddMetaData(object = pbmc, metadata = percent.BD, col.name = "percent.BD")

percent.DD <- pmax(pbmc@assays$RNA@data["EMCN", ],pbmc@assays$RNA@data["KRT18", ])/Matrix::colSums(pbmc@assays$RNA@data[c("EMCN","KRT18"), ])
percent.DD[is.na(percent.DD)] <- 1
percent.DD[(pbmc@assays$RNA@data["EMCN", ] < 2) & (pbmc@assays$RNA@data["KRT18", ] < 2)] <- 1
pbmc <- AddMetaData(object = pbmc, metadata = percent.DD, col.name = "percent.DD")

percent.ED <- pmax(pbmc@assays$RNA@data["CD68", ],pbmc@assays$RNA@data["KRT18", ])/Matrix::colSums(pbmc@assays$RNA@data[c("CD68","KRT18"), ])
percent.ED[is.na(percent.ED)] <- 1
percent.ED[(pbmc@assays$RNA@data["CD68", ] < 2) & (pbmc@assays$RNA@data["KRT18", ] < 2)] <- 1
pbmc <- AddMetaData(object = pbmc, metadata = percent.ED, col.name = "percent.ED")

percent.FD <- pmax(pbmc@assays$RNA@data["CD68", ],pbmc@assays$RNA@data["EMCN", ])/Matrix::colSums(pbmc@assays$RNA@data[c("CD68","EMCN"), ])
percent.FD[is.na(percent.FD)] <- 1
percent.FD[(pbmc@assays$RNA@data["CD68", ] < 2) & (pbmc@assays$RNA@data["EMCN", ] < 2)] <- 1


percent.AD <- pmax(pbmc@assays$RNA@data["PTPRC", ],pbmc@assays$RNA@data["KRT19", ])/Matrix::colSums(pbmc@assays$RNA@data[c("PTPRC","KRT19"), ])
percent.AD[is.na(percent.AD)] <- 1
percent.AD[(pbmc@assays$RNA@data["PTPRC", ] < 2) & (pbmc@assays$RNA@data["KRT19", ] < 2)] <- 1
pbmc <- AddMetaData(object = pbmc, metadata = percent.AD, col.name = "percent.AD")

percent.GD <- pmax(pbmc@assays$RNA@data["PTPRC", ],pbmc@assays$RNA@data["ALB", ])/Matrix::colSums(pbmc@assays$RNA@data[c("PTPRC","ALB"), ])
percent.GD[is.na(percent.GD)] <- 1
percent.GD[(pbmc@assays$RNA@data["PTPRC", ] < 2) & (pbmc@assays$RNA@data["ALB", ] < 2)] <- 1
pbmc <- AddMetaData(object = pbmc, metadata = percent.GD, col.name = "percent.GD")

pbmc <- AddMetaData(object = pbmc, metadata = percent.FD, col.name = "percent.FD")
VlnPlot(pbmc,features = c('percent.CD','percent.DD','percent.BD','percent.ED','percent.FD','percent.mt','percent.GD','percent.AD'),ncol = 4,pt.size = 0)
FeatureScatter(pbmc,'nFeature_RNA','nCount_RNA')
pbmc <-  subset(pbmc,percent.ED > 0.9 & percent.CD > 0.90 & percent.BD > 0.9 & percent.DD > 0.9 & percent.FD > 0.9 & percent.AD > 0.9 & percent.GD > 0.9 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 &nCount_RNA < 40000)

### second
library(DoubletFinder)
library(dplyr)
library(stringr)
FindDoublets <- function(library_id, seurat_aggregate) {
  rnaAggr <- seurat_aggregate
  seurat_obj <- subset(rnaAggr, idents = library_id)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj)
  # ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  DimPlot(seurat_obj)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_kidney <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = F)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK,
                                      nExp = round(0.05*length(seurat_obj@active.ident)), 
                                      reuse.pANN = FALSE, sct = F)
  
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  
  #p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +geom_bar(stat = "identity") +  ggtitle(paste0("pKmax=",pK)) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  #p3 <- FeaturePlot(seurat_doublets, features = pANN)
  
  # create a df of barcodes and doublet designations
  df_doublet_barcodes <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  return(df_doublet_barcodes)
}


# take an aggregated snRNA seurat object and a list library_id and find doublets. return a df of doublet barcodes
# send DimPlot and FeaturePlot of doublets for each library to here("plots")
Idents(pbmc) <- "orig.ident"
rnaAggr <- pbmc
list.doublet.bc <- lapply(levels(rnaAggr), function(x) {FindDoublets(x, seurat_aggregate = rnaAggr)})
doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1") # this is the barcode column
table(doublet_id) # quantify total doublet vs. singlet calls (expect ~5% doublets)

# add doublet calls to aggregated snRNA object as doublet_id in meta.data slot
rnaAggr <- AddMetaData(rnaAggr,doublet_id)

# filter out doublets prior to snRNA preprocessing
Idents(rnaAggr) <- "doublet_id"
saveRDS(rnaAggr,'./doubletid.lung2.rds')

###
scRNA  <- subset(doublet_id,idents = "Singlet")
scRNAlist <- SplitObject(pbmc2,split.by = 'orig.ident')
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))
scRNA <- RunFastMNN(object.list = scRNAlist)
scRNA <- RunUMAP(scRNA, reduction = "mnn", dims = 1:30)
scRNA <- FindNeighbors(scRNA, reduction = "mnn", dims = 1:30)
scRNA <- FindClusters(scRNA,resolution = 0.8)
DimPlot(scRNA, pt.size=0.1,label = T)  
pbmc <- scRNA
pbmc$batch <-substr(pbmc$orig.ident, 1,1)
library(harmony)
pbmc  <- subset(doubletid.lung2,idents = "Singlet")
percent.BD <- pmax(pbmc@assays$RNA@data["CD3D", ],pbmc@assays$RNA@data["CD68", ])/Matrix::colSums(pbmc@assays$RNA@data[c("CD3D","CD68"), ])
percent.BD[is.na(percent.BD)] <- 1
percent.BD[(pbmc@assays$RNA@data["CD3D", ] < 2) & (pbmc@assays$RNA@data["CD68", ] < 2)] <- 1
pbmc <- AddMetaData(object = pbmc, metadata = percent.BD, col.name = "percent.BD")
pbmc2 <-  subset(pbmc,  percent.BD > 0.99 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 &nCount_RNA < 40000)
pbmc <- CreateSeuratObject(pbmc2@assays$RNA@counts[rowSums(as.matrix(pbmc2@assays$RNA@counts) > 0) > 10,],meta.data = pbmc2@meta.data)
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <-FindVariableFeatures(pbmc)
pbmc <- RunPCA(pbmc,features = VariableFeatures(pbmc),npcs = 30,verbose = F)

pbmc <- pbmc %>% 
  RunHarmony('orig.ident',plot_convergence=TRUE)
harmony_embeddings <- Embeddings(pbmc,'harmony')
ElbowPlot(pbmc)
pbmc <- pbmc %>%
  RunTSNE(reduction='harmony',dims=1:30)%>%
  RunUMAP(reduction='harmony',dims=1:30)%>%
  FindNeighbors(reduction='harmony',dims=1:30) %>%
  FindClusters(resolution=0.8) %>%
  identity()

pbmc <- pbmc %>% 
  FindNeighbors(reduction='harmony',dims=1:30) %>%
  FindClusters(resolution=0.7) %>%
  identity()


DimPlot(pbmc,reduction = 'umap',label = T)
pbmc3.markers <- FindAllMarkers(pbmc3,logfc.threshold = log(1.5),test.use = 'wilcox')


library(dplyr)
top10 <- pbmc3.markers %>% filter(avg_log2FC>0 & p_val<0.05) %>% group_by(cluster) %>% top_n(wt = avg_log2FC,5)
pbmc3$batch <- plyr::mapvalues(x = pbmc3$orig.ident,from  = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                               to = c('CIP(-)', 'CIP(-)', 'CIP(-)', 'Healthy', 'Healthy', 'Healthy', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)' ))
pbmc3$batch <- factor( pbmc3$batch,levels = c( 'Healthy',  'CIP(-)',  'CIP(+)' ))

pbmc3$sample <- plyr::mapvalues(x = pbmc3$orig.ident,from  = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                                to = c('C1', 'C2', 'C3', 'H1', 'H2', 'H3', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ))
pbmc3$sample <- factor(x = pbmc3$sample,levels  = c('H1', 'H2', 'H3' ,'C1', 'C2', 'C3', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ))
pbmc$cluster <- plyr::mapvalues(pbmc$seurat_clusters,from = c('22','6','19','18','14','7','24','25','1','8','4','23','12','21','16','9','11','20','5','3','13','10','0','2','15','17'),to = c('1','1','2','3','4','5','5','5','6','6','7','7','8','9','10','11','12','13','14','15','15','16','17','17','17','17'))
pbmc$cluster <- factor(pbmc$cluster,levels  = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'))
Idents(pbmc) <- 'cluster'
DimPlot(pbmc3,reduction = 'umap',label = T,group.by = 'cluster')
pbmc3 <- pbmc [,sample(row.names(pbmc@meta.data),30000)]

p1 <- DimPlot(pbmc3,cols = sample.color,label = T)
save(file = './pbmc.Rds',pbmc,pbmc3,pbmc3.markers,top10)

library(reshape2)
plot_gene_matrix <- function(data, gene, anno, col){
  
  gene_data<- data[gene,]
  
  violin_data <- data.frame(cluster=anno[,'celltype'],t(gene_data[,rownames(anno)]),check.names = F)
  violin_data <- melt(violin_data, id='cluster')
  names(violin_data) <- c('cluster','gene','value')
  
  ggplot(data=violin_data,aes(x=cluster,y=value,fill=cluster))+geom_violin()+coord_flip()+facet_grid(cluster~gene,scales = 'free')+
    scale_fill_manual(values = col)+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_blank())+
    theme(strip.text.x = element_text(colour = "black", size=rel(0.8)),strip.background.y =element_rect(fill=c('red','green')),strip.text.y = element_blank())+labs(y='log2 normalized value')
}


plot_gene_matrix(pbmc3@assays$RNA@data,gene = c('SCGB3A1','CAPS','LCN2','CD3D','CD8A','CD4','FOXP3','KLRF1','XCL2','JCHAIN','CD79A','MKI67','CLEC9A','CD1C','LAMP3','S100A8','APOE','C1QA','MARCO'),anno = pbmc3@meta.data,col = rev(sample.color))

VlnPlot(pbmc3 ,features =  c('SCGB3A1','CAPS','LCN2','CD3D','CD8A','CD4','FOXP3','KLRF1','XCL2','JCHAIN','CD79A','MKI67','CLEC9A','CD1C','LAMP3','S100A8','APOE','C1QA','MARCO'),cols =   rev(sample.color),pt.size = 0)

pbmc$cluster <- plyr::mapvalues(pbmc$seurat_clusters,from = c('22','6','19','18','14','7','24','25','1','8','4','23','12','21','16','9','11','20','5','3','13','10','0','2','15','17'),
                                to = c('1','1','2','3','4','5','5','5','6','18','7','7','8','9','10','11','12','13','14','15','15','16','17','17','17','17'))

pbmc$cluster <- plyr::mapvalues(pbmc$seurat_clusters,from =c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'),to = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'))
pbmc <- pbmc [,!pbmc$cluster %in%'5']
pbmc$celltype <- plyr::mapvalues(pbmc$cluster,from =c('1','2','3','4', '6','7','8','9','10','11','12','13','14','15','16','17','18'),to = c('Ciliated cells 1','Ciliated cells 2','mDC1-CLEC9A','Secretory cells', 'CD8','CD4-Naive','Treg','Plasma','B cell','Cycling','mDC2-CD1C','pDC-LAMP3','Mono/Mac1','Mono/Mac2','Mono/Mac3','Mono/Mac4','NK'))
pbmc$celltype <- factor(pbmc$celltype,levels = c('Ciliated cells 1','Ciliated cells 2','Secretory cells', 'CD8','CD4-Naive','Treg','NK','Plasma','B cell','Cycling','mDC1-CLEC9A','mDC2-CD1C','pDC-LAMP3','Mono/Mac1','Mono/Mac2','Mono/Mac3','Mono/Mac4'))
Idents(pbmc) <- 'celltype'

genes <- read.table(file = 'D:/code/ligand_receptor.txt',header = T)
cytokines <- read.csv(file ='D:\\code\\common_data\\CD分子查询\\cytokine.csv',header = T)


Idents(myeloid2) <- 'batch'
myeloid2 <- FindVariableFeatures(myeloid2)
markers <- mye.markers
markers <- FindAllMarkers(myeloid2,features = VariableFeatures(myeloid2),test.use = 'wilcox',logfc.threshold = log(1.5))
myeloid2$batch <- factor(myeloid2$batch,levels = c('Healthy','CIP(-)','CIP(+)'))
markers.sig_Ligand <- subset(markers, markers$gene %in% genes$Ligand)
markers.sig_Ligand <- subset(markers.sig_Ligand, markers.sig_Ligand$gene %in% cytokines$gene)
markers.sig_Ligand %>% group_by(cluster) %>% filter(avg_log2FC >0)  %>% top_n(40, -p_val_adj) -> top10_Ligand

markers.sig_receptor <- subset(markers, markers$gene %in% genes$Receptor)
markers.sig_receptor <- subset(markers.sig_receptor, markers.sig_receptor$gene %in% cytokines$gene)
markers.sig_receptor %>% group_by(cluster) %>% filter(avg_log2FC >0)  %>% top_n(10, -p_val_adj) -> top15_receptor


gene_data <- data.frame(t(as.matrix(myeloid2@assays$RNA@data)[c('IL10RB','IL4R','LTBR','TNFRSF1A','IFNGR1','IFNGR2','IL10RA','IL1R2','TGFBR1'),]),cluster=myeloid2@meta.data$batch,check.names = F)
average_data <- aggregate(.~cluster, gene_data, mean)
cluster_name <- average_data[,1]
average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
rownames(average_data) <- cluster_name
average_data <- t(average_data)
phmat1 <- t(scale(t(average_data)))
phmat1[phmat1> 2] <- 2
phmat1[phmat1 < -2] <- -2
#degs
library(pheatmap)
pheatmap(phmat1,cluster_rows = F,cluster_cols = F,fontsize = 8,colorRampPalette(colors = c('#789FBA','white','#882F21'))(100))



Idents(Tcell2) <- 'batch'
Tcell2$batch <- factor(Tcell2$batch,levels = c('Healthy','CIP(-)','CIP(+)'))
Tcell2 <- FindVariableFeatures(Tcell2)
markers <- FindAllMarkers(Tcell2,features = VariableFeatures(Tcell2),test.use = 'wilcox')

markers.sig_Ligand <- subset(markers, markers$gene %in% genes$Ligand)
markers.sig_Ligand <- subset(markers.sig_Ligand, markers.sig_Ligand$gene %in% cytokines$gene)
markers.sig_Ligand %>% group_by(cluster) %>% filter(avg_log2FC >0)  %>% top_n(10, -p_val_adj) -> top10_Ligand

markers.sig_receptor <- subset(markers, markers$gene %in% genes$Receptor)
markers.sig_receptor <- subset(markers.sig_receptor, markers.sig_receptor$gene %in% cytokines$gene)
markers.sig_receptor %>% group_by(cluster) %>% filter(avg_log2FC >0)  %>% top_n(10, -p_val_adj) -> top10_receptor


gene_data <- data.frame(t(as.matrix(myeloid2@assays$RNA@data)[c(),]),cluster=myeloid2@meta.data$batch,check.names = F)
average_data <- aggregate(.~cluster, gene_data, mean)
cluster_name <- average_data[,1]
average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
rownames(average_data) <- cluster_name
average_data <- t(average_data)
phmat1 <- t(scale(t(average_data)))
phmat1[phmat1> 2] <- 2
phmat1[phmat1 < -2] <- -2
#degs
library(pheatmap)
pheatmap(phmat1,cluster_rows = F,cluster_cols = F,fontsize = 8,colorRampPalette(colors = c('#789FBA','white','#882F21'))(100))

annot <- CD8@meta.data
cell_num <- table(as.character(annot$orig.ident))
num_table <- table(annot$orig.ident, annot$celltype)[names(cell_num),]
num_table <- data.frame(num_table/as.numeric(cell_num))
names(num_table) <- c('sample','cluster','proportion')
library(stringr)
num_table$batch <- plyr::mapvalues(x = num_table$sample,from  = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                                   to = c('CIP(-)', 'CIP(-)', 'CIP(-)', 'Healthy', 'Healthy', 'Healthy', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)' ))
#加载包
library(ggpubr)
library(reshape)
compar <- list(c("CIP(+)", "Healthy" ),c("CIP(-)", "Healthy" ), c('CIP(+)','CIP(-)'))
num_table$batch <- factor(num_table$batch,levels = c('Healthy', 'CIP(-)','CIP(+)'))
ggboxplot(num_table, x = "batch", y = "proportion", color = "batch",palette = "jco", add = "jitter", facet.by = "cluster", short.panel.labs = FALSE)+ stat_compare_means(label="p.format",method = 't.test',comparisons = compar)+theme(legend.position="right",axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(ncol =8,'cluster',scales = 'free')+labs(y='Cell Proportion (%)')
ggboxplot(num_table, x = "cluster", y = "proportion", color = "batch",  add = "jitter",palette = batch.color) +theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(x="", y='Cell Proportion (%)') +theme(legend.position="right")+ stat_compare_means(label="p.format",method = 't.test',comparisons = compar)
ggboxplot(num_table, x = "cluster", y = "proportion", color = "batch",  add = "jitter",palette = batch.color) +theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(x="", y='Cell Proportion (%)') +theme(legend.position="right")+ stat_compare_means(label="p.format",method = 't.test',comparisons = compar)

ggbarplot(num_table, x = "cluster", y = "proportion", add = "mean_se",fill = 'batch',
          color = "batch", palette = batch.color, 
          position = position_dodge(0.8)) + theme_bw()+ stat_compare_means(comparisons = compar, method = 't.test',label = "p.format" ) +theme(legend.position="right",axis.text.x = element_text(angle = 45, hjust = 1))

#### pie plot
pie_data <- as.data.frame(table(pbmc@meta.data$celltype))
labels = paste(pie_data$Var1,':',round(pie_data$Freq/sum(pie_data$Freq)*100),'%',"(", pie_data$Freq, ")",sep = '')
pie(x = pie_data$Freq, labels = labels, col = my_cor2)


ph_colors_sub <- list(cluster= setNames(pal_dolphin[1:length(unique(CST_sub@meta.data$cluster))], sort(unique(CST_sub@meta.data$cluster))),
                      subnames = setNames(pal_dolphin[1:length(unique(CST_sub@meta.data$subnames))], sort(unique(CST_sub@meta.data$subnames))))



### myeloid
myeloid <- pbmc  [, pbmc$celltype %in% c('mDC1-CLEC9A','mDC2-CD1C','pDC-LAMP3','Mono/Mac1','Mono/Mac2','Mono/Mac3','Mono/Mac4','Cycling')]
myeloid <- myeloid [, !myeloid$seurat_clusters %in% c('2','4','7','8','9','11','14')]
myeloid <- CreateSeuratObject(myeloid@assays$RNA@counts[rowSums(as.matrix(myeloid@assays$RNA@counts) > 0) > 30,],meta.data = myeloid@meta.data)
myeloid <- NormalizeData(myeloid)
myeloid <- CellCycleScoring(myeloid,s.features =G1S,g2m.features = G2M )
myeloid <- ScaleData(myeloid,vars.to.regress  = c('G2M.Score','S.Score','nCount_RNA'))
myeloid <-FindVariableFeatures(myeloid)
myeloid <- RunPCA(myeloid,features = VariableFeatures(myeloid),npcs = 30,verbose = F)
library(harmony)
myeloid <- myeloid %>% 
  RunHarmony('orig.ident',plot_convergence=TRUE)
harmony_embeddings <- Embeddings(myeloid,'harmony')
ElbowPlot(myeloid)
myeloid <- myeloid %>%
  RunTSNE(reduction='harmony',dims=1:30)%>%
  RunUMAP(reduction='harmony',dims=1:30)%>%
  FindNeighbors(reduction='harmony',dims=1:30) %>%
  FindClusters(resolution=0.4) %>%
  identity()

myeloid <- myeloid %>% 
  FindNeighbors(reduction='harmony',dims=1:30) %>%
  FindClusters(resolution=0.2) %>%
  identity()


DimPlot(myeloid,reduction = 'umap',label = T)
myeloid.markers <- FindAllMarkers(myeloid,logfc.threshold = log(1.5),test.use = 'wilcox')
library(dplyr)
top10 <- myeloid.markers %>% filter(avg_log2FC>0 & p_val<0.05) %>% group_by(cluster) %>% top_n(wt = avg_log2FC,6)
myeloid$batch <- plyr::mapvalues(x = myeloid$orig.ident,from  = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                                 to = c('CIP(-)', 'CIP(-)', 'CIP(-)', 'Healthy', 'Healthy', 'Healthy', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)' ))


myeloid2$cluster <- plyr::mapvalues(myeloid2$seurat_clusters,from = c('9','7','10','3','6','4','5','1','0','11','8','2'),to = c('1','2','3','4','5','6','7','7','7','7','8','9'))
myeloid2$cluster <- factor(myeloid2$cluster,levels = c('1','2','3','4','5','6','7','8','9'))
myeloid2$celltype <- plyr::mapvalues(myeloid2$cluster,from = c('1','2','3','4','5','6','7','8','9'),to = c('DC-CLEC9A','DC-CD1C','DC-LAMP3','Mono1','Mono2','Mono/Nac','Mac1','Mac-Cycling','Mac2'))
myeloid2$celltype <- factor(myeloid2$celltype,levels = c('DC-CLEC9A','DC-CD1C','DC-LAMP3','Mono1','Mono2','Mono/Nac','Mac1','Mac-Cycling','Mac2'))
Idents(myeloid2) <- 'celltype'
DimPlot(myeloid2, pt.size=0.1,label = T ) 
myeloid2 <- FindClusters(myeloid2,resolution = 0.6)
myeloid2$celltype <- plyr::mapvalues(myeloid2$seurat_clusters,from = c('12','9','13','6','5','8','2','11','0','4','3','7','1','10'),to = c('DC-CLEC9A','DC-CD1C','DC-LAMP3','Mono1','Mono2','Mono/Mac1','Mono/Mac2','Mono/Mac2','Mac1','Mac1','Mac1','Mac1','Mac2','Mac-Cycling'))
myeloid2$celltype <- factor(myeloid2$celltype,levels =  c('DC-CLEC9A','DC-CD1C','DC-LAMP3','Mono1','Mono2','Mono/Mac1','Mono/Mac2','Mac1','Mac2','Mac-Cycling'))

mac.sig <- read.delim2(file = 'D:/code/mac_sig.txt')
myeloid2 <- AddModuleScore(myeloid2,features = list(c('CCL5', 'CCR7', 'CD40', 'CD86', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'IL1A',
                                                      'IL1B', 'IL6', 'IRF1', 'IRF5',  'KYNU','S100A9','S100A9','S100A12'),
                                                    c('CCL4', 'CCL13', 'CCL18', 'CCL20', 'CCL22', 'CD276', 'CLEC7A', 'CTSA', 'CTSB', 'CTSC', 'CTSD','MRC1','APOE','APOC1','C1QA','C1QB','C1QC',
                                                      'FN1', 'IL4R', 'IRF4', 'LYVE1', 'MMP9', 'MMP14', 'MMP19', 'MSR1', 'TGFB1', 'TGFB2', 'TGFB3', 'TNFSF8', 'TNFSF12', 'VEGFA', 'VEGFB', 'VEGFC')),
                           name = c('M1.score','M2.score'))
geneset <- readxl::read_excel('D:/code/common_data/mmc2.xlsx')
geneset2 <- readxl::read_excel('D:/code/common_data/mmc5.xlsx')
geneset  <- na.omit(geneset )
myeloid2 <- AddModuleScore(myeloid2,features = list( geneset$`Activated DC`, geneset$`Migratory DC`,geneset2$M1,geneset2$M2,geneset2$Angiogenesis,geneset2$Phagocytosis ),name = c('Activated.score','Migratory.score','M1.score','M2.score','Angiogenesis.score','Phagocytosis.score'))
source('D:/code/function.R')
VlnPlot(myeloid2,features = c('Activated.score1','Migratory.score2','M1.score1','M2.score2'),cols = my36colors,ncol = 2,pt.size = 0)
anno <- myeloid2@meta.data[,c("Activated.score1" ,   "Migratory.score2" ,   "M1.score3"    ,       "M2.score4"    ,       'celltype')]
gsva_aggr <- aggregate(.~celltype,anno,mean)
rownames(gsva_aggr) <- gsva_aggr$celltype
gsva_aggr <- gsva_aggr[,-1]
phmat <- t(scale(t(gsva_aggr)))
phmat <- gsva_aggr
phmat [phmat>0.2] <- 0.2
phmat [phmat < -0.2] <- -0.2
pheatmap::pheatmap(t(phmat),color = colorRampPalette(c("blue","white",  "red" ))(100),show_colnames = T)



#加载包
library(ggpubr)
annot <- myeloid2@meta.data
cell_num <- table(as.character(annot$orig.ident))
num_table <- table(annot$orig.ident, annot$celltype)[names(cell_num),]
num_table <- data.frame(num_table/as.numeric(cell_num))
names(num_table) <- c('sample','cluster','proportion')
library(stringr)
Tcell2$batch <- plyr::mapvalues(x = Tcell2$orig.ident,from  = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                                to = c('CIP(-)', 'CIP(-)', 'CIP(-)', 'Healthy', 'Healthy', 'Healthy', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)' ))
#加载包
library(ggpubr)
library(reshape)
compar <- list(c("CIP(+)", "Healthy" ),c("CIP(-)", "Healthy" ), c('CIP(+)','CIP(-)'))
num_table$batch <- factor(num_table$batch,levels = c('Healthy', 'CIP(-)','CIP(+)'))
ggboxplot(num_table, x = "batch", y = "proportion", color = "batch",palette = "jco", add = "jitter", facet.by = "cluster", short.panel.labs = FALSE)+ stat_compare_means(label="p.format",method = 't.test',comparisons = compar)+theme(legend.position="right",axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(ncol =7,'cluster',scales = 'free')+labs(y='Cell Proportion (%)')
ggboxplot(num_table, x = "cluster", y = "proportion", color = "batch",  add = "jitter",palette = batch.color) +theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(x="", y='Cell Proportion (%)') +theme(legend.position="right")+ stat_compare_means(label="p.format",method = 't.test',comparisons = compar)

gene_data <- data.frame(t(as.matrix(myeloid2@assays$RNA@data)[myeloid2.top10$gene,]),cluster=myeloid2@meta.data$celltype,check.names = F)
average_data <- aggregate(.~cluster, gene_data, mean)
cluster_name <- average_data[,1]
average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
rownames(average_data) <- cluster_name
average_data <- t(average_data)
phmat1 <- t(scale(t(average_data)))
phmat1[phmat1> 2] <- 2
phmat1[phmat1 < -2] <- -2
#degs
library(pheatmap)
pheatmap(t(phmat1),cluster_rows = F,cluster_cols = F,fontsize = 8,colorRampPalette(colors = c('#789FBA','white','#882F21'))(100))




allTF_human <- read.csv("D:/code/TF_data/allTF.human_v2.txt", header = F)
allmarkers_HEC_sig_TF <- subset(CD4.markers, gene %in% allTF_human$V1)
allmarkers_HEC_sig_TF %>% group_by(cluster)%>% filter(avg_log2FC > 0) %>% top_n(10, -p_val_adj) -> top10_TF
pheatmap(CST_HEC@data[top10_TF$gene, order(CST_HEC@meta.data[, "HEC_cluster", drop = F])],
         annotation_colors = ph_colors_HEC, fontsize_row = 12, cluster_rows = F, 
         cluster_cols = F, annotation_col = CST_HEC@meta.data[,"HEC_cluster", drop = F],
         show_colnames = F,  color = colorRampPalette(colors = c("darkblue","white","red"))(100),
         border_color = NA, annotation_legend = T, main = "TF")

surface_gene <- read.table(file = "D:/code/TF_data/uniprot_reviewed_surface_marker.txt")
allmarkers_HEC_sig_sur <- subset(CD4.markers, gene %in% surface_gene$x)
allmarkers_HEC_sig_sur %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% top_n(10, -p_val_adj) -> top10_sur
pheatmap(CST_HEC@data[top10_sur$gene, order(CST_HEC@meta.data[, "HEC_cluster", drop = F])],
         annotation_colors = ph_colors_HEC, fontsize_row = 12, cluster_rows = F, 
         cluster_cols = F, annotation_col = CST_HEC@meta.data[,"HEC_cluster", drop = F],
         show_colnames = F,  color = colorRampPalette(colors = c("darkblue","white","red"))(100),
         border_color = NA, annotation_legend = T, main = "surface")



Tcell2$celltype <- plyr::mapvalues(Tcell2$seurat_clusters,from = c('0','1','2','10','8','12','9','11','6','3','7','4','5'),to = c('CD8-Cyto','CD8-Trm-ZNF683','CD8-Exha','CD8-IFIT','CD8-Cycling-G2M','CD8-Cycling-S','Mait','NKT','NK','CD4-CD40LG','CD4-Exha','CD4-Naive-IL7R','Treg'))
Tcell2$celltype <-factor(Tcell2$celltype,levels = c('CD8-Cyto','CD8-Trm-ZNF683','CD8-Exha','CD8-IFIT','CD8-Cycling-G2M','CD8-Cycling-S','Mait','NKT','NK','CD4-CD40LG','CD4-Exha','CD4-Naive-IL7R','Treg'))
Idents(Tcell2) <- 'celltype'
Tcell2$batch <- plyr::mapvalues(x = Tcell2$orig.ident,from  = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                                to = c('CIP(-)', 'CIP(-)', 'CIP(-)', 'Healthy', 'Healthy', 'Healthy', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)' ))


### CD4
CD4 <- Tcell2 [, Tcell2$celltype %in% c('CD4-CD40LG','CD4-Exha','CD4-Naive-IL7R','Treg')]
CD4 <- CD4 [,! CD4$seurat_clusters %in% c('5')]

CD4 <- CreateSeuratObject(CD4@assays$RNA@counts[rowSums(as.matrix(CD4@assays$RNA@counts) > 0) > 20,],meta.data = CD4@meta.data)
CD4 <- NormalizeData(CD4)
G1S <- toupper(c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"))
G2M <- toupper(c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3" ,"Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3", "Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"))
CD4 <- CellCycleScoring(CD4,s.features =G1S,g2m.features = G2M )
CD4 <- ScaleData(CD4)
CD4 <-FindVariableFeatures(CD4)
CD4 <- RunPCA(CD4,features = VariableFeatures(CD4),npcs = 30,verbose = F)
library(harmony)
CD4 <- CD4 %>% 
  RunHarmony('orig.ident',plot_convergence=TRUE)
harmony_embeddings <- Embeddings(CD4,'harmony')
ElbowPlot(CD4)
CD4 <- CD4 %>%
  RunTSNE(reduction='harmony',dims=1:20)%>%
  RunUMAP(reduction='harmony',dims=1:20)%>%
  FindNeighbors(reduction='harmony',dims=1:20) %>%
  FindClusters(resolution=0.4) %>%
  identity()

CD4 <- CD4 %>%  
  FindNeighbors(reduction='harmony',dims=1:30) %>%
  FindClusters(resolution=1.1) %>%
  identity()


DimPlot(CD4,reduction = 'umap',label = T)
CD4$celltype <- plyr::mapvalues(CD4$seurat_clusters,from = c('6','3','5','0','4','1','2'),to = c('CD4-Cycling','CD4-Exha','CD4-preExha','CD4-CD40LG','CD4-IFIT','CD4-Naive-IL7R', 'Treg'))
CD4$celltype <- factor(CD4$celltype,levels  = c('CD4-Cycling','CD4-Exha','CD4-preExha','CD4-CD40LG','CD4-IFIT','CD4-Naive-IL7R', 'Treg'))
Idents(CD4) <- 'celltype'
CD4.markers <- FindAllMarkers(CD4,logfc.threshold = log(1.5),test.use = 'wilcox')
library(dplyr)
CD4top10 <- CD4.markers %>% filter(avg_log2FC>0 & p_val<0.05) %>% group_by(cluster) %>% top_n(wt = avg_log2FC,10)

allTF_human <- read.csv("D:/code/TF_data/allTF.human_v2.txt", header = F)
surface_gene <- read.table(file = "D:/code/TF_data/uniprot_reviewed_surface_marker.txt")
allmarkers_HEC_sig_TF <- subset(CD4.markers, gene %in% allTF_human$V1)
allmarkers_HEC_sig_sur <- subset(CD4.markers, gene %in% surface_gene$x)
allmarkers_HEC_sig_TF %>% group_by(cluster)%>% filter(avg_log2FC > 0) %>% top_n(10, -p_val_adj) -> top10_TF
allmarkers_HEC_sig_sur %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% top_n(10, -p_val_adj) -> top10_sur

gene_data <- data.frame(t(as.matrix(CD4@assays$RNA@data)[CD4top10$gene,]),cluster=CD4@meta.data$celltype,check.names = F)
average_data <- aggregate(.~cluster, gene_data, mean)
cluster_name <- average_data[,1]
average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
rownames(average_data) <- cluster_name
average_data <- t(average_data)
phmat1 <- t(scale(t(average_data)))
phmat1[phmat1> 2] <- 2
phmat1[phmat1 < -2] <- -2
#degs
library(pheatmap)
load('D:/lab_cor.Rdata')
pheatmap(phmat1,cluster_rows = F,cluster_cols = F,fontsize = 8,colorRampPalette(colors = c('#789FBA','white','#882F21'))(100))

CD4$sample <- plyr::mapvalues(x = CD4$orig.ident,from  = c('C1', 'C2', 'C3', 'C4',    'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                              to = c('C1', 'C2', 'C3', 'H1',    'H3', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ))
CD4$sample <- factor(x = CD4$sample,levels  = c('H1',     'H3' ,'C1', 'C2', 'C3', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ))
source('D:/code/function.R')
barplot_interaction2(CD4,x = 'sample',y = 'celltype',color = my_cor1)



CD4$batch <- plyr::mapvalues(x = CD4$orig.ident,from  = c('C1', 'C2', 'C3', 'C4',   'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                             to = c('CIP(-)', 'CIP(-)', 'CIP(-)', 'Healthy',   'Healthy', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)' ))
#加载包
#加载包
library(ggpubr)
annot <- CD8@meta.data
cell_num <- table(as.character(annot$orig.ident))
num_table <- table(annot$orig.ident, annot$celltype)[names(cell_num),]
num_table <- data.frame(num_table/as.numeric(cell_num))
names(num_table) <- c('sample','cluster','proportion')
library(stringr)
num_table$batch <- plyr::mapvalues(x = num_table$sample,from  = c('C1', 'C2', 'C3', 'C4',   'C5',    'C6', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' ),
                                   to = c('CIP(-)', 'CIP(-)', 'CIP(-)', 'Healthy', 'Healthy',   'Healthy', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)', 'CIP(+)' ))
#加载包
library(ggpubr)
library(reshape)
compar <- list(c("CIP(+)", "Healthy" ),c("CIP(-)", "Healthy" ), c('CIP(+)','CIP(-)'))
num_table$batch <- factor(num_table$batch,levels = c('Healthy', 'CIP(-)','CIP(+)'))
ggboxplot(num_table, x = "batch", y = "proportion", color = "batch",palette = "jco", add = "jitter", facet.by = "cluster", short.panel.labs = FALSE)+ stat_compare_means(label="p.format",method = 't.test',comparisons = compar)+theme(legend.position="right",axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(ncol =4,'cluster',scales = 'free')+labs(y='Cell Proportion (%)')
ggboxplot(num_table, x = "cluster", y = "proportion", color = "batch",  add = "jitter",palette = batch.color) +theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(x="", y='Cell Proportion (%)') +theme(legend.position="right")+ stat_compare_means(label="p.format",method = 't.test',comparisons = compar)

Naive.gene <- c('CCR7','TCF7','LEF1','SELL')
Resident.gene <- c('CD69','RUNX3','NR4A1')
# Inhibitory.gene <- c('TIGIT','CTLA4','LAG3','PDCD1','HAVCR2')
Cytokines <- c('IL2','IL17A','LAMTOR3')
Cyto.gene <- c('GNLY','GZMB','GZMK','GZMA','IFNG','PRF1','NKG7')
Costimulatory <- c('TNFRSF14','TNFRSF9','CD28','ICOS')
Treg <- c('FOXP3','IL2RA')
CellType <- c('KLRF1','KLRD1','CD4','CD3D','CD3E','CD8A','CD8B')
TF <- c('TBX21','ZNF683','ZEB2','ID2','EOMES','HIF1A','TCF7','IRF8')
gene <- c(Naive.gene,Resident.gene,Inhibitory.gene,Cytokines,Cyto.gene,Costimulatory,Treg,CellType,TF)
gene_data <- data.frame(t(as.matrix(Tcell2@assays$RNA@data)[gene,]),cluster=Tcell2@meta.data$celltype,check.names = F)
average_data <- aggregate(.~cluster, gene_data, mean)
cluster_name <- average_data[,1]
average_data <- apply(average_data[,2:ncol(average_data)],2,as.numeric)
rownames(average_data) <- cluster_name
average_data <- t(average_data)
phmat1 <- t(scale(t(average_data)))
phmat1[phmat1> 2] <- 2
phmat1[phmat1 < -2] <- -2
#degs
library(pheatmap)
pheatmap(t(gsva_aggr),cluster_rows = F,cluster_cols = F,fontsize = 8,colorRampPalette(colors = c('#4D67B0','#789FBA','white','brown','#882F21'))(100))


### CD4
CD4 <- Tcell2 [, Tcell2$celltype %in% c('CD4-CD40LG','CD4-Exha','CD4-Naive-IL7R','Treg')]
CD4 <- CD4 [,! CD4$seurat_clusters %in% c('4','5')]

CD4 <- CreateSeuratObject(CD4@assays$RNA@counts[rowSums(as.matrix(CD4@assays$RNA@counts) > 0) > 20,],meta.data = CD4@meta.data)
CD4 <- NormalizeData(CD4)
G1S <- toupper(c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"))
G2M <- toupper(c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3" ,"Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3", "Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"))
CD4 <- CellCycleScoring(CD4,s.features =G1S,g2m.features = G2M )
CD4 <- ScaleData(CD4)
CD4 <-FindVariableFeatures(CD4)
CD4 <- RunPCA(CD4,features = VariableFeatures(CD4),npcs = 30,verbose = F)
library(harmony)
CD4 <- CD4 %>% 
  RunHarmony('orig.ident',plot_convergence=TRUE)
harmony_embeddings <- Embeddings(CD4,'harmony')
ElbowPlot(CD4)
CD4 <- CD4 %>%
  RunTSNE(reduction='harmony',dims=1:30)%>%
  RunUMAP(reduction='harmony',dims=1:30)%>%
  FindNeighbors(reduction='harmony',dims=1:30) %>%
  FindClusters(resolution=0.3) %>%
  identity()
library(dplyr)
CD8 <- CD8 %>%  
  FindNeighbors(reduction='harmony',dims=1:30) %>%
  FindClusters(resolution=0.7) %>%
  identity()


DimPlot(CD8,reduction = 'umap',label = T)

CD4$celltype <- plyr::mapvalues(CD4$seurat_clusters,from = c('2','7','0','9','10','11','4','3','2','8','6','5'),to = c('CD4-Exha','CD4-Exha','CD4-CD40LG','CD4-CD40LG','CD4-Naive-IL7R','CD4-Naive-IL7R','CD4-Naive-IL7R','CD4-Naive-IL7R','CD4-Naive-IL7R','CD4-Naive-IL7R','Treg_TNFRSF9(-)','Treg_TNFRSF9(+)'))
CD4$celltype <- factor(CD4$celltype,levels  = c('CD4-Exha', 'CD4-CD40LG','CD4-Naive-IL7R','Treg_TNFRSF9(-)','Treg_TNFRSF9(+)'))
Idents(CD4) <- 'celltype'
CD4.markers <- FindAllMarkers(CD4,logfc.threshold = log(1.5),test.use = 'wilcox')
library(dplyr)
CD4top10 <- CD4.markers %>% filter(avg_log2FC>0 & p_val<0.05) %>% group_by(cluster) %>% top_n(wt = avg_log2FC,8)


CD8$celltype <- plyr::mapvalues(CD8$seurat_clusters,from = c('0','13','8','1','12','11','2','10','3','4','9','5','7','6'),to = c('CD8-ZNF683','CD8-ZNF683','CD8-ZNF683','CD8-Cyto','CD8-Cyto','CD8-Cyto2','CD8-Cyto2','CD8-Cyto2','CD8-Exha','CD8-IFIT','CD8-Exha2','CD8-Exha2','CD8 Cycling S','CD8 Cycling G2M'))
CD8$celltype <- factor(CD8$celltype,levels  =  c('CD8-ZNF683', 'CD8-Cyto','CD8-Cyto2', 'CD8-Exha','CD8-IFIT','CD8-Exha2', 'CD8 Cycling S','CD8 Cycling G2M'))

CD8$celltype <- plyr::mapvalues(CD8$celltype, from  =  c('CD8-ZNF683', 'CD8-Cyto','CD8-pre-Exha', 'CD8-Exha','CD8-Exha2','CD8-IFIT','CD8 Cycling S','CD8 Cycling G2M'),to  =  c('CD8-ZNF683', 'CD8-APOC1','CD8-Cyto-GZMK', 'CD8-Cyto-PRF1','CD8-Exha-CXCL13','CD8-IFIT', 'CD8 Cycling S','CD8 Cycling G2M'))

Idents(CD8) <- 'celltype'
CD8.markers <- FindAllMarkers(CD8,logfc.threshold = log(1.5),test.use = 'wilcox')
library(dplyr)
CD8top10 <- CD8.markers %>% filter(avg_log2FC>0 & p_val<0.05) %>% group_by(cluster) %>% top_n(wt = avg_log2FC,30)


sue <- merge(CD4,c(CD8,myeloid2))
a <- row.names(pbmc3@meta.data) [!row.names(pbmc3@meta.data) %in% row.names(sue@meta.data)]
sue2 <- pbmc3 [, a ]

library(Seurat)
library(monocle)
CD8 <- CD4 [,CD4$celltype %in% c('CD4-Exha',    'CD4-preExha'  ,   'CD4-CD40LG', 'CD4-Naive-IL7R',           'Treg' )]
pd <- CD8@meta.data [,c( 'celltype' ),drop=F]
pd <- new('AnnotatedDataFrame',pd)
HSMM <- newCellDataSet(as.matrix(CD8@assays$RNA@data),phenoData = pd)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table, mean_expression >0.12 &
                           dispersion_empirical >= 1.0  * dispersion_fit)$gene_id
cell_cycle<-read.csv("D:/cell_cycle_genes.csv")
ordering_genes <- setdiff(ordering_genes, cell_cycle$`Cell cycle related genes(304)`)

HSMM<- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM,max_components=2)
HSMM <- orderCells(HSMM, reverse=FALSE,root_state = 3)
pData(HSMM)
load('D:/lab_cor.Rdata')
p1 <- plot_cell_trajectory(HSMM,color_by = "celltype",show_branch_points = F,cell_size = 0.7)+theme_bw()+scale_color_manual(values = my_cor1)
p2 <- plot_cell_trajectory(HSMM,color_by = 'State',show_branch_points = F,cell_size = 0.7)+theme_bw()+scale_color_manual(values = my_cor1)
HSMM <- orderCells(HSMM, root_state = 3)
p3 <- plot_cell_trajectory(HSMM,color_by = 'Pseudotime',show_branch_points = F,cell_size = 0.7)+theme_bw()+scale_color_gradientn(colours =  colorRampPalette(rev(brewer.pal(10, "Spectral")))(99)) +
  scale_size(range = c(3,3))
CombinePlots(list(p1,p2,p3),ncol = 1)
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
plotdf=pData(HSMM)
library(ggridges)
ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(2.5,5),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+theme(panel.grid = element_blank())


Time_diff <- differentialGeneTest(HSMM[ordering_genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff$gene_short_name <- row.names(Time_diff)
#MT <- grep('^MT-',x = row.names(Time_diff),value = T)
#Time_diff <- Time_diff [!Time_diff$gene_short_name %in% 'MT',]
library(dplyr)
Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
library(RColorBrewer)
Time_genes <- Tcell.markers.top$gene[81:215] [!duplicated(Tcell.markers.top$gene[81:215])]
plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))
plot_genes_branched_heatmap(HSMM[c(Time_genes ),], branch_point = 1,num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))

#### Branch Pattern by monocle2
BEAM_res <- BEAM(HSMM, branch_point = 1, cores =4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res$gene_short_name <- rownames(BEAM_res)
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
library(RColorBrewer)
library(dplyr)
Time_genes <- top_n(BEAM_res, n = 120, desc(qval)) %>% pull(gene_short_name) %>% as.character()
plot_genes_branched_heatmap(HSMM[c(Time_genes,'TREM2','SPP1',  'MARCO'),], branch_point = 1,num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))
plot_genes_branched_heatmap(HSMM[ Time_genes,  ], branch_point = 1,num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))

allTF_human <- read.csv("D:/code/TF_data/allTF.human_v2.txt", header = F) 
BEAM_res2 <- BEAM_res [intersect(BEAM_res$gene_short_name,allTF_human$V1),]
tf <- top_n(BEAM_res2, n = 150, desc(qval)) %>% pull(gene_short_name) %>% as.character()
plot_genes_branched_heatmap(HSMM[tf,], branch_point = 1,num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))

ggsave("Time_heatmap.pdf", p, width = 5, height = 10)

plot_pseudotime_heatmap(HSMM[c(row.names(BEAM_res)[1:100] ),], num_clusters = 3, cores=4, show_rownames=TRUE)


###cellchat
anno1 <- CD4@meta.data[,c('celltype','batch'),drop=F]
anno2 <- myeloid2@meta.data[,c('celltype','batch'),drop=F]
anno3 <- CD8@meta.data[,c('celltype','batch'),drop=F]
anno <- rbind(anno1,anno3,anno2)
anno <- anno [row.names(anno) %in% row.names(pbmc3@meta.data),drop=F,]
pbmc4<- pbmc3[,row.names(anno)]
library(CellChat)
library(Seurat)
load('./pbmc3.Rds')
pbmc <- merge(CD4,c(CD8,myeloid2))

data.input = as.matrix(pbmc@assays$RNA@data) # normalized data matrix
identity = data.frame(group = pbmc@meta.data$cluster, row.names = names(pbmc@meta.data$cluster)) # create a dataframe consisting of the cell labels
unique(identity$group)
cellchat <- createCellChat(object = data.input, group.by = "cluster",meta = pbmc@meta.data )
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
save('./cellchat.manuscript.Rds',cellchat)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

genes <- cellchat@netP
pathways.show.all <- cellchat@netP$pathways
pathways.show <- c("SPP1") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
par(mfrow=c(1,2))
netVisual_aggregate(cellchat, signaling = '', layout = "circle")
par(mfrow=c(3,4))
netVisual_aggregate(cellchat, signaling = 'PD-L1', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'PDL2', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'TIGIT', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'CD80', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'CD86', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'CD226', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'CXCL', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'CD40', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'SPP1', layout = "chord")
netVisual_aggregate(cellchat, signaling = 'IFN-II', layout = "chord")

netVisual_chord_gene(cellchat, sources.use = 15, targets.use = c(10), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = c(1), targets.use = c(5), legend.pos.x = 15)
netVisual_chord_gene(cellchat, sources.use = c(3,4,5,6), targets.use = c(14:18), legend.pos.x = 15)
netVisual_chord_gene(cellchat, sources.use = c( 10,11,13), targets.use = c(14:18), legend.pos.x = 15)
plotGeneExpression(cellchat, signaling = "COLLAGEN")
netVisual_chord_gene(cellchat, sources.use = c(9,10,11), targets.use = c(15),legend.pos.x = 8)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3), targets.use = c(4:5), legend.pos.x = 8)
netVisual_chord_gene(cellchat, targets.use = c(12), sources.use = c(1:6), legend.pos.x = 8)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = 'TNF', geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
netVisual_individual(cellchat, signaling =  'TNF',  pairLR.use = LR.show, layout = "chord")



# Hierarchy plot
vertex.receiver = seq(1,2) # a numeric vector
netVisual_individual(cellchat, signaling = c('IFN-II','TNF','CXCL'),  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = c('IFN-II','TNF','CXCL'), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = c( 'CXCL' ), width = 8, height = 2.5, font.size = 10)


p1 <- netVisual_bubble(cellchat, sources.use = c(3,2,4,5,7,10,11,12,13), targets.use =  c(16:20), signaling = c('IFN-II','CD40','PD-L1','PDL2','TIGIT' ,'CD80','CD86','CD226'), remove.isolate = FALSE)
p2 <- netVisual_bubble(cellchat, targets.use = c(3,2,4,5,7,10,11,12,13),  sources.use=  c(16:20), signaling = c('IFN-II', 'CD40','PD-L1','PDL2','TIGIT' ,'CD80','CD86','CD226'), remove.isolate = FALSE)

p3 <- netVisual_bubble(cellchat, sources.use = c(3,2,4,5,7,10,11,12,13), targets.use =  c(16:20), signaling = c('CCL','CXCL'), remove.isolate = FALSE)
p4 <- netVisual_bubble(cellchat, targets.use = c(3,2,4,5,7,10,11,12,13),  sources.use=  c(16:20), signaling = c('CCL','CXCL' ), remove.isolate = FALSE)


####
data_cpDB <- sue@assays$RNA@counts[rowSums(as.matrix(sue@assays$RNA@counts) > 0) > 20,row.names(anno)]
# eg = bitr(rownames(data_cpDB), fromType="SYMBOL", toType=c("ENSEMBL"), OrgDb="org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
anno <- anno[sample(row.names(anno),10000),drop=F,]
data_cpDB <- data.input
eg <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(data_cpDB), columns = c("SYMBOL", "ENSEMBL"), keytype="SYMBOL")
eg <- eg[!duplicated(eg$SYMBOL), ]
eg <- na.omit(eg)
data_cpDB <- data_cpDB[rownames(data_cpDB) %in% eg$SYMBOL, ]
dim(eg);dim(data_cpDB); table(rownames(data_cpDB) == eg$SYMBOL)
rownames(data_cpDB) <- plyr::mapvalues(rownames(data_cpDB), from = eg$SYMBOL, to = eg$ENSEMBL)
anno <- annot [intersect(colnames(data_cpDB),row.names(annot)),]
anno <- anno[,'celltype',drop=F]
colnames(anno) <- "cell_type"
getwd()
dir.create('./cellphonedb2')
write.csv(as.matrix(data_cpDB), file = "./cellphonedb2/all_cpDB_count.csv")
write.csv(anno, file = "./cellphonedb2/all_cpDB_meta.csv")
# cellphonedb method statistical_analysis ./all_cpDB_meta.csv ./all_cpDB_count.csv --iterations=10 --threads=2




library(Seurat)
library(monocle)
myeloid <- myeloid2 [,!myeloid2$celltype %in% c('DC-CLEC9A','Mac2', 'Mac-Cycling','Mac1')]
pd <- myeloid@meta.data [,c('batch','celltype'),drop=F]
pd <- new('AnnotatedDataFrame',pd)
HSMM <- newCellDataSet(as.matrix(myeloid@assays$RNA@data),phenoData = pd)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.1 &
                           dispersion_empirical >= 1.5 * dispersion_fit)$gene_id
cell_cycle<-read.csv("D:/cell_cycle_genes.csv")
ordering_genes <- setdiff(ordering_genes, cell_cycle$`Cell cycle related genes(304)`)

HSMM<- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM,max_components=2)
HSMM2 <- orderCells(HSMM2, reverse=FALSE)
pData(HSMM)
load('D:/lab_cor.Rdata')
p1 <- plot_cell_trajectory(HSMM,color_by = "celltype",show_branch_points = F,cell_size = 0.7)+theme_bw()+scale_color_manual(values = c('grey','grey','orange','orange', 'purple','purple' ))
p2 <- plot_cell_trajectory(HSMM,color_by = 'State',show_branch_points = F,cell_size = 0.7)+theme_bw()+scale_color_manual(values = my_cor1)
HSMM <- orderCells(HSMM, root_state = 6)
p3 <- plot_cell_trajectory(HSMM,color_by = 'Pseudotime',show_branch_points = F,cell_size = 0.7)+theme_bw()+scale_color_gradientn(colours =  colorRampPalette(rev(brewer.pal(10, "Spectral")))(99)) +
  scale_size(range = c(3,3))
CombinePlots(list(p1,p2,p3),ncol = 1)
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
plotdf=pData(HSMM)
library(ggridges)
ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+theme(panel.grid = element_blank())


Time_diff <- differentialGeneTest(HSMM[ordering_genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff$gene_short_name <- row.names(Time_diff)
#MT <- grep('^MT-',x = row.names(Time_diff),value = T)
#Time_diff <- Time_diff [!Time_diff$gene_short_name %in% 'MT',]
library(dplyr)
Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
library(RColorBrewer)
plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))
plot_genes_branched_heatmap(HSMM[c(Time_genes,'MARCO'),], branch_point = 1,num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))

#### Branch Pattern by monocle2
BEAM_res <- BEAM(HSMM, branch_point = 1, cores =4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res$gene_short_name <- rownames(BEAM_res)
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
library(RColorBrewer)
library(dplyr)
Time_genes <- top_n(BEAM_res, n = 150, desc(qval)) %>% pull(gene_short_name) %>% as.character()
plot_genes_branched_heatmap(HSMM[Time_genes,], branch_point = 1,num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))

allTF_human <- read.csv("D:/code/TF_data/allTF.human_v2.txt", header = F) 
BEAM_res2 <- BEAM_res [intersect(BEAM_res$gene_short_name,allTF_human$V1),]
tf <- top_n(BEAM_res2, n = 150, desc(qval)) %>% pull(gene_short_name) %>% as.character()
plot_genes_branched_heatmap(HSMM[tf,], branch_point = 1,num_clusters=4, show_rownames=T,hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))

ggsave("Time_heatmap.pdf", p, width = 5, height = 10)

plot_pseudotime_heatmap(HSMM[c(row.names(BEAM_res)[1:100] ),], num_clusters = 3, cores=4, show_rownames=TRUE)


data <- pbmc.myeloid3@assays$RNA@data[,rownames(pData(HSMM))]
pData(HSMM)$CellType <- pData(HSMM)$celltype2

pData(HSMM)$NOTCH3 <-  as.matrix(data)["NOTCH3", ]
pData(HSMM)$THY1 <-  as.matrix(data)["THY1", ]
pData(HSMM)$CDKN2A <-  as.matrix(data)["CDKN2A", ]
pData(HSMM)$TP63 <-  as.matrix(data)["TP63", ]
pData(HSMM)$SOX2 <-  as.matrix(data)["SOX2", ]
pData(HSMM)$TFF3 <-  as.matrix(data)["TFF3", ]
library(RColorBrewer)
plot_cell_trajectory(HSMM,color_by = "Havcr1",cell_size = 0.4)+  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))+theme(legend.position = "right")
plot_cell_trajectory(HSMM,color_by = "TP63",cell_size = 0.4)+  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))+theme(legend.position = "right")

p1 <- plot_cell_trajectory(HSMM,color_by = "NOTCH3",cell_size = 0.4)+  scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(6))+theme(legend.position = "right")
p2 <- plot_cell_trajectory(HSMM,color_by = "THY1",cell_size = 0.4)+  scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(6))+theme(legend.position = "right")
p3 <- plot_cell_trajectory(HSMM,color_by = "SPRR1B",cell_size = 0.4)+  scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(6))+theme(legend.position = "right")
p4 <- plot_cell_trajectory(HSMM,color_by = "CEACAM5",cell_size = 0.4)+  scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(6))+theme(legend.position = "right")
p5 <- plot_cell_trajectory(HSMM,color_by = "CDKN2A",cell_size = 0.4)+  scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(6))+theme(legend.position = "right")
p6 <- plot_cell_trajectory(HSMM,color_by = "TFF3",cell_size = 0.4)+  scale_color_gradientn(colours = colorRampPalette(c('grey','red'))(6))+theme(legend.position = "right")

library(GSEABase)
library(GSVA)
geneset2 <- getGmt('D:/code/symbols/h.all.v6.2.symbols.gmt')
pbmc4 <- pbmc3[,sample(row.names(pbmc3@meta.data),10000)]
system.time(gsva_es_all <- gsva(as(pbmc4@assays$RNA@data[rowSums(as.matrix(pbmc4@assays$RNA@counts) > 0) > 20,],"matrix"), geneset2, mx.diff=1,parallel.sz =20))
save(file = './all.gsva.Rds',gsva_es_all)
## fibroblast merge
phmat <- t(scale(t(gsva_es_h)))
phmat [phmat>1.5] <- 1.5
phmat [phmat < -1.5] <- -1.5

anno <- EPI@meta.data
gsva_data <- data.frame(t(phmat[,rownames(anno)]),check.names = F)
gsva_data$cluster <- anno$celltype2
gsva_aggr <- aggregate(.~cluster,gsva_data,mean)
rownames(gsva_aggr) <- gsva_aggr$cluster
gsva_aggr <- gsva_aggr[,-1]
colnames(gsva_aggr) <- gsub(pattern = "HALLMARK_",replacement = "",x=colnames(gsva_aggr))
pheatmap::pheatmap(t(gsva_aggr),color = colorRampPalette(c('navy','white','red'))(100),show_colnames = T)
pheatmap(t(phmat1),cluster_rows = F,cluster_cols = F,fontsize = 8,colorRampPalette(colors = c('#789FBA','white','#882F21'))(100))

####
library(pheatmap)
pheatmap(phmat[, rownames(anno)[order(anno$cluster)]], show_colnames = F, cluster_cols = F,show_rownames = T,clustering_method = 'complete',
         color = colorRampPalette(c('navy','white','red'))(100),annotation_col =anno, treeheight_row = 0)

pheatmap(phmat[, rownames(anno)[order(anno$res.0.3)]], show_colnames = F, cluster_cols = F,show_rownames = T,clustering_method = 'complete',
         color = colorRampPalette(c('navy','white','red'))(100),annotation_col = anno, treeheight_row = 0)

pheatmap(phmat[, rownames(anno)[order(anno$orig.ident)]], show_colnames = F, cluster_cols = F,show_rownames = T,clustering_method = 'complete',
         color = colorRampPalette(c('navy','white','red'))(100),annotation_col = anno, treeheight_row = 0)

###scenic
library(SCENIC)
library(GENIE3)
library(AUCell)
library(foreach)
dir.create('SCENIC')
setwd('./SCENIC/')
dir.create('int')
##准备细胞meta信息
cellInfo <- data.frame(Tcell2@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype")] <- "celltype"
cellInfo <- cellInfo[,c("sample","celltype")]
saveRDS(cellInfo, file="./int/cellInfo.Rds")
##准备表达矩阵
#为了节省计算资源，随机抽取1000个细胞的数据子集
subcell <- sample(colnames(Tcell2),2000)
scRNAsub <- Tcell2[,subcell]
saveRDS(scRNAsub, "scRNAsub.rds")
TCGAexprrMat <- as.matrix(scRNAsub@assays$RNA@counts)
##设置分析环境
mydbDIR <- "D:/cisTarget_databases/"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather",
           "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "lung")
saveRDS(scenicOptions, "./int/scenicOptions.rds")


##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(TCGAexprrMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(TCGAexprrMat), 
                           minSamples = ncol(TCGAexprrMat) * 0.01)
TCGAexprrMat_filtered <- TCGAexprrMat[genesKept, ]
##计算相关性矩阵
runCorrelation(TCGAexprrMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
TCGAexprrMat_filtered_log <- log2(TCGAexprrMat_filtered+1)
runGenie3(TCGAexprrMat_filtered_log, scenicOptions, nParts = 20)
##推断共表达模块
runSCENIC_1_coexNetwork2modules(scenicOptions)
##推断转录调控网络（regulon）
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "lung")
runSCENIC_2_createRegulons(scenicOptions)
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量

##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
TCGAexprrMat_all <- as.matrix(Tcell2@assays$RNA@counts)
TCGAexprrMat_all <- log2(TCGAexprrMat_all+1)
runSCENIC_3_scoreCells(scenicOptions,  TCGAexprrMat_all)

runSCENIC_4_aucell_binarize(scenicOptions,  TCGAexprrMat_all)
##导入原始regulonAUC矩阵

getwd()
scenic <-  readRDS('./myeloid/SCENIC/int/4.2_binaryRegulonActivity_nonDupl.Rds')
dim(scenic)
table(row.names(pbmc.tumor2@meta.data) %in% colnames(scenic))
scenic <- scenic[,rownames(myeloid2@meta.data) %in% colnames(scenic) ]
###sample cell
anno <- myeloid2@meta.data [colnames(scenic),"celltype",drop=F]
cor_data <- cor(t(scenic),method = 'spearman')
hist(rowSums(cor_data>0.2,na.rm = T))
filter_regulon <- rownames(cor_data)[rowSums(cor_data>0.5,na.rm = T)>3]
scenic_data <- data.frame(t(scenic),check.names = F)
scenic_data$cluster <- anno$celltype
scenic_aggr <- aggregate(.~cluster,scenic_data,mean)
cell_num <- table(anno$celltype)
for (i in 1:nrow(scenic_aggr)){
  cl <- scenic_aggr$cluster[i]
  num <- cell_num[cl]
  
  scenic_aggr[i,2:ncol(scenic_aggr)]<- scenic_aggr[i,2:ncol(scenic_aggr)]/num
  
}

rownames(scenic_aggr) <- scenic_aggr$cluster
scenic_aggr <- scenic_aggr[,2:ncol(scenic_aggr)]
phmat <- t(scale(t(scenic_aggr)))
phmat [phmat>1.5] <- 1.5
phmat [phmat < -1.5] <- -1.5
library(pheatmap)
pheatmap(t(phmat),show_rownames = T,show_colnames = T,cluster_rows = T,cluster_cols = F)


ggData <- cbind(myeloid2@meta.data,myeloid2@reductions$umap@cell.embeddings)
data.frame(ggData,IRF1=scenic[grep('IRF1',rownames(scenic)), ]) %>% ggplot(aes(x=UMAP_1, y=UMAP_2,col=IRF1))+geom_point(size=0.1)+
  scale_color_gradientn(colours = colorRampPalette(c('white','red'))(10))+theme_classic()

ggpubr::ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11),ncol = 4,nrow = 3)




expressed_gene <- rownames(pbmc.tumor2@assays$RNA@data[rowSums(as.matrix(pbmc.tumor2@assays$RNA@data) > 1) > 10, ])
diff_he <- differentialGeneTest(HSMM[expressed_gene,], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 4)
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
library(dplyr)
diff_he_sig <- subset(diff_he, qval < 0.01)
diff_he_sig$gene_short_name <- row.names(diff_he_sig)
diff_he_sig <- diff_he_sig[!duplicated(diff_he_sig$gene_short_name), ]
diff_he_sig <- tibble::rownames_to_column(diff_he_sig) %>% dplyr::arrange(plyr::desc(-qval))
Time_genes <- top_n(diff_he_sig, n = 250, desc(qval)) %>% pull(gene_short_name) %>% as.character()

p <- plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters = 3,cores = 1,show_rownames = F, trend_formula  = "~sm.ns(Pseudotime, df=3)",return_heatmap = T)
p

t <- as.data.frame(cutree(p$tree_row, k=3))
colnames(t) <- "Cluster"
t$Gene <- rownames(t)
genes_c1 <- rownames(t[t$Cluster == 1, ])
genes_c2 <- rownames(t[t$Cluster == 2, ])
genes_c3 <- rownames(t[t$Cluster == 3, ])

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
engo <-     enrichGO(gene          = genes_c2,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1)
View(engo@result)
dotplot(engo)+ scale_y_discrete( labels = function(x) str_wrap(x, width = 60))

pbmc.tumor2 <- pbmc.tumor2 [,!pbmc.tumor2$CellType %in% 'EPI_3']
pbmc.tumor2@meta.data$Pseudotime <- pData(HSMM)$Pseudotime
data_smoothed <- t(rowSmooth(t(pbmc.tumor2@data[genes_c2, ]), n = 8, do.scale = T, proportion = 20))
data_smoothed <- MinMax(data_smoothed, min = -1.5, max = 1.5)
pheatmap(data_smoothed[, order(pbmc.tumor2@meta.data[, "Pseudotime"])],
         annotation_colors = ph_colors, fontsize_row = 10, cluster_rows = F, 
         cluster_cols = F, annotation_col = pbmc.tumor2@meta.data[,c("stage","Pseudotime"), drop = F],show_rownames = T,
         show_colnames = F,  color = colorRampPalette(colors = c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"))(100),
         border_color = NA, annotation_legend = T, gaps_row = (c(8,16)))

pbmc.tumor2 <- AddModuleScore(pbmc.tumor2, features = list(genes_c1),  name = "cell_cycle")
pbmc.tumor2 <- AddModuleScore(pbmc.tumor2, features = list(genes_c2), name = "Metaplasima")
pbmc.tumor2 <- AddModuleScore(pbmc.tumor2, features = list(genes_c3), name = "EMT")

ggdata <- pbmc.tumor2@meta.data

ggplot(ggdata, aes(x = Pseudotime)) + 
  geom_smooth(aes(y = cell_cycle1), method = "auto", size = 1.5, color = "#30A9DE") +
  geom_smooth(aes(y = Metaplasima1), method = "auto", size = 1.5, color = "#A593E0") +
  geom_smooth(aes(y = EMT1), method = "auto", size = 1.5, color = "#C16200") + 
  labs(x= "Pseudotime",  y = "Features") + 
  annotate("text", x=1.2, y=0.25, label="cell cycle", fontface="italic", colour="#30A9DE", size=5)+
  annotate("text", x=7.5, y=0.5, label="Metaplasima", fontface="italic", colour="#A593E0", size=5)+
  annotate("text", x=16, y=0.8, label="EMT", fontface="italic", colour="#C16200", size=5)+
  theme(axis.text = element_blank(), axis.ticks = element_blank(),legend.title=element_blank(),legend.position = c(0.9, 0.9))



## HE correlated to SPRR3
### genes correlated with SPRR3
# source("D:/data/bzj2018human 小胚数据/TACS_analysis_plot.R")
# TACS <- TACS(dge = pbmc_DA, gene1 = "SPRR3", gene2 = "CD44", cutoffs = c(0, 0), return_val = "all", density = T, col = ph_colors$DAcluster,
#              facet_by = "DAcluster", num_genes_add = 100, genesets_predetermined = F, dge_reference = pbmc_DA)
# TACS$plot
cor_matrix <- t(as.matrix(pbmc.tumor2@assays$RNA@data))
cor_matrix <- cor_matrix[, colSums(cor_matrix > 0) > 3]
# cor_matrix_SPRR3 <- cor(TACS$plot_df[, "SPRR3_score", drop = F], cor_matrix, method = "spearman")
cor_matrix_SPRR3 <- cor(cor_matrix[,"SPRR3"], cor_matrix, method = "spearman")
cor_matrix_SPRR3 <- t(cor_matrix_SPRR3)
# matrix_magic_t_scaled <- t(apply(cor_matrix, 2, scale))
# matrix_magic_t_scaled <- MinMax(matrix_magic_t_scaled, min = -2, max = 2)
tmp <- cor_matrix_SPRR3[cor_matrix_SPRR3 > 0.4 | cor_matrix_SPRR3 < -0.4, ]
length(tmp)
sort(tmp)
genes_cor_SPRR3 <- sort(cor_matrix_SPRR3[,1], decreasing = T)[2:40]
genes_uncor_SPRR3 <- sort(cor_matrix_SPRR3[,1], decreasing = F)[1:40]

ggdata_cor <- data.frame(order = 1:79, gene_cor = c(names(genes_cor_SPRR3),rev(names(genes_uncor_SPRR3))), 
                         Correlation_coefficient = c(genes_cor_SPRR3, rev(genes_uncor_SPRR3)))
library(ggrepel)
ggplot(ggdata_cor, mapping = aes(x = order, y = Correlation_coefficient,label = gene_cor)) + 
  geom_bar(width = 0.9, stat = "identity", fill = "lightblue") +
  # geom_label_repel()+
  labs(x = NULL, y = "Correlation_coefficient", title = "Genes related with SPRR3") +  
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(), 
        strip.background = element_rect(fill = NA, colour = NA),
        strip.placement = "inside",
        strip.text = element_text(face = "bold.italic")) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.5) +
  geom_text(size = 4,fontface="bold.italic", angle = 90)


RPS_gene <- grep("RP", names(sort(cor_matrix_SPRR3[,1], decreasing = T))[1:500], value = T)


library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
engo <-     enrichGO(gene          = names(sort(cor_matrix_SPRR3[,1], decreasing = F))[1:100],
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1)
View(engo@result)
write.csv(engo@result, file = "./cor_SPRR3_GO.csv")
write.csv(engo@result, file = "./uncor_SPRR3_GO.csv")


mypie_plot <- function(pbmc, group.by = NULL, split.by = NULL, col.group = NULL, ncol = NULL){
  pie_data <- as.data.frame(prop.table(x = table(pbmc@meta.data[[group.by]], pbmc@meta.data[[split.by]]), margin = 2))
  colnames(pie_data) <- c("cluster", "group",  "proportion")
  pie_data$cluster <- factor(pie_data$cluster, levels = levels(pbmc@meta.data[[group.by]]))
  pie_data$group <- factor(pie_data$group, levels = levels(pbmc@meta.data[[split.by]]))
  ggplot(data=pie_data, mapping=aes_string(x="1", y = "proportion", fill="cluster"))+
    geom_bar(stat="identity", width=1,position='stack', size=1)+
    coord_polar("y", start=0)+
    facet_wrap(~group, ncol = ncol)+
    scale_fill_manual(values=col.group)+
    #geom_text(stat="identity",aes_string(y="proportion", label = paste(round(proportion*100, 1), "%")), size=3, position=position_stack(vjust = 0.5)) +
    theme_minimal()+ 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(size = 12))
}

load('./all.gsva.Rds')
anno1 <- CD4@meta.data[,'celltype',drop=F]
anno2 <- myeloid2@meta.data[,'celltype',drop=F]
anno3 <- CD8@meta.data[,'celltype',drop=F]
anno <- rbind(anno1,anno3,anno2)
load('./pbmc3.Rds')
gsva <- data.frame(t(gsva_es_all[,intersect(row.names(anno),colnames(gsva_es_all))]))
anno2 <- anno[intersect(row.names(anno), colnames(gsva_es_all)),drop=F,]
tpm <- cbind(anno2,gsva)
colnames(tpm) <- gsub(pattern = "HALLMARK_",replacement = "",x=colnames(tpm))
library(ggpubr)
source('D:/code/function.R')
p1 <- ggplot(data =tpm, aes(x=celltype, y=TNFA_SIGNALING_VIA_NFKB,fill=celltype ))+ geom_violin() +scale_fill_manual(values =my36colors) + stat_compare_means() +theme(axis.text.x = element_text(angle = 45, hjust = 1)) +theme(legend.position="none") +geom_boxplot(width=0.1,position=position_dodge(0.9))
p2 <- ggplot(data =tpm, aes(x=celltype, y=IL6_JAK_STAT3_SIGNALING,fill=celltype ))+ geom_violin() +scale_fill_manual(values =my36colors) + stat_compare_means() +theme(axis.text.x = element_text(angle = 45, hjust = 1)) +theme(legend.position="none") +geom_boxplot(width=0.1,position=position_dodge(0.9))
p1+p2



annot <- myeloid2@meta.data
cell_num <- table(as.character(annot$orig.ident))
num_table <- table(annot$orig.ident, annot$celltype)[names(cell_num),]
num_table <- data.frame(num_table/as.numeric(cell_num))
names(num_table) <- c('sample','cluster','proportion')

library(tidyr)
wide3<-spread(num_table,cluster,proportion)
row.names(wide3) <- wide3$sample
wide3 <- wide3 [,-1]
wide <- read.csv('./w.csv')
row.names(wide) <- wide$X
wide <- wide[,-1]
library(ggpubr)
p1 <- ggscatter(wide, x = "DC.LAMP3", y = "CD4.Exha.CXCL13",
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson",   label.sep = "\n"))

p2 <- ggscatter(wide, x = "DC.LAMP3", y = "CD4.IFIT",
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson",   label.sep = "\n"))
p3 <- ggscatter(wide, x = "DC.LAMP3", y = "Mono1",
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson",   label.sep = "\n"))









anno1 <- anno [anno$batch == 'CIP(-)',]
anno2 <- anno [anno$batch == 'CIP(+)',]
anno3 <- anno [anno$batch == 'Healthy',]
save(file='./cellchatprepare.Rds',End,Fib,pbmc.myeloid2,pbmc.Tcell2,pbmc.tumor2,anno)

data.input1 =  pbmc4@assays$RNA@data[,rownames(anno1)] # normalized data matrix
cellchat1 <- createCellChat(object = data.input1 , meta = anno1 , group.by = "celltype")
cellchat1 <- setIdent(cellchat1, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat1@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human# use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB
cellchat1@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat1 <- subsetData(cellchat1) # This step is necessary even if using the whole database
future::plan("multiprocess", workers =1) # do parallel
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
# project gene expression data onto PPI network (optional)
cellchat1 <- projectData(cellchat1, PPI.human)
cellchat1 <- computeCommunProb(cellchat1)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
groupSize <- as.numeric(table(cellchat1@idents))


data.input2=  pbmc4@assays$RNA@data[,rownames(anno2)] # normalized data matrix
cellchat2 <- createCellChat(object = data.input2 , meta = anno2 , group.by = "celltype")
cellchat2 <- setIdent(cellchat2, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat2@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human# use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB
cellchat2@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat2 <- subsetData(cellchat2) # This step is necessary even if using the whole database
future::plan("multiprocess", workers =1) # do parallel
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
# project gene expression data onto PPI network (optional)
cellchat2 <- projectData(cellchat2, PPI.human)
cellchat2 <- computeCommunProb(cellchat2)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
groupSize <- as.numeric(table(cellchat2@idents))

pbmc5 <- pbmc3[,pbmc3$batch%in%'Healthy']
data.input3=  pbmc5@assays$RNA@data # normalized data matrix
cellchat3 <- createCellChat(object = data.input3 , meta = pbmc5@meta.data[,'celltype',drop=F] , group.by = "celltype")
cellchat3 <- setIdent(cellchat3, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat3@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human# use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB
cellchat3@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat3 <- subsetData(cellchat3) # This step is necessary even if using the whole database
future::plan("multiprocess", workers =1) # do parallel
cellchat3 <- identifyOverExpressedGenes(cellchat3)
cellchat3 <- identifyOverExpressedInteractions(cellchat3)
# project gene expression data onto PPI network (optional)
cellchat3 <- projectData(cellchat3, PPI.human)
cellchat3 <- computeCommunProb(cellchat3)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat3 <- filterCommunication(cellchat3, min.cells = 10)
cellchat3 <- computeCommunProbPathway(cellchat3)
cellchat3 <- aggregateNet(cellchat3)
groupSize <- as.numeric(table(cellchat3@idents))



object.list <- list(CIP1 = cellchat1, CIP2= cellchat2, Healthy= cellchat3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(3,1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(3,1,2), measure = "weight")
gg1 + gg2

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)+coord_flip()
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


library(ComplexHeatmap)
i=1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

pathways.show <- c("TNF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(4), targets.use = c(7,10,11,15,19,20,21,22),slot.name = "netP", title.name = paste0("Signaling pathways sending from EPI - ", names(object.list)[i]), legend.pos.x = 10)
}



## chemokines
data <- read.table('D:/cpf/chemokines.txt')
data$group <- substr(row.names(data),1,2)
library(ggpubr)
data$group <- plyr::mapvalues(data$group,from = c('NO','CI'),to = c('CIP(-)','CIP(+)'))
data$group <- factor(data$group,levels =  c('CIP(-)','CIP(+)'))
p1 <- ggboxplot(data , x = "group", y = "CXCL10", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p2 <- ggboxplot(data , x = "group", y = "CCL21", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p3 <- ggboxplot(data , x = "group", y = "CXCL13", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p4 <- ggboxplot(data , x = "group", y = "CCL27", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p5 <- ggboxplot(data , x = "group", y = "CXCL5", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p6 <- ggboxplot(data , x = "group", y = "CXCL11", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p7 <- ggboxplot(data , x = "group", y = "CXCL24", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p8 <- ggboxplot(data , x = "group", y = "CXCL26", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p9 <- ggboxplot(data , x = "group", y = "CX3CL1", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p10 <- ggboxplot(data , x = "group", y = "CXCL6", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p11 <- ggboxplot(data , x = "group", y = "GM_CSF", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p12 <- ggboxplot(data , x = "group", y = "CXCL1", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p13 <- ggboxplot(data , x = "group", y = "CXCL10", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p14 <- ggboxplot(data , x = "group", y = "CXCL2", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p15 <- ggboxplot(data , x = "group", y = "CCL1", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p16 <- ggboxplot(data , x = "group", y = "IFN_γ", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p17 <- ggboxplot(data , x = "group", y = "IL_1β", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p18 <- ggboxplot(data , x = "group", y = "IL_2", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p19 <- ggboxplot(data , x = "group", y = "IL_4", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p20 <- ggboxplot(data , x = "group", y = "IL_6", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p21 <- ggboxplot(data , x = "group", y = "CXCL8", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')


p22 <- ggboxplot(data , x = "group", y = "IL_10", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p23 <- ggboxplot(data , x = "group", y = "IL_16", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p24 <- ggboxplot(data , x = "group", y = "CXCL11", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p25 <- ggboxplot(data , x = "group", y = "CCL8", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p26 <- ggboxplot(data , x = "group", y = "CCL2", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p27 <- ggboxplot(data , x = "group", y = "CCL7", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p28 <- ggboxplot(data , x = "group", y = "CCL13", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p29 <- ggboxplot(data , x = "group", y = "CCL22", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p30 <- ggboxplot(data , x = "group", y = "MIF", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p31 <- ggboxplot(data , x = "group", y = "CXCL9", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p32 <- ggboxplot(data , x = "group", y = "CCL3", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')


p33 <- ggboxplot(data , x = "group", y = "CCL15", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p34<- ggboxplot(data , x = "group", y = "CCL20", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p35 <- ggboxplot(data , x = "group", y = "CCL19", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p36 <- ggboxplot(data , x = "group", y = "CCL23", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p37 <- ggboxplot(data , x = "group", y = "CXCL16", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p38 <- ggboxplot(data , x = "group", y = "CXCL12", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p39 <- ggboxplot(data , x = "group", y = "CCL17", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p40 <- ggboxplot(data , x = "group", y = "CCL25", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')
p41 <- ggboxplot(data , x = "group", y = "TNFa", palette  = "jco",color = 'group') +theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_compare_means()+ scale_y_log10() +theme(legend.position = 'none')

ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,  p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,p40,p41,nrow = 4,ncol = 10)
