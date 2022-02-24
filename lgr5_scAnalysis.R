library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(SeuratDisk)
library(tidyverse)
library(readxl)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)


#Do not re run SCT
pc <- LoadH5Seurat("pavan_sc_cellid_umap.h5Seurat")


Idents(pc) <- "gs_prediction"
#subset crypt cells and perform differential expression
crypt <- subset(pc, idents = "Crypt cells")
DefaultAssay(crypt) <- "RNA"
lgr5_crypt <- subset(crypt, Lgr5 > 0)
lgr5_crypt$orig.ident <- as.factor(lgr5_crypt$orig.ident)
levels(lgr5_crypt$orig.ident)
DimPlot(lgr5_crypt)
levels(lgr5_crypt$orig.ident) <- rep(c('control', 'gvhd', 'naive'), each = 2)
DimPlot(lgr5_crypt, group.by = 'orig.ident')

crypt$orig.ident <- as.factor(crypt$orig.ident)
levels(crypt$orig.ident)
DimPlot(crypt)
levels(crypt$orig.ident) <- rep(c('control', 'gvhd', 'naive'), each = 2)
DimPlot(crypt, group.by = 'orig.ident')


DefaultAssay(lgr5_crypt) <- "SCT"
Idents(lgr5_crypt) <- 'orig.ident'
condition.markers <- FindAllMarkers(lgr5_crypt, logfc.threshold = 0, min.pct = .25, test.use = 'MAST')
condition.markers <- filter(condition.markers, p_val_adj < .05)
table(condition.markers$cluster)
write.csv(condition.markers,'crypt_condition_markers.csv', quote = F, row.names = F)

#DefaultAssay(lgr5_crypt) <- "SCT"
Idents(lgr5_crypt) <- 'orig.ident'
gvhd.genes <- FindMarkers(lgr5_crypt, ident.1 = 'gvhd', ident.2 = 'control', 
                          logfc.threshold = 0, min.pct = 0.25, test.use = 'MAST')

gvhd.genes <- arrange(gvhd.genes, avg_log2FC)
gvhd.genes <- rownames_to_column(gvhd.genes, var = "gene")
gvhd.genes <- filter(gvhd.genes, p_val_adj < 1e-1)
dim(gvhd.genes)

write.csv(gvhd.genes, "lgr5cryptCellsDE.csv", quote = F, row.names = F)
with(gvhd.genes, plot(-log10(p_val_adj)~ avg_log2FC))

head(gvhd.genes)
pdf('VlnPlot_DEgenes_crypt_lgr5.pdf', width = 12)
VlnPlot(lgr5_crypt, c('Reg3b', 'Reg3g','Lgals2','Calm1','Olfm4'))
dev.off()

##################

gvhd.genes <- read.csv('lgr5cryptCellsDE.csv')

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

k <- keys(txdb, keytype = "GENEID")
ensembl <- useMart(biomart = "ensembl",dataset="mmusculus_gene_ensembl")
map <- getBM(attributes=c('entrezgene_id','mgi_symbol'),
             filters = 'entrezgene_id',
             values = k,
             mart = ensembl)
head(map)

gids <- merge(gvhd.genes,map,by.x='gene',by.y='mgi_symbol')

glist.down <- gids[gids$avg_log2FC < 0,]
glist.up <- gids[gids$avg_log2FC > 0,]

GOdf <- data.frame(Entrez=c(glist.down$entrezgene_id, glist.up$entrezgene_id),
                   group = c(rep("down",length(glist.down$entrezgene_id)),rep("up",length(glist.up$entrezgene_id))))

head(GOdf)
#GO enrichment analysis
formula_res <- compareCluster(Entrez~group,data=GOdf,
                              fun="enrichGO",
                              ont = "BP",
                              pvalueCutoff  = 1e-10,
                              OrgDb = 'org.Mm.eg.db',
                              pAdjustMethod = "BH")

dotplot(formula_res, showCategory = 10, title = "GO Enrichment Analysis")


go_data <- formula_res@compareClusterResult


et_genes <- unlist(strsplit(go_data[go_data$Description == "electron transfer activity",]$geneID,'/'))
pt_genes <- unlist(strsplit(go_data[go_data$Description == "proton transmembrane transporter activity",]$geneID,'/'))
atp_genes <- unlist(strsplit(go_data[go_data$Description == go_data$Description[6],]$geneID,'/'))
pc_genes <- unlist(strsplit(go_data[go_data$Description == "proton channel activity",]$geneID,'/'))


oxphos_test <- unlist(strsplit(go_data[go_data$Description == "oxidative phosphorylation",]$geneID,'/'))
combined <- c(et_genes,pt_genes,atp_genes,pc_genes)
oxphos_go <- unique(combined)

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
oxphos_all <- getBM(attributes=c('entrezgene_id','mgi_symbol','go_id'),
                   filters = 'go', values = 'GO:0006119', mart = ensembl)

#oxphos_all[oxphos_all$go_id == 'GO:0006119',]

oxphos_all <- read_excel('oxphos.xlsx')

gids <- map[map$entrezgene_id %in% oxphos_test,]$mgi_symbol

oxphos_genes_sig <- intersect(gids,oxphos_all$Symbol)

write.table(gids,'oxphos_sig_final.txt',row.names=F,quote=F)

check <- gvhd.genes[gvhd.genes$gene %in% oxphos_genes_sig,]
