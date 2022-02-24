# BiocManager::install("ChIPseeker")
# BiocManager::install("clusterProfiler")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("ReactomePA")
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(org.Mm.eg.db)

# Load data
samplefiles <- list.files(pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("G34R","WT")

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-5000, 5000), verbose=FALSE)

toMatch <- c("Exon","Intron")
#function to extract peaks in promoter sites
extractGB <- function(w){
  ind <- grep(paste(toMatch,collapse="|"), w$annotation)
  return(w[ind, ])
}

extractIG <- function(w){
  ind <- grep(pattern = "Intergenic", w$annotation)
  return(w[ind, ])
}

extractPromoter <- function(w){
  ind <- grep(pattern = "Promoter", w$annotation)
  return(w[ind, ])
}

peakAnnoList <- lapply(peakAnnoList, function(x) extractPromoter(as.data.frame(x)))

lapply(peakAnnoList, function(i) dim(i))

# Create a list with genes from each sample
genes <- lapply(peakAnnoList, function(i) unique(i$geneId))
lapply(genes, function(i) length(i))

#keep only genes that are uniquely enriched in each condition
glist.wt <- setdiff(genes$WT, genes$G34R)
glist.t <- setdiff(genes$G34R, genes$WT)

write.table(as.numeric(glist.wt),"H3K27me3_G34R_WTMunique_promoter.txt",sep="\n",row.names=FALSE)
write.table(as.numeric(glist.t),"H3K27me3_WT_G34Runique_promoter.txt",sep="\n",row.names=FALSE)






GOdf <- data.frame(Entrez=c(glist.wt, glist.t),
                   group = c(rep("WT",length(glist.wt)),rep("G34R",length(glist.t))))

head(GOdf)
#GO enrichment analysis
formula_res <- compareCluster(Entrez~group,data=GOdf,
                              fun="enrichGO", 
                              OrgDb="org.Mm.eg.db",
                              pvalueCutoff  = 0.05,
                              pAdjustMethod = "BH")

dotplot(formula_res, showCategory = 20, title = "GO Enrichment Analysis")


