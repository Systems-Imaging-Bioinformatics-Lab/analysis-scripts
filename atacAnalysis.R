#load libraries
#BiocManager::install('GenomicRanges')
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#BiocManager::install("org.Mm.eg.db")
#install.packages('vsn')
library("vsn")
library(ggplot2)
library(ChIPseeker)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DESeq2)
library(dplyr)
library(tibble)
library(purrr)
library(GenomicRanges)
library(GenomicFeatures)
library("biomaRt")
library(readxl)
library(ReactomePA)
library(stringr)
library(clusterProfiler)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#----------------------------------------------------------------
#get required supporting data

k <- keys(txdb, keytype = "GENEID")
ensembl <- useMart(biomart = "ensembl",dataset="mmusculus_gene_ensembl")
map <- getBM(attributes=c('entrezgene_id','mgi_symbol'),
             filters = 'entrezgene_id',
             values = k,
             mart = ensembl)
head(map)

#----------------------------------------------------------------
#read counts data 
peak_signal <- read.table("gvhd_all_readCounts.tab",
                          sep = "\t", header = T, comment.char = "")
names(peak_signal)[1] <- "chr"
#peak_signal$chr <- as.factor(paste("chr", peak_signal$chr, sep = ""))
peak_signal$start <- as.numeric(peak_signal$start); peak_signal$end <- as.numeric(peak_signal$end)
head(peak_signal)


#read metadata
cdat <- read.csv("combined_metadata.csv")
cts <- as.matrix(peak_signal[, 4:ncol(peak_signal)])
cdat$samples <- colnames(cts)
all(cdat$sample == colnames(cts))

all(rownames(coldata) == colnames(cts))

#prepare dds object for generating pca plot
dds <- DESeqDataSetFromMatrix(cts, cdat, ~ condition)
dds
vsd <- vst(dds)
meanSdPlot(assay(vsd))

#----------------------------------------------------------------
#generate pca plot
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = "pcaPlotMergedPeak.pdf", width = 8, height = 6, onefile = F)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

#----------------------------------------------------------------
#perform differential analysis
# cdat <- filter(cdat, cell %in% c("LNCaP","V16D") & condition %in% c("E2F1", "EV"))
# cts <- cts[, cdat$sample]
# dim(cts)
# 
# dds <- DESeqDataSetFromMatrix(cts, cdat, ~ cell + condition)
dds <- estimateSizeFactors(dds)
#means <- rowMeans(counts(dds, normalized=TRUE))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#filtering
#extract count matrix, calculate row means

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, alpha = .05, 
               contrast=c("condition","gvhd","control"))
summary(res)
plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup="condition")

#----------------------------------------------------------------
#filter results
res_tb <- res %>% data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble()
# res_tb <- res_tb[order(res_tb$padj), ]
res.sig <- res_tb %>% dplyr::filter(padj < 5e-2)
dim(res.sig)
sum(res.sig$log2FoldChange > 1)
up_sites <- filter(res.sig, log2FoldChange > 0)
down_sites <- filter(res.sig, log2FoldChange < 0)

#----------------------------------------------------------------
#enrichment analysis for upregulated sites
#first create genomic ranges object for peaks
#sites <- down_sites

sites <- up_sites

sites$gene
# peaks1 <- peak_signal[6:11]
# peaks2 <- peak_signal[19:24]
# peaks3 <- peak_signal[1:3]
# peaks <- cbind(peaks3,peaks1,peaks2)

peaks <- peaks[sites$gene,]
head(peaks)

peaks <- peak_signal[sites$gene, 1:3]
head(peaks)

gr <- GRanges(seqnames = peaks$chr, 
              ranges = IRanges(peaks$start, peaks$end),
              "score" = sites$log2FoldChange,
              "padj" = sites$padj,
              "baseMean" = sites$baseMean)

gr <- sort(gr)
gr

# pdf(file = "widthUpregRegions.pdf",
#     width = 6,height = 6, onefile = F)
# hist(width(gr))
# dev.off()

pdf(file = "widthDownRegions.pdf",
    width = 6,height = 6, onefile = F)
hist(width(gr))
dev.off()


#----------------------------------------------------------------
#motif enrichment analysis at upregualted sites
df <- data.frame(gr)[, 1:3]
df$uniqID <- 1:nrow(df)
write.table(df, "downregRegions_gvhd.bed", sep = "\t", quote = F, row.names = F, col.names = F)

#----------------------------------------------------------------
#annotate peaks
peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 1500),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
# 
# pdf(file = "annotatePieUpregRegionsLNCap.pdf",
#     width = 6,height = 6, onefile = F)
# plotAnnoPie(peakAnno)
# dev.off()

pdf(file = "annotatePieDownregRegionsLNCap.pdf",
    width = 6,height = 6, onefile = F)
plotAnnoPie(peakAnno)
dev.off()

gene <- seq2gene(gr, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
head(gene)
symbol <- filter(map, entrezgene_id %in% gene)
dim(symbol)
write.table(symbol$mgi_symbol, "downregulated_genesGVHD.txt", quote = F, row.names = F)


length(which(hyper_genes_atac$V1 %in% symbol$mgi_symbol))
#overlap with genes of interest
targetSet <- read.csv("AR-gene_signature.csv")
head(targetSet)
noquote(formatC(intersect(targetSet$Gene.Symbol, symbol$hgnc_symbol)))

#----------------------------------------------------------------
#look for enriched pathways
x <- enrichPathway(gene, organism = "mouse")

pdf(file = "enrichPathwayUpregSignals.pdf",
    width = 8,height = 6, onefile = F)
p <- dotplot(x)
p <- p + scale_y_discrete(labels=function(x) str_wrap(x, width=30))
p
dev.off()

#----------------------------------------------------------------
#get data from all peaks for pre-ranked GSEA
gr <- GRanges(seqnames = peak_signal$chr, 
              ranges = IRanges(peak_signal$start, peak_signal$end),
              "score" = res_tb$log2FoldChange,
              "padj" = res_tb$padj,
              "baseMean" = res_tb$baseMean)
gr <- sort(gr)
gr
peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 1500),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
anno <- peakAnno@anno

ps <- anno[grep("Promoter", anno$annotation), ]
ps <- distinct(data.frame(ps), SYMBOL, .keep_all = T)
dim(ps)
# 
# df <- data.frame(ps)[, 1:3]
# df$uniqID <- 1:nrow(df)
# write.table(df, "MR42Ddistal_intergenic.bed", sep = "\t", quote = F, row.names = F, col.names = F)
# 
# with(ps, plot(-log10(padj) ~ score))
# ps$rank <- -log10(ps$padj) * ps$score
# ps <- ps[order(ps$rank, decreasing = T), ]
# head(ps)
# 
# #this is my dataset for performing pre-ranked gsea
# write.table(ps, "data_preRankedGSEA.txt", sep = "\t", quote = F)