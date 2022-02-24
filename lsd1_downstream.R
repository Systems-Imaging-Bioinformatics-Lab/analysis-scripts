#Motif analysis, visulaizations for LSD1 project
library(readxl)
#BiocManager::install('chipenrich')
library(chipenrich)
#BiocManager::install("universalmotif")
library(universalmotif)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(rstatix)

#Motif analysis
rnaseq <- read_excel('DEGs.xlsx',sheet=2)

rnaseq_upreg <- subset(rnaseq,rnaseq$FC.easy >= 1.5 & rnaseq$gene_biotype == 'protein_coding')

write.table(rnaseq_upreg$gene,'upreg_rnaseq.txt',quote=F,row.names = F,col.names = F)

data(peaks_E2F4, package = 'chipenrich.data')
head(peaks_E2F4)


data <- read.table('hyper_peaks_lncap_sp2509_dmso_48hr_polyenrich.bed')

write.table(data,'hyper_peaks_lncap_sp2509_dmso_48hr_polyenrich.bed',sep='\t',row.names = F, col.names = F, quote=F)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



jas <- read_jaspar('JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt')


map <- data.frame()

for(i in 1:length(jas)){
  names <- c(jas[[i]]@name,jas[[i]]@altname)
  map <- rbind(map,names)
}

colnames(map) <- c('id','name')


test <- convert_motifs(jas , class = "universalmotif-universalmotif")

write_homer(test,'homer_motifs.motif',overwrite = T)


known_motifs <- read.table('knownResults_hypo_silsd1_ntc_sig.txt',skip=1)

final <- merge(known_motifs,map,by.x='V1',by.y='id')

all_motifs <- final[order(final$V3, decreasing = F),]

cols <- c('MotifID','MotifName','ConsensusP-value','Log P-value','q-value (Benjamini)',	'# of Target Sequences with Motif(of 554)', '% of Target Sequences with Motif','# of Background Sequences with Motif(of 48304)','% of Background Sequences with Motif','name')

colnames(all_motifs) <- cols

write.table(all_motifs,'known_motifs_hypo_silsd1_ntc_sig_JASPAR.txt',sep='\t',quote=F, row.names = F)


hyper <- read.csv('known_motifs_hyper_silsd1_ntc_sig_JASPAR.txt',sep='\t')
hypo <-  read.csv('known_motifs_hypo_silsd1_ntc_sig_JASPAR.txt',sep='\t')

hyper_rank <- hyper[,c('MotifID','ConsensusP.value','name')]
hypo_rank <- hypo[,c('MotifID','ConsensusP.value','name')]

merged <- merge(hyper_rank,hypo_rank,by.x='MotifID',by.y='MotifID')

merged[,c(2,4)] <- log10(merged[,c(2,4)])

merged$diff <- merged$ConsensusP.value.y - merged$ConsensusP.value.x

merged_sorted <- merged[order(merged$diff,decreasing = T),]
merged_sorted$rank <- c(1:nrow(merged_sorted))


options(ggrepel.max.overlaps = Inf)

merged_sorted$direction <- as.factor(ifelse(merged_sorted$diff > 2, "siLSD1", ifelse(merged_sorted$diff < -2,"NTC","n.s")))

p <- ggplot(merged_sorted, aes(x=rank, y=diff))+ geom_point(mapping = aes(color=direction) )+ geom_text_repel(data=subset(merged_sorted, abs(diff) > 2),
            aes(rank,diff,label=name.x)) + labs(x='Motif',y='Global Rank(Differential log p-value)')

write.table(merged_sorted,'ranked_motifs_diffp_silsd1_ntc.txt',sep='\t',quote=F, row.names = F)

########
#Visualizations

sp2509_peaks <- read.csv('SP2509_48_peak_percentage.txt',sep='\t')
dmso_peaks <- read.csv('DMSO_48_peak_percentage.txt',sep='\t')

sp2509_peaks_split <- sp2509_peaks %>%
  separate(col = Feature.Frequency, into=c("1","2","perc"),sep= " ") %>%
  unite(col = "annotation", 1, 2, sep='_')

n_sp2509_peaks <- 138343
  
n_dmso_peaks <- 142823

dmso_peaks_split <- dmso_peaks %>%
  separate(col = Feature.Frequency, into=c("1","2","perc"),sep= " ") %>%
  unite(col = "annotation", 1, 2, sep='_')

pvals <- c()
for(i in 1:nrow(sp2509_peaks_split)){
  perc1 <- round(as.numeric(sp2509_peaks_split[i,'perc']),2)
  perc2 <- round(as.numeric(dmso_peaks_split[i,'perc']),2)
  n1 <- round(perc1 * n_sp2509_peaks / 100)
  n2 <- round(perc2 * n_dmso_peaks / 100)
  A <- matrix(c(n1, n2, n_sp2509_peaks-n1, n_dmso_peaks-n2), nrow = 2)
  print(sp2509_peaks_split[i,'annotation'])
  res <- fisher.test(A)
  print(res$p.value)
}


bp <- ggbarplot(
  dmso_peaks_split , x = "annotation", y = "perc", palette = "npg"
)

dmso_peaks_split$perc <- as.numeric(dmso_peaks_split $perc) 

final_df <- rbind(dmso_peaks_split ,sp2509_peaks_split)
final_df$cond <- c(rep('DMSO',11),rep('sp2509',11))

final_df$perc <- as.numeric(final_df$perc)

p <- ggplot(data=final_df, aes(x=annotation, y=perc, fill = cond)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
 
p + scale_fill_manual(values=c('#999999','#E69F00')) + theme(axis.text.x = element_text(angle = 90)) +
  geom_text(data=final_df, aes(x=1,y=5,label="***"), size=4) +
  geom_text(data=final_df, aes(x=9,y=20,label="***"), size=4)



















