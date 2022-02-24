library(stringr)

go.analysis <- read.table("down_analysis.txt",
                          skip = 11, header = T, sep = "\t")

#only looking at pathways with > 10 genes
go <- go.analysis[go.analysis[,3] > 10, ]
# go <- go[order(go$upload_1..FDR.), ]
#top 10 most enriched pathways

# go<-go.analysis
go <- go[1:10,]

go <- go[order(go$upload_1..FDR., decreasing = T),]



labels <- as.character(go$GO.biological.process.complete)
labels <- sapply(strsplit(labels, " (GO", fixed = T),'[',1)

labels <- str_wrap(labels, width = 40)

pdf(file = "GVHD_depleted_pathways.pdf", 
    width = 30, onefile = F, pointsize = 16,
    height = 20)
    

linch <-  max(2*strwidth(labels, "inch"), na.rm = TRUE)
par(mai=c(1.2,linch, 1.5,1.2))


cntrs <- barplot(go[,3],col = 'black',
                 horiz = T, #width = 20, xlim = c(0, 15),
                 names.arg = labels, space = 1,
                 cex.names = 1.75, las = 1.5,  
                 font.axis = 2, cex.axis = 1.5)

mtext(text = "Gene Counts", side = 1, line = 2.5, font = 2, cex = 2)
ylim0 <- par()$usr[3:4] 
par(new = T)
plot.new()
xrng <- range(-log10(go$upload_1..FDR.))
plot.window(ylim = ylim0, xlim = c(xrng[1]-2, xrng[2]+2),
            yaxs = "i")
points(x = -log10(go$upload_1..FDR.), y = cntrs, cex = 2.5, col = "goldenrod1",pch = 19)
lines(x = -log10(go$upload_1..FDR.), y = cntrs, col = "goldenrod1",lwd = 4)
axis(side = 3, font = 2, cex.axis = 2)
mtext(text = "-log10 (FDR)", side = 3, line = 3, font = 2, cex = 2)
dev.off()



