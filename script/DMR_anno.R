suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ChIPseeker))

argv <- commandArgs(T)
input_bed <- argv[1]
sqlite <- argv[2]
gene_anno_file <- argv[3]
output_xls <- argv[4]
output_pie <- argv[5]

peak <- readPeakFile(input_bed)
txdb <- loadDb(sqlite)
ganno <- read.table(gene_anno_file, sep='\t', header = T, quote = '', row.names = 1)

gene_anno <- function(x) {
  peakAnno <- annotatePeak(x, TxDb=txdb, tssRegion = c(-2000,0))
  anno <- as.data.frame(peakAnno)
  anno <- merge(anno, ganno, by.x=31, by.y=0, all.x=T)
  return(list(peakAnno=peakAnno, anno=anno))
}

peak_anno <- gene_anno(peak) 
write.table(peak_anno$anno, output_xls, sep='\t', quote=F, row.names=F)

pdf(output_pie, width=7, height=5)
plotAnnoPie(peak_anno$peakAnno)
dev.off()
