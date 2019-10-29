argv <- commandArgs(T)
gene_anno_file <- argv[1]
ganno <- read.table(gene_anno_file, sep='\t', header = T, quote = '', row.names = 1, comment.char = '')

## hyper
geneId <- read.table('anno/hyper_geneId')$V1
hyper_files <- dir(path = 'anno', pattern = 'hyper_list_anno.xls', full.names = T)
hyper_dfm <- as.data.frame(matrix(0, nrow = length(geneId), ncol = length(hyper_files)))
rownames(hyper_dfm) <- geneId
colnames(hyper_dfm) <- sub('anno/', '', sub('_hyper_list_anno.xls', '', hyper_files))

for (i in 1:ncol(hyper_dfm)) {
  anno <- read.table(hyper_files[i], header = T, sep = '\t', quote = '')
  anno <- anno[anno$annotation != 'Distal Intergenic', ]
  hyper_dfm[rownames(hyper_dfm) %in% anno$geneId, i] <- hyper_dfm[rownames(hyper_dfm) %in% anno$geneId, i] + 1
}

hyper_dfm <- merge(hyper_dfm, ganno, by.x=0, by.y=0, all.x=T)
colnames(hyper_dfm)[1] <- 'geneId'
hyper_dfm <- hyper_dfm[order(rowSums(hyper_dfm[2:length(hyper_files)]), decreasing = T), ]
write.table(hyper_dfm, 'anno/hyper_DMR_target_genes_1yes_0no.xls', sep='\t', quote=F, row.names=F)

## hypo
geneId <- read.table('anno/hypo_geneId')$V1
hypo_files <- dir(path = 'anno', pattern = 'hypo_list_anno.xls', full.names = T)
hypo_dfm <- as.data.frame(matrix(0, nrow = length(geneId), ncol = length(hypo_files)))
rownames(hypo_dfm) <- geneId
colnames(hypo_dfm) <- sub('anno/', '', sub('_hypo_list_anno.xls', '', hypo_files))

for (i in 1:ncol(hypo_dfm)) {
  anno <- read.table(hypo_files[i], header = T, sep = '\t', quote = '')
  anno <- anno[anno$annotation != 'Distal Intergenic', ]
  hypo_dfm[rownames(hypo_dfm) %in% anno$geneId, i] <- hypo_dfm[rownames(hypo_dfm) %in% anno$geneId, i] + 1
}

hypo_dfm <- merge(hypo_dfm, ganno, by.x=0, by.y=0, all.x=T)
colnames(hypo_dfm)[1] <- 'geneId'
hypo_dfm <- hypo_dfm[order(rowSums(hypo_dfm[2:length(hypo_files)]), decreasing = T), ]
write.table(hypo_dfm, 'anno/hypo_DMR_target_genes_1yes_0no.xls', sep='\t', quote=F, row.names=F)
