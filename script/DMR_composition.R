library(RColorBrewer)
colset <- brewer.pal(3, 'Set2')

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

dfm <- read.table(input)
cmp <- dfm$V2/sum(dfm$V2)
names(cmp) <- dfm$V1

pdf(output, wid=6, hei=6)
pie(cmp, col = colset, border = 'white', clockwise = T, labels = paste0(round(cmp*100, 0), '%'), cex.lab = 3)
legend('topright', names(cmp), fill = colset, bty = 'n', border = NA)
dev.off()
