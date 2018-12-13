output <- commandArgs(T)[1]

pdf(output, wid=7, hei=7)
vn <- venn(x, show.plot = F)
y <- sapply(attr(vn, 'intersections'), length)
names(y) <- gsub(':', '&', names(y))
  
col_set=brewer.pal(4, 'Set3')
print(plot(euler(y), quantities=T, legend=T, fills=list(fill=col_set, alpha=1)))
dev.off()
