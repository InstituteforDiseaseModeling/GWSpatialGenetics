library(viridis)
library(gplots)
library(rafalib)
# read in the clustering results, add a new column to the metadata
g.f <- snakemake@input[['genotype']]
metadata.f <- snakemake@input[['metadata']]
metadata.outf <- snakemake@output[['new_metadata']]
cuth <- snakemake@params[['hclust_height']]
outdir <- snakemake@params[['outdir_figures']]
# read genotype matrix
g <- read.table(g.f, sep='\t', quote='', header=T)
# discard extra columns
g <- t(g[,3:ncol(g)])
# missing set to NA
g[g=="."] <- NA
# convert to numeric matrix
g <- data.matrix(as.data.frame(g, stringsAsFactors = F)) -1
# bad samples should be removed already
# g <- g[!(rownames(g) %in% remove.samples),]

# do hclust and cut tree
# using manhattan distance here
h <- hclust(dist(g, method = 'manh'))
clusters <- cutree(h, h=cuth)

# Add some plotting functionality
# a hclust plot with the cutoff as a line
if (!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}
pdf(file.path(outdir, 'clustering_dendrogram.pdf'), height=5, width=15, pointsize = 8)
mypar()
plot(h)
abline(h=cuth, col='firebrick', lty=2, lwd=2)
dev.off()

pdf(file.path(outdir, 'variant_heatmap.pdf'), height=12, width=12, pointsize = 8)
heatmap.2(t(g), Rowv = NA, trace='none',
          col=viridis(32), na.color = 'firebrick',
          hclustfun = function(x) h,
          key.xlab = "0=ref, 1=alt, red=MISSING",
          key.title = "Variant calls",
          lhei = c(1,6),
          margins=c(4,4))

dev.off()

# add clusters to metadata file
metadata <- read.table(metadata.f, sep='\t', quote='', header=T)
metadata$new_cluster <- clusters[metadata$name]

# write new version
write.table(metadata, metadata.outf, sep='\t', quote=F, row.names = F, col.names = T)
