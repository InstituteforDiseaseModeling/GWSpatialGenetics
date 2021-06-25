# read in the clustering results, add a new column to the metadata
g.f <- snakemake@input[['genotype']]
metadata.f <- snakemake@input[['metadata']]
metadata.outf <- snakemake@output[['new_metadata']]
cuth <- snakemake@params[['hclust_height']]

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

# add to metadata file
metadata <- read.table(metadata.f, sep='\t', quote='', header=T)
metadata$new_cluster <- clusters[metadata$name]

# write new version
write.table(metadata, metadata.outf, sep='\t', quote=F, row.names = F, col.names = T)

