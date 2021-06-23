# read in the clustering results, add a new column to the metadata
cdhit.f <- snakemake@input[['clusters']]
metadata.f <- snakemake@input[['metadata']]
metadata.outf <- snakemake@output[['new_metadata']]

# read the cluster table
cluster.table <- read.table(cdhit.f, sep=">", header=F, fill=T)

# it's in a god-awful format, so parse it in a loop
cluster.list <- list()
i <- 1
while (i <= nrow(cluster.table)){
    # start of a new cluster
    if (cluster.table[i,1]==""){
        cluster.name <- strsplit(cluster.table[i,2], split=' ')[[1]][2]
        cluster.list[[cluster.name]] <- c()
    } else {
        sample.name <- strsplit(cluster.table[i,2], split="\\.")[[1]][1]
        cluster.list[[cluster.name]] <- c(cluster.list[[cluster.name]], sample.name)
    }
    i <- i+1
}


# add to metadata file
metadata <- read.table(metadata.f, sep='\t', quote='', header=T)
rownames(metadata) <- metadata$name
metadata$cdhit_cluster <- NA
for (c in names(cluster.list)){
    for (s in cluster.list[[c]]){
        metadata[s, 'cdhit_cluster'] <- c
    }
}
# write new version
write.table(metadata, metadata.outf, sep='\t', quote=F, row.names = F, col.names = T)

