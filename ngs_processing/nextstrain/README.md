## Interactive visualization with Nextstrain
[Nextstrain](https://nextstrain.org/) is an open-source visualization tool for pathogen genome data. We've built a pipeline to transform the variant call data into a format ready for visualization with Nextstrain. To use this, first install the program and dependencies in a new conda environment from the provided yaml file:
```
conda env create -f nextstrain_env.yml

# check the install worked
conda activate nextstrain
nextstrain check-setup --set-default
# the output from the above command should end with 
Setting default environment to native.
```
Then, copy the configuration file `nextstrain/config_nextstrain.yaml` to your working directory and change the parameters. You can control the method used for tree building and the output directory. 

Launch the snakemake workflow: 
```
# activate conda environment if necessary
conda activate nextstrain
snakemake -s nextstrain/nextstrain.snakefile --configfile nexstrain/config_nextstrain.yaml -j 1
```

Then, fire up the `auspice` viewer tool with the auspice path in your output directory:
```
auspice view --datasetDir nextstrain_iqtree/auspice/
# You may see an error message about "XX is not a valid address."
# I was able to fix this by appending the following to the command:
HOST="localhost"; auspice view --datasetDir nextstrain_iqtree/auspice/
```
Direct your browser to the link provided. The dataset is named 'GW' by default.

### Alternative clustering in Nexstrain
As an alternative to the tree-building options in Nexstrain, I built a simple clustering method that works directly on the variant call data. This method first decomposes bilallelic sites into respective single-allele variants, then does hierarchical clustering with the ward.D2 method on the manhattan distance between genotype vectors. The output of the clustering will be in the output directory "clustering_figures". You can set a height on the hierarchical clustering to cut the tree and produce a number of discrete clusters, which can be visualized on the nextstrain tree by selecting the "Color By" parameter as "new_cluster" in the web browser. 

If you want to try different clustering cutoffs, simply delete "outdir/clustering_figures/clustering_dendrogram.pdf", change the "hclust_height" parameter in the config, and re-run the Snakemake command. Only the last few rules will need to run, and it will be very quick.
