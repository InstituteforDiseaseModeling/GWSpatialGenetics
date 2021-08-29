# Pipeline for processing next-generation sequencing (fastq) files from raw sequencing.

The Institute of Disease Modeling has been part of an interdisciplinary  collaboration with Elizabeth Thiele (Vassar) and James Cotton (Wellcome Sanger Institute) to maximize the value of epidemiological and genetics data to understand Guinea worm transmission in Chad. Preliminary analyses by IDM has shown whole mitochondrial genome data can give higher resolution information about genetic relatedness in a population than the current three-locus method. 
This pipeline is intended to process NGS sequencing from whole mitochondrial DNA to variants, with additional features to count primers and coverage. 


## Download publicly available genomes:
I downloaded the data corresponding to https://www.biorxiv.org/content/10.1101/808923v1. The data are available from the European Nucleotide Archive project ERP117282. This analysis focused solely on the samples collected in Chad as specified in Supplementary Table 1.

I attempted to follow the pipeline from the preprint as best as possible. The GEM masking program was difficult to get running and mitochondrial mask regions were defined manually by the authors, thus some undesirable regions may be included in the SNPs. Instead of including a known variants file from masking, I followed GATK4 guidelines to call variants, filter HQ (QD > 30), and use those in the base quality recalibration step until convergence.

## Set-up
Pipelines are wrapped into Conda virtual environments and organized by Snakemake. These must be run on a Linux based system.

Step 1: Download and set up Miniconda (https://docs.conda.io/en/latest/miniconda.html)
Step 2: Run install to create a virtual environment with the necessary programs.
```
conda env create -f ngs_align.yml
```

## Configure the pipeline
All settings for the pipeline live in the "configfile," an example is provided at `configGW_mtAll.yaml`. Edit this file to change the output directory and other settings. All sequencing file inputs depend on the metadata file specified; look at `metadata_example.tsv` for a template. 

##### Update 08/2021: Create metadata files (optional)
A customizable script has been included to generate a metadata file from a directory path and sample key file. The script handles naming schemes for the Cornell and Qiagen sequencing centers. The sample file for Qiagen must contain a minimum of `sample_number` and `sample` tab delimited columns for matching. If no sample file is provided, the script will rename samples [s1...sN] up to the number of samples in the directory. 

Note: Relevant sample information is in the Cornell name, sample key file is not required. 

```
snakemake -s /path/to/git/clone/scripts/generate_manifest.smk --config fastq_dir='/path/to/raw_reads_dir' output_file='/path/to/metadata.tsv' sample_key='/path/to/sample_key.txt'
```
If not providing a sample key, leave the option in the command line or provide a dummy file. If the sample key file does not exist, it will default to the [s1...sN] naming scheme.


## Run the pipeline
You can now run all steps of the pipeline with a single command. This will run everything from quality control to variant calling. Change the number of cores and jobs here to fit your machine. The `--use-conda` flag has been added to source a custom conda environment for the final rule in the pipeline, creation of a primer QC report via R markdown. 
```
snakemake -s path/to/git/clone/gw_processing.snakefile --configfile path/to/project/yaml/configGW_mtAll.yaml --cores 8 --jobs 8 --use-conda
```
Base quality score recalibration (BQSR) will take place with two iterations by default. The BQSR comparisons to the initial round of variant calling will guide if 2 iterations are sufficient to normalize read errors.

That's it! At the end you will have VCF files with genotypes at positions. The final, filtered callset will be in a folder `04_variant_calls_final` and will be filtered to remove samples with less than 50% breadth at a depth of 5 reads, and variants that had mostly missing data (controlled by the "max_missing parameter" in the config)

## Combine multiple batches and join variant calling
If you have multiple runs and want to combine all samples to enable better joint variant calling, use `gw_joint_calling.snakefile` and `config_joint_calling.yaml`. This pipeline takes in `.g.vcf.gz` files from each sample in each batch, and combines them into a GATK GenomicsDB. Variant calling is then done on the whole set. The configfile needs a file with a list of g.vcf files, one per line, an output directory, and the reference fasta file. Then, call the pipeline like so:
```
snakemake -s path/to/git/clone/gw_joint_calling.snakefile --configfile path/to/project/yaml/config_joint_calling.yaml --cores 8 --jobs 8
```

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
