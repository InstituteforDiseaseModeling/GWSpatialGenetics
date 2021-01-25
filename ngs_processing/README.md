# Pipeline for processing next-generation sequencing (fastq) files from raw sequencing.

The Institute of Disease Modeling has been part of an interdisciplinary  collaboration with Elizabeth Thiele (Vassar) and James Cotton (Wellcome Sanger Institute) to maximize the value of epidemiological and genetics data to understand Guinea worm transmission in Chad. Preliminary analyses by IDM has shown whole mitochondrial genome data can give higher resolution information about genetic relatedness in a population than the current three-locus method. 
This pipeline is intented to process NGS sequencing from whole mitochondrial DNA to variants, with additional feautres to count primers and coverage. 


# Download publicly available genomes:
I downloaded the data corresponding to https://www.biorxiv.org/content/10.1101/808923v1. The data are available from the European Nucleotide Archive project ERP117282. This analysis focused solely on the samples collected in Chad as specified in Supplementary Table 1.

I attempted to follow the pipeline from the preprint as best as possible. The GEM masking program was difficult to get running and mitochondrial mask regions were defined manually by the authors, thus some undesirable regions may be included in the SNPs. Instead of including a known variants file from masking, I followed GATK4 guidelines to call variants, filter HQ (QD > 30), and use those in the base quality recalibration step until convergence.

# Set-up
Pipelines are wrapped into Conda virtual environments and organized by Snakemake. These must be run on a Linux based system.

Step 1: Download and set up Miniconda (https://docs.conda.io/en/latest/miniconda.html)
Step 2: Run install to create a virtual environment with the necessary programs.
```
conda env create -f ngs_align.yml
```

# Configure the pipeline
All settings for the pipeline live in the "configfile," an example is provided at `configGW_mtAll.yaml`. Edit this file to change the output directory and other settings. All sequencing file inputs depend on the metadata file specified; look at `metadata_example.tsv` for a template. 

# Run the pipeline
You can now run all steps of the pipeline with a single command. This will run everything from quality control to variant calling. Change the number of cores and jobs here to fit your machine. 
```
snakemake -s path/to/git/clone/wg_processing.snakefile --configfile path/to/project/yaml/project.yaml --cores 8 --jobs 8
```
Base quality score recalibration (BQSR) will take place with two iterations by default. The BQSR comparisons to the initial round of variant calling will guide if 2 iterations are sufficient to normalize read errors.

That's it! At the end you will have VCF files with genotypes at positions. 
