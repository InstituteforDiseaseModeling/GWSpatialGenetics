# Pipeline for processing next-generation sequencing (fastq) files from raw sequencing.

The Institute of Disease Modeling has been part of an interdisciplinary  collaboration with Elizabeth Thiele (Vassar) and James Cotton (Wellcome Sanger Institute) to maximize the value of epidemiological and genetics data to understand Guinea worm transmission in Chad. Preliminary analyses by IDM has shown whole mitochondrial genome data can give higher resolution information about genetic relatedness in a population than the current three-locus method. 

This pipeline is can process NGS sequencing from whole mitochondrial DNA to variants, however default settings gave been updated to process samples that have been processed using an amplicon panel. 

**2024/08: EICC_Version branch contains files updated by Robert Arthur at the Emory Integrated Computational Core. They are not tested by IDM, for any issues with run errors, contact Robert Arthur directly at robert.arthur@emory.edu.**

#### Pipeline updates as of October 2021:
* Single step single batch processing and joint genotyping with samples in existing data base
* Detection and merging of duplicate biological samples by exact name matching across sequencing batches. 
* Addition of base recalibration with known variants or iterative ecalibration using user specified samples
* Optional deduplication of aligned reads 

#### Pipeline updates as of August 2021:
* Addition of rule that creates metadata files. Will be updated for different sequencing centers as needed. 


## Config parameters

#### WGS versus amplicon defaults
The procedure for processing samples sequenced from standard whole genome sequencing compared to amplicons are mininal. Amplicon amplified samples do not undergo deduplication and produce additional quality reports to verify the success of primers.

#### Base quality score recalibration (BQSR) variants
In instances where a set of known variants is preferred to recalibrate sequencing errors, the file will check the existance of known variants provided in config. If a file is not provided, it will proceed to recalibrate using the provided samples in the metadata file with the interation specified in the config. BQSR will take place with two iterations by default. The BQSR comparisons to the initial round of variant calling will guide if 2 iterations are sufficient to normalize read errors.

**Note: Cumulative genotype databases are generated separately for known or iterative calls from each batch.** 

For Guinea worm, we provided a set of known variants from a set of 54 publicly available worm samples (Durrant et al. 2020, PLOS NTDs, ENA PRJEB1236). For more information on run specifications for these provided variants, see the bsqr_comparison.Rmd document in `./scripts.`  

This pipeline follows the steps in the above publication as best as possible. The GEM masking program was difficult to get running and mitochondrial mask regions were defined manually by the authors, thus some undesirable regions may be included in the SNPs. Instead of including a known variants file from masking, we followed GATK4 guidelines to call and filter variants. The default parameters are 


## Set-up
Pipelines are wrapped into Conda virtual environments and organized by Snakemake. These must be run on a Linux based system.

1. Download and set up Miniconda (https://docs.conda.io/en/latest/miniconda.html)
2. Run install to create a virtual environment with the necessary programs.
```
conda env create -f ngs_align.yaml
```

## Configure the pipeline
All settings for the pipeline live in the "configfile," an example is provided at `configGW_mtAll.yaml`. Edit this file to change the output directory and other settings. All sequencing file inputs depend on the metadata file specified; look at `metadata_example.tsv` for a template. 

#### Create metadata files (optional)
A customizable script has been included to generate a metadata file from a directory path and sample key file. The script handles naming schemes for the Cornell and Qiagen sequencing centers. The sample file for Qiagen must contain a minimum of `sample_number` and `sample` comma separated columns for matching. If no sample file is provided, the script will rename samples `[s1...sN]` up to the number of samples in the directory. 

Note: Relevant sample information is in the Cornell name, sample key file is not required. 

```
snakemake -s /path/to/git/clone/scripts/generate_manifest.smk --config fastq_dir='/path/to/raw_reads_dir' output_file='/path/to/metadata.tsv' sample_key='/path/to/sample_key.csv'
```
If not providing a sample key, leave the option in the command line or provide a dummy file. If the sample key file does not exist, it will default to the `[s1...sN]` naming scheme. This naming scheme is beneficial to run a single round, but is not recommended as it will merge different, unrelated samples together in the joint calling round. 


## Run the pipeline 

You can now run all steps of the pipeline with a single command. This will run everything from quality control to variant calling. Change the number of cores and jobs here to fit your machine. The `--use-conda` flag has been added to source a custom conda environment for the final rule in the pipeline, creation of a primer QC report via R markdown. 

**Note: For updates that allow for joint calling across batches as a single step, this pipleline assumes all processed batches are placed in a single main directory without additional nesting,** for example `/home/{.?}/guinea_worm/{batch_name}`. The cumulative batch genotype calls will be available in `/home/{.?}/guinea_worm/joint_calls_{recal/known}`.

The script currently assumes that the batch name is the last directory provided in the config file. If the user wishes to to add further nesting to output file, the update must be made in `wgs_functions.smk` on line 21. The user must specify the split structure to capture above, where `FILE_STRUCTURE[0]` would return `/home/{.?}/guinea_worm/{batch_name}` and `FILE_STRUCTURE[1]` would return `{batch_name}`. 

```
snakemake -s path/to/git/clone/gw_processing.snakefile --configfile path/to/project/yaml/config_example.yaml --cores 8 --jobs 8 --use-conda -R check_samples
```

Depending on the BQSR options above, this pipeline takes in `.g.vcf.gz` files from each sample in each batch, and combines them into a GATK GenomicsDB Variant filter paramters match hard filtering suggestions from GATK while also removing variants that had mostly missing data (controlled by the "max_missing parameter" in the config). The final unfiltered and filtered calls, with all previous batches, can be found in `joint_calls_{recal/known}/vcf_files/{batch_name}.vcf.gz`.   

Individual batch calls are still run by default when using batch reaclibration, and the unfiltered and filtered VCFs can be found in `/home/{.?}/guinea_worm/{batch_name}/haplocall_{max_iteration}/joint_genotype[Filtered].vcf.gz`. For known variants, update the `wgs_functions.smk` user_output function to include the currently commented out file to run that rule.  

