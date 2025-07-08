# Pipeline for processing next-generation sequencing (fastq) files from raw sequencing.

The Institute of Disease Modeling has been part of an interdisciplinary  collaboration with Elizabeth Thiele (Vassar) and James Cotton (Wellcome Sanger Institute) to maximize the value of epidemiological and genetics data to understand Guinea worm transmission in Chad. Preliminary analyses by IDM has shown whole mitochondrial genome data can give higher resolution information about genetic relatedness in a population than the current three-locus method. 

This pipeline is can process NGS sequencing from whole mitochondrial DNA to variants, however default settings gave been updated to process samples that have been processed using an amplicon panel. 


## Config parameters

#### WGS versus amplicon defaults
The procedure for processing samples sequenced from standard whole genome sequencing compared to amplicons are minimal. Amplicon amplified samples do not undergo deduplication and produce additional quality reports to verify the success of primers.

#### Base quality score recalibration (BQSR) variants

In instances where a set of known variants is preferred to recalibrate sequencing errors, the file will check the existence of known variants provided in config. If a file is not provided, it will proceed to recalibrate using the provided samples in the metadata file with the iteration specified in the config. BQSR will take place with two iterations by default. The BQSR comparisons to the initial round of variant calling will guide if 2 iterations are sufficient to normalize read errors.

For Guinea worm, we provided a set of known variants from a set of 54 publicly available worm samples (Durrant et al. 2020, PLOS NTDs, ENA PRJEB1236) and some early amplicon batches. For more information on run specifications for these provided variants, see the bsqr_comparison.Rmd document in `./scripts.`  

This pipeline follows the steps in the above publication as best as possible. The GEM masking program was difficult to get running and mitochondrial mask regions were defined manually by the authors, thus some undesirable regions may be included in the SNPs. Instead of including a known variants file from masking, we followed GATK4 guidelines to call and filter variants. The guideline parameters are set as default parameters in the config file. 


## Set-up
Pipelines are wrapped into Conda virtual environments and organized by Snakemake. These must be run on a Linux based system.

1. Download and set up Miniconda (https://docs.conda.io/en/latest/miniconda.html)
2. Set up Mamba for faster package dependency resolution when creating environments.
```
conda install -n base -c conda-forge mamba
```
3. Run install to create a virtual environment with the necessary programs.
```
mamba env create -f ngs_align_updated.yaml
```

## Configure the pipeline
All settings for the pipeline live in the "configfile," an example is provided at `configGW_mtAll.yaml`. Edit this file to change the output directory and other settings. All sequencing file inputs depend on the metadata file specified; look at `metadata_example.tsv` for a template. 

#### Create metadata files (Optional)
A customizable script has been included to generate a metadata file from a directory path and sample key file. The script handles naming schemes for the Cornell and Qiagen sequencing centers. The sample file for Qiagen must contain a minimum of `sample_number` and `sample` comma separated columns for matching. If no sample file is provided, the script will rename samples `[s1...sN]` up to the number of samples in the directory. 

Note: Relevant sample information is in the Cornell name, sample key file is not required. 

```
mamba activate ngs_align_updated
$GIT='/path/to/git/clone/ngs_processing'
snakemake -s $GIT/scripts/generate_manifest.smk --config fastq_dir='/path/to/raw_reads_dir' output_file='/path/to/metadata.tsv' sample_key='/path/to/sample_key.csv'
```
If not providing a sample key, leave the option in the command line or provide a dummy file. If the sample key file does not exist, it will default to the `[s1...sN]` naming scheme. This naming scheme is beneficial to run a single round, but is not recommended as it will merge different, unrelated samples together in the joint calling round. 


## Run the pipeline 

### Single batch processing with known variants 

You can now run all steps of the pipeline with a single command. This command will run everything from quality control to variant calling. Change the number of cores and jobs here to fit your machine. 

**Note: For updates that allow for joint variant calling across batches as a single step, this pipleline assumes all processed batches are placed in a single main directory without additional nesting,** for example `/home/{.?}/guinea_worm/{batch_name}`. The cumulative batch genotype calls will be available in `/home/{.?}/guinea_worm/joint_calls_{known/discovery}`.

The script currently assumes that `"batch_{.?}"` is the structure of individual batch names. The output files will be nested in an `out/` directory (will not duplicate this pattern if specified in the config output path). The parent directory for the project is specified as the directory structure before `batch_{.?}`. If the user wishes to to add further nesting to output file, the update must be made in `wgs_functions.smk` on line 21.  

```
mamba activate ngs_align_updated
$GIT='/path/to/git/clone/ngs_processing'
snakemake -s $GIT/gw_processing.snakefile --configfile path/to/project/yaml/config_example.yaml --cores 8 --jobs 8 -R check_samples
```

Depending on the BQSR options above, this pipeline takes in `.g.vcf.gz` files from each sample in each batch, and combines them into a GATK GenomicsDB Variant filter parameters match hard filtering suggestions from GATK while also removing variants that had mostly missing data (controlled by the "max_missing parameter" in the config). The final unfiltered and filtered calls, with all previous batches, can be found in `joint_calls_known/vcf_files/{batch_name}.vcf.gz`.   


### Multibatch processing for variant discovery

Since Guinea worm is not an organism with a known set of high-quality variants, it is recommended to run as many specimens as possible in variant discovery. Specifically, it may improve the inclusion of lower frequency alleles since the quality of specific variants across the mitochondrial genome can vary between samples and may cause variants to be filtered from a smaller set of samples. For this project, we ran a set of high-quality variants discovered from a handful of published data and pilot amplicon methodology batch data for BSQR of individual batches with the idea to run iterative discovery when the historical backlog of samples was complete. 

The input for multibatch variant discovery is similar to a single batch, with a few minor updates:

1. The sample mapping files, concatenated across batches. The script automatically looks for sample alignment file existence based on the expected directory structure mentioned above. Technical sample replicates, by exact name, will be idnetified and merged into the /home/{recalibration_batch_name}/
2. Config file should not contain a path for `known_variants`. Leave as a blank list " `''`" to trigger the discovery recalibration version of the pipeline.

Iterative discovery calibration with the existing batches starts with the original alignment files. If files do not exist, users will need to run those specific samples as single batches through the alignment step only. To do this, add the `–until merge_add_groups` tag to the Snakemake command for a single batch run with missing samples. 

Some files are created directly in the single batch directories:
- `/home/{.?}/guinea_worm/{original_batch_name}/out/02_align/recalibrate/bsqr_{iteration}/{sample}_{iteration}Iter.bam`
- `/home/{.?}/guinea_worm/{original_batch_name}/out/02_align/recalibrate/summary/{sample}_{iteration}_bqsrCovariates.pdf`
- `/home/{.?}/guinea_worm/{original_batch_name}/out/03_variant_calls/discovery/bsqr_{iteration}/{sample}.g.vcf.gz`

Joint call directory:
If sample duplicates exist, merged BAM files from individual batch alignments and subsequent steps are found in the same file structure as individual batches but convert `{original_batch_name}` to `{recalibration_batch_name}`.

Regardless of technical duplicates or not, the iterative base recalibration variant calls with all user specified samples are found in `/home/{.?}/guinea_worm/joint_calling_discovery/{recalibration_batch_name}/out`. For each iteration, the relevant variant files files are found in the subdirectories `bsqr_{interation}/jointGenotypeFiltered.vcf`. The final iteration sample GVCF list is automatically copied to `/home/{.?}/guinea_worm/joint_calling_known/` for input to run with future batches, that will now use the same high-quality variants from this iteration for base recalibration. 

**Note: When a discovery recalibration is called, old GVCF lists are archived so that the latest version of individual sample calls are used for joint calling with future batches. The old batch GVCF list files, determined by timestamp relative to the copied recalibration GVCF list are moved to `/home/{.?}/guinea_worm/joint_calling_known/preBSQR_{recalibration_batch_name}/archive`.**


### Pipeline version updates 

#### July 2025:
* Updated Conda environment and latest versions of packages that follow Anaconda’s 2024 licensing changes for open source use ([see licensing terms](https://legal.anaconda.com/policies/en/?name=terms-of-service#anaconda-terms-of-service) that require any organization with 200 or more employees or contractors to purchase a paid license to use Anaconda's software). Packages in the environment pull only from conda-forge and bioconda channels.
* Deprecation of the Conda environment used solely to create batch summary mtDNA reports in RMarkdown. The necessary packages are now in a unified Conda environment. 
* Integration of multiple batches into the discovery recalibration calling version of the pipeline, with a new dependency on existing alignment files for user specified samples. 
* Improved resolution of technical replicates across sequencing batches. Instead of including the original file and the new file that has been merged with samples with identical names, it now only keeps the merged version of the sample for variant calling. 
* Improved code organization by moving functions that apply to specific rules to the respective Snakemake files instead of in wgs_functions.smk.   


#### October 2021:
* Single step single batch processing and joint genotyping with samples in existing data base
* Detection and merging of duplicate biological samples by exact name matching across sequencing batches. 
* Addition of base recalibration with known variants or iterative calibration using user specified samples
* Optional deduplication of aligned reads 

#### August 2021:
* Addition of rule that creates metadata files. Will be updated for different sequencing centers as needed. 
