## mtDNA Analyses Files

This folder contains scripts for mtDNA analyses and visualizations. The libraries required to run the analysis file is not included in the conda environments for processing (Snakemake pipeline) or NextStrain visualization, however the Nextstrain environment does have a majority of the required packages. All necessary packages should download and load automatically however user experience may vary. We recommend for users to run the pipeline manually to confirm the necessary packages download and load correctly before using the command line version of the script. 

Outputs from mtDNA_analyses.R can be input into the NextStrain visualization tool, found in the parent directory.

Necessary files:
- Diploid VCF: Diploid VCF generated from the Guinea worm NGS tiling amplicon panel processing pipeline.
- Metadata: Metadata for the sequenced samples in the VCF file. Minimum fields include sampleID (Vassar generated), country, year, and host for standard analyses.

See the other optional command line options for filtering with 

```
mtDNA_analyses.R --help
```

### Versions

04/2024: Added VCF processing that was previously part of the NextStrain instance to the analysis file. This includes creating a haploid genotype from the diploid VCF (to catch heteroplasmy or nuclear gene duplicates), filtering variant positions based on sample quality, and filtering samples with lower sequence quality at variant positions. It also organizes several iterations of analysis code and visualization into a single script. Metadata checks are only subject to previously known inconsistencies with data (ie multiple formats for GPS), more may be needed in the future. 