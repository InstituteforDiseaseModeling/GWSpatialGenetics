# mtDNA Analyses Files

This folder contains scripts for mtDNA analyses and visualizations. The libraries required to run the analysis file is not included in the conda environments for processing (Snakemake pipeline) or NextStrain visualization, however the Nextstrain environment does have a majority of the required packages. All necessary packages should download and load automatically however user experience may vary. We recommend for users to run the pipeline manually to confirm the necessary packages download and load correctly before using the command line version of the script. 

Outputs from mtDNA_analyses.R can be input into the NextStrain visualization tool, found in the parent directory.

Necessary files:
- Diploid VCF: Diploid VCF generated from the Guinea worm NGS tiling amplicon panel processing pipeline.
- Metadata: Metadata for the sequenced samples in the VCF file. Minimum fields include sampleID (Vassar generated), country, year, and host for standard analyses.

See the other optional command line options for filtering with 

```
mtDNA_analyses.R --help
```


## Output Files 

### Text files for flexible analyses
**GenomicSamplesMetadata_Database_{batch_name}.tsv**

Assumed the suffix used by the Genomics Working Group for metadata base name, but any metadata name file will receive a {batch_name} appendix. The script runs minor formatting checks on a few columns including host, GPS coordinates, and emergence date. The following additional columns are added.

* filt_site_missing_prop: The proportion of variant positions within a barcode that are missing.  
* excluded_for_analysis: Whether the sample is eligible for inclusion in analyses based on the proportion of missing barcode positions. Default threshold of position missingness set to 10% of the barcode.    
* successfully_sequenced: Confirmation of whether a sample did generate a sequencing file that was not filtered during earlier processing steps.    
* sequence: Barcode for each sample created from high quality variant positions.         
* group: Random numerical assignment for unique barcodes in order of appearance in the VCF.    
* frequency: The count for each group in the total dataset of available sequences.      
* amplicon_barcode: Same as group, but ordered by the frequency, i.e. the most frequent barcodes is now 1.         
* amplicon: Same as amplicon barcode for common barcodes, but groupings for less commonly observed barcodes for visualization ease. Generally, groupings for less frequent barcodes are observed <10, <5, or once in the total dataset. 


**{batch_name}_jointHaploidFilter.vcf.gz**

VCF file with high quality variants and samples created from user defined parameters. The reduction in samples and variants from the parameters are summarized in **{batch_name}_filterSummary.tsv**. 

**{batch_name}_filterSummary.tsv**

Summary of **{batch_name}_jointHaploidFilter.vcf.gz** properties, including the  number of samples, variants, and filtering parameters. 

**{batch_name}_relatedness.txt**

Long format of all pairwise comparisons of sequenced samples that have passed QC. Note, to reduce file size only the lower matrix of pairwise comparisons is saved. For the total length of the barcode for the pairwise relatedness denominator, refer to **{batch_name}_filterSummary.txt**. This is the predominant file used for analyses. 

- Var1: Vassar assigned labID of one sample in the pairwise comparison 
- Var2: Vassar assigned labID of one sample in the pairwise comparison 
- All: The dissimilarity of two sequences with all positions from the barcode created of high quality variants. 
- rmNA: The dissimilarity of two sequences, excluding the missing positions for either sequence in the pair from the barcode created of high quality variants.  
- missing: The proportion of missing positions within a pair of sequences for interpretation of the rmNA parameter. 
- (Optional) meters: The haversine distance between two samples, as provided in the metadata file. 

**{batch_name}_CountryBarcodeOverlapCounts.tsv**

A list of amplicon barcodes that are shared between countries within a year for further investigation of putative links between cases. Includes the number of specimen samples from each country that share that amplicon. 

**{batch_name}_CountryBarcodeOverlapSpecimens.tsv**

Long format of **{batch_name}_CountryBarcodeOverlapCounts.tsv** specifically including which samples by Vassar LabID are included for each amplicon barcode by year. The kinship_consideration_group is a dummy grouping variable for potentially automating kinship analsyes from microsatellite data.  


### Plots 

All plots from the script are automatically assigned a date (Y/M/D) of creation as a prefix. Wildcards of country below can also be replaced by specific groupings, such as user defined foci, villages, or epidemiological relevant cases. Flexibility to generate plots with more distinct groupings for the types of plots created is easily achieved since the plotting code is flexible to specify samples of interest. 

- barcode{country}.png: A barplot of barcodes found in a country over time, from the 'amplicon' column of the updated metadata. Barcode colors are consistent throughout country plots, and may require updating for instances where 80+ unique barcodes after groupings may exist.

- barcodeProportion{country}.png: A two panel barplot of barcodes found in each country, but scaled by proportion per year. Accompanied by a plot of the number of specimens sequences per year colored by the infection host. 

- pariwiseHexByDist_{country}.png: A scatter plot of genetic similarity and distance to broadly observe the relationship between there two variables. To limit plotting time and output file size, the frequency of points is reduced into a heatmap of hexagons. 

- simDist_{Sountry}.png: Similar to pairwiseHexByDist, but shows pairs as points as colors by same or +/- 1 year. Accompanied by plot panels with sample counts that successfully sequenced and counts by infection host. 

- relatednessCumulativeDensity.png: A cumulative density plot of pairwise relatedness within countries and between all countries. Country inclusion varies by specimen counts. 

- relatednessDensity.png: A smoothed density plot of pairwise relatedness within countries and between all countries. Country inclusion varies by specimen counts. 

- relatednessMatrix_{country/group}.png: Pairwise related matrices (upper matrix) and variant position availability for confidence (lower matrix). Accompanied by a plot of the number of specimens sequences per year colored by the infection host and the barcode grouping that would match the accompanying Nextstrain visualization. 

![Matrix](AmpliconProtocolSlide.png "General advice for interpretation of the genetic similarity and variant position missingness.")


## Versions
07/2024: First version of documentation for output files added. 

04/2024: Added VCF processing that was previously part of the NextStrain instance to the analysis file. This includes creating a haploid genotype from the diploid VCF (to catch heteroplasmy or nuclear gene duplicates), filtering variant positions based on sample quality, and filtering samples with lower sequence quality at variant positions. It also organizes several iterations of analysis code and visualization into a single script. Metadata checks are only subject to previously known inconsistencies with data (ie multiple formats for GPS), more may be needed in the future. 