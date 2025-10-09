# *Dracunculus medinensis* mitochondrial amplicon sequence processing 

Author: Jessica Ribado (jessica.ribado@gatesfoundation.org)

Date: 05/2025

## Identifying *Dracunculus medinensis* mitochondrial genome variants 

Single nucleotide polymorphism (SNP) variants were identified from the Dracunculus medinensis mitochondrial genome following the GATK v4.1.4 Best Practices pipeline for non-model organisms, with minor modifications as previously described [Van der Auwera et al 2013, Ribado et al 2021]. Briefly, bases with a mapping quality score below 20 were trimmed from the ends of short reads prior to alignment to the D. medinensis mitochondrial reference genome (version ENA JN555591.1) using BWA v0.7.17 [Li et al 2009]. Variant calling was performed on all aligned reads using the HaplotypeCaller tool with a ploidy setting of 2 to enable detection of heterozygous sites indicative of potential heteroplasmy or contamination. To improve variant call accuracy, a high-confidence set of SNPs was generated through a bootstrapping approach and used to recalibrate base quality scores due to variable amplicon sequencing depths. SNPs were filtered using the VariantFiltration tool with the following parameters: QD < 2.0, FS > 60.0, MQ < 40.0, MQRankSum < −12.5, ReadPosRankSum < −8.0, and DP < 10. Only SNPs were retained using SelectVariants with the --select-type-to-include SNP flag. Final base quality score recalibration was conducted using the BaseRecalibrator and ApplyBQSR tools iteratively (2 rounds) until convergence of variant calls was confirmed with visual inspection of quality plots.

<del>

**High quality variant pipeline prior to July 2025** A set of high-quality variants was first identified from an initial batch of 1178 samples representing *D. medinensis* specimens collected from dog cases or human infections from 10 countries and 12 years, but enriched for dog infections from Chad and human infections from South Sudan (Table S1). Some specimens were duplicated across sequencing batches for batch comparison. These samples averaged 71,584 (33-390-122,660) uniquely aligned reads and were used to establish the representative, high-confidence **279** variants for recalibration. For each subsequent sequencing batch, joint genotyping was performed using GATK’s GenotypeGVCFs tool. A consolidated list of all current and prior samples was generated using the GenomicsDBImport tool to ensure cumulative joint calling across batches.

Table S1: Geographical and spatial representation of *D. medinensis* specimen samples included in the identification on high-quality single nucleotide variants for GATK base recalibration. 

|      | 2006 | 2008 | 2009 | 2011 | 2012 | 2013 | 2014 | 2015 | 2016 | 2017 | 2018 | 2019 |
|------|------|------|------|------|------|------|------|------|------|------|------|------|
| Angola  | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 1    | 5    |
| Burkina Faso  | 2    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    |
| Cameroon  | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 1    |
| Côte d'Ivoire  | 2    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    |
| Chad  | 0    | 0    | 0    | 0    | 4    | 7    | 0    | 7    | 132  | 701  | 30   | 38   |
| Ethiopia  | 0    | 0    | 0    | 0    | 0    | 7    | 6    | 15   | 3    | 8    | 29   | 4    |
| Mali  | 0    | 1    | 0    | 1    | 0    | 7    | 22   | 2    | 0    | 0    | 16   | 3    |
| Niger  | 0    | 1    | 1    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    |
| South Sudan  | 0    | 0    | 0    | 0    | 0    | 14   | 87   | 4    | 2    | 0    | 14   | 3    |
| Sudan  | 0    | 0    | 0    | 0    | 0    | 1    | 0    | 0    | 0    | 0    | 0    | 0    |

</del>


**High quality variant pipeline after July 2025** High quality variants were identified using the complete library of Guinea worm specimens with sufficient quality DNA (17706 specimens, including 14 lab infections from ferrets and 6 historical specimens from human cases of unknown country origin). All samples were jointly genotyped using the GenotypeGVCFs tool in each recalibration round. The final recalibration resulted in 1,908 high-quality variants representing the diversity of specimens from 11 countries across 17 years for base recalibration. 

Table S1: Geographical and spatial representation of *D. medinensis* specimen samples included in the identification on high-quality single nucleotide variants for GATK base recalibration. 

|      | 2006 | 2008 | 2009 | 2011 | 2012 | 2013 | 2014 | 2015 | 2016 | 2017 | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 |
|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| Angola  | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 1    | 4    | 1    | 0    | 7    | 87   | 37   | 
| Burkina Faso  | 2    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 
| Cameroon  | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 1    | 6    | 9    | 50   | 391  | 585  | 
| Central African Republic   | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 1    | 1    | 1    | 0    | 
| Côte d'Ivoire  | 2    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    |
| Chad  | 0    | 0    | 0    | 1    | 8    | 8    | 29   | 287  | 1926 | 1306 | 1955 | 4166 | 2937 | 1229 | 1061 | 694  | 51   | 
| Ethiopia  | 0    | 0    | 0    | 0    | 0    | 16   | 9    | 42   | 23   | 65   | 42   | 52   | 79   | 5    | 11   | 4    | 3    | 
| Mali  | 0    | 1    | 0    | 1    | 0    | 9    | 35   | 3    | 5    | 0    | 18   | 7    | 11   | 17   | 68   | 54   | 28   | 
| Niger  | 0    | 1    | 1    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 
| South Sudan  | 0    | 0    | 0    | 0    | 0    | 14   | 100  | 8    | 2    | 0    | 15   | 15   | 1    | 4    | 14   | 2    | 24   | 
| Sudan  | 0    | 0    | 0    | 0    | 0    | 1    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 
| TAMU Lab  | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 14   | 
| Unknown  | 1    | 0    | 1    | 1    | 3    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 0    | 


The final joint genotyping file was filtered using the same criteria applied to identify high-confidence SNPs as described above. Additional custom quality control metrics that assess sequencing coverage, coverage evenness, and amplicon primer performance were automatically generated and summarized in a per-batch PDF report. All quality control and processing steps were implemented within a reproducible Snakemake workflow which is publicly available at: https://github.com/InstituteforDiseaseModeling/GWSpatialGenetics/ngs_processing.


## Additional resources and notes 

* GATK Guidance for filtering variants on non-model organisms for the version of the pipeline: https://gatk.broadinstitute.org/hc/en-us/articles/360036346192--Tool-Documentation-Index. Guidance was updated in 2024 to be VSQR instead of BSQR + Variant filtration. Since there is a high threshold for keeping variants in the pipeline and in determining variants within the barcode I do not see this as something we need to update in the processing pipeline, especially since The Broad is exploring new tools to replace VSQR.

* GATK BSQR on non-model organisms: https://evodify.com/gatk-in-non-model-organism/

* GATK variant filtering criteria: 
    - Definitions: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    - For non-model organisms: https://gatk.broadinstitute.org/hc/en-us/articles/360035532412-Can-t-use-VQSR-on-non-model-organism-or-small-dataset

* To get sample names from the known variant set - individual variant information per sample was removed to reduce file size 

    ```zcat /home/jribado/git/GWSpatialGenetics/ngs_processing/input_files/gw_known_amplicon.vcf.gz | grep ^#CHROM | sed 's/\t\t*/\n/g' | sed '1,9d' > bsqr_hq_samples.txt```

* To count of high-quality variants included in the BSQR set

    ```zcat GWSpatialGenetics/ngs_processing/input_files/gw_known_amplicon.vcf.gz | grep ^ENA | wc -l```


## References
Van der Auwera, G.A., Carneiro, M.O., Hartl, C., et al. (2013). From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Current Protocols in Bioinformatics, 43(1), 11.10.1–11.10.33. https://doi.org/10.1002/0471250953.bi1110s43

Ribado, J. V., Li, N. J., Thiele, E., Lyons, H., Cotton, J. A., Weiss, A., Ouakou, P. T., Moundai, T., Zirimwabagabo, H., Guagliardo, S. A. J., Chabot-Couture, G., & Proctor, J. L. (2021). Linked surveillance and genetic data uncovers programmatically relevant geographic scale of Guinea worm transmission in Chad. PLOS Neglected Tropical Diseases, 15(5), e0009609. https://doi.org/10.1371/journal.pntd.0009609

Li, H. and Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324


