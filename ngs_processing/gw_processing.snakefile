'''
Processing next-generation sequencing files for variant calling
Author - Jessica Ribado
Date - 2019/12

This pipeline follows Johannes Köster's template (https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling) based on the GATK 4.0 Best Practices for calling variants from sequencing data.
(https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165)

This pipeline was modified from the outline above to incorporate variant calls for non-model organisms. It allows for the iteration of discovering high quality variants that are used to recalibrate variant calls across the genome. For organisms with known variants (ex. malaria), this pipeline is essentially equivalent to Johannes'.
'''

################################################################################
#  call python modules and functions
include: "scripts/wgs_functions.smk"
include: "scripts/reference_index.smk"
# include: "scripts/gem_mask.smk"
include: "scripts/qc.smk"
include: "scripts/alignment.smk"
include: "scripts/primer_check.smk"
include: "scripts/variant_calling_gvcf.smk"
include: "scripts/variant_filtering.smk"

################################################################################
rule all:
    input:
#        expand(join(PROJECT_DIR,  "01_processing/00_qc_reports/{step}_multiqc/multiqc_report.html"), step=['pre', 'post']),
        expand(join(PROJECT_DIR, "01_processing/02_align/{ena_id}_pairAligned.bam"), ena_id=ENA_ID),
        expand(join(PROJECT_DIR, "01_processing/02_align/coverage/{ena_id}_pairAligned_coverage.txt"), ena_id=ENA_ID),
        expand(join(PROJECT_DIR, "01_processing/02_align/primer_counts/{ena_id}_primerCounts_summary.txt"), ena_id=ENA_ID),
        expand(join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.bam"), sample=SAMPLE_PREFIX, iteration=config['iteration']),
        expand(join(PROJECT_DIR, "01_processing/02_align/recalibrate/{sample}_{iteration}Iter_bsqrCovariates.pdf"), sample=SAMPLE_PREFIX, iteration= list(range(1,config['iteration']+1))),
        expand(join(PROJECT_DIR, "01_processing/03_variant_calls/jointGenotype_{iteration}Iter_filtered.{extension}"),  iteration=list(range(1,config['iteration']+1)), extension=['vcf.gz', 'txt'])
        #join(PROJECT_DIR, "01_processing/03_variant_calls/bcftools/all_merged.txt")
