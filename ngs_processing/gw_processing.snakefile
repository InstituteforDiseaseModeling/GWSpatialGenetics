'''
Processing next-generation sequencing files for variant calling
Author - Jessica Ribado
Date - 2019/12

This pipeline follows Johannes Köster's template 
(https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling)
 based on the GATK 4.0 Best Practices for calling variants from sequencing data.
(https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165)

This pipeline was modified from the outline above to incorporate variant calls for 
non-model organisms. It allows for the iteration of discovering high quality variants 
that are used to recalibrate variant calls across the genome. For organisms with known 
variants (ex. malaria), this pipeline is essentially equivalent to Johannes'.
'''

################################################################################
#  call python modules and functions
include: "scripts/wgs_functions.smk"
# include: "scripts/gem_mask.smk"
if isfile(VARIANT_FILE):
    include: "scripts/reference_index.smk"
    include: "scripts/qc.smk"
    include: "scripts/alignment.smk"
    include: "scripts/primer_check.smk"
    include: "scripts/variant_calling_gvcf_knownVariants.smk"
    include: "scripts/merge_batch_replicates.smk"
else:
    include: "scripts/variant_calling_gvcf_discovery.smk"
include: "scripts/joint_calling.smk"
include: "scripts/variant_filtering.smk"


################################################################################
rule all:
    input:
        user_outputs()
