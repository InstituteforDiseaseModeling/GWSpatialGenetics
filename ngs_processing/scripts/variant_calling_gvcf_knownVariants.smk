import os, re
import numpy as np
import pandas as pd

# Rules that mirror interative recalibration calling, with different output files
###############################################################################
rule calculate_apply_known_base_recalibration:
    input:
        bam = rules.merge_add_groups.output.bam,
        known_vcf = config['variant_calling']['known_variants']
    output: 
        table = join(PROJECT_DIR, "02_align/recalibrate/known_variants/{sample}.table"),
        bam = join(PROJECT_DIR, "02_align/recalibrate/known_variants/{sample}.bam")
    shell: """
        gatk BaseRecalibrator \
            --reference {REF_FILE} \
            --input {input.bam} \
            --output {output.table} \
            --known-sites {input.known_vcf}
        gatk ApplyBQSR \
            --reference {REF_FILE} \
            --input {input.bam} \
            --bqsr-recal-file {output.table} \
            --output {output.bam}
    """

###############################################################################
rule check_base_known_recalibration:
    input: rules.calculate_apply_known_base_recalibration.output.table
    output: join(PROJECT_DIR, "02_align/recalibrate/{sample}_known_bqsrCovariates.pdf")
    shell: """
        gatk AnalyzeCovariates \
            -bqsr {input} \
            -plots {output}
    """    

################################################################################
rule call_haplotypes_known:
    input:
        bam = rules.calculate_apply_known_base_recalibration.output.bam,
        bam_idx = rules.merge_add_groups.output.bai, # not directly used in rule
        refdict = rules.create_ref_dict.output, # not directly used in rule
        faidx = rules.build_ref_faidx.output    # not directly used in rule
    output: join(PROJECT_DIR, "03_variant_calls", "known_variants", "{sample}.g.vcf.gz")
    threads: 2
    params:
        # ploidy = config['variant_calling']['ploidy']
        ploidy = 2
    shell: """
        gatk HaplotypeCaller \
            --reference {REF_FILE} \
            --input {input.bam} \
            --output {output} \
            -ploidy {params.ploidy} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads}
    """

################################################################################
rule gvcf_known_file:
    input: expand(join(PROJECT_DIR, "03_variant_calls", "known_variants", "{sample}.g.vcf.gz"), sample = list(set(metadata['sample'])))
    output: join(PROJECT_DIR, "03_variant_calls", "known_variants", "genomicsDBimport_gvcf_list.txt")
    params:     
        gvcf_dir = join(PROJECT_DIR, "03_variant_calls", "known_variants")
    run:
        gvcf_list_df(str(params.gvcf_dir), str(output))  


################################################################################
# additional rules for optional single batch joint genotype calling
################################################################################
rule batch_joint_genotyping:
    input: rules.GenomicsDBImport.output
    output: join(PROJECT_DIR, "03_variant_calls", "known_variants", "joint_genotype.vcf.gz")
    params: 
        db_dir = join(PROJECT_DIR, "03_variant_calls", "known_variants", "genomicsDB"),
    shell: """
        gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant gendb://{params.db_dir} \
            --output {output} \
    """

################################################################################
rule filter_batch_vcf:
    input:  rules.batch_joint_genotyping.output
    output: join(PROJECT_DIR, "03_variant_calls", "known_variants", "joint_genotypeAnnotated.vcf.gz")    
    params:
        qd_score = config['variant_filter']['quality_depth'],
        fisher   = config['variant_filter']['fisher_strand'],
        map_qual = config['variant_filter']['mapping_quality'],
        map_root = config['variant_filter']['mapping_rootsq'],
        map_rank = config['variant_filter']['mapping_rank'],
        depth    = config['variant_filter']['read_depth_min']
    shell: """
        gatk VariantFiltration  \
            --reference {REF_FILE} \
            --variant {input} \
           --filter-expression \"QD < {params.qd_score}\" --filter-name \"QD\" \
           --filter-expression \"FS > {params.fisher}\" --filter-name \"FS\" \
           --filter-expression \"MQ < {params.map_qual}\" --filter-name \"MQ\" \
           --filter-expression \"MQRankSum < {params.map_root}\" --filter-name \"MQRank\" \
           --filter-expression \"ReadPosRankSum < {params.map_rank}\" --filter-name \"ReadPosRank\"  \
           --filter-expression \"DP < {params.depth}\" --filter-name \"depth\"  \
           --genotype-filter-expression "isHet == 1" \
           --genotype-filter-name "isHetFilter" \
           --output {output}  
    """

################################################################################
rule snp_batch_vcf:
    input:  rules.filter_batch_vcf.output
    output: join(PROJECT_DIR, "03_variant_calls", "known_variants", "joint_genotypeFiltered.vcf.gz")
    params: 
        max_missing = 1 - config['variant_filter']['max_missing'] 
    shell: """
        gatk SelectVariants \
            --reference {REF_FILE} \
            --variant {input} \
            --select-type-to-include SNP \
            --exclude-filtered \
            --max-nocall-fraction {params.max_missing} \
            --output {output} 
    """

