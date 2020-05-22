wildcard_constraints:
    sample="[^/]+",
    iteration="[0-9]+",
    n="[0-9]+",

################################################################################
rule call_haplotypes:
    input:
        bam = recursive_bam,
        bam_idx = rules.merge_bam_index.output # not directly used in rule
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz")
    threads: 8
    params:
        ploidy = config['variant_calling']['ploidy']
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
rule combine_gvcf:
    input: expand(join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz"), sample=SAMPLE_PREFIX, iteration =config['iteration'])
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/combined.g.vcf")
    params:
        gvcf_string = expand("--variant " + join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz"), sample=SAMPLE_PREFIX, iteration =config['iteration'])
    shell: """
		gatk CombineGVCFs \
            --reference {REF_FILE} \
            {params.gvcf_string} \
            --output {output} \
    """

################################################################################
rule joint_genotyping:
    input: rules.combine_gvcf.output
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/joint_genotype.vcf.gz")
    shell: """
		gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant {input} \
            --output {output} \
    """

################################################################################
rule snp_discovery_vcf:
    input:  rules.joint_genotyping.output
    output: temp(join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/joint_genotypeFiltered_SNPs.vcf.gz"))
    shell: """
		gatk SelectVariants \
            --reference {REF_FILE} \
            --variant {input} \
            --select-type-to-include SNP \
            --output {output} \
    """

################################################################################
rule filter_discovery_vcf:
    input:  rules.snp_discovery_vcf.output
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/joint_genotypeFiltered.vcf.gz")
    params:
        qd_score = config['variant_filter']['quality_depth'],
        fisher   = config['variant_filter']['fisher_strand'],
        map_qual = config['variant_filter']['mapping_quality'],
        map_root = config['variant_filter']['mapping_rootsq'],
        map_rank = config['variant_filter']['mapping_rank']
    shell: """
		gatk VariantFiltration  \
            --reference {REF_FILE} \
            --variant {input} \
            --output {output} \
            --filter-expression \"QD < {params.qd_score} || FS > {params.fisher} || MQ < {params.map_qual} || MQRankSum < {params.map_root} || ReadPosRankSum < {params.map_rank}\" \
            --filter-name "gatk_hardFilt" \
    """

###############################################################################
rule calculate_base_recalibration:
    input:
        bam = recursive_bam,
        discovery_vcf = rules.filter_discovery_vcf.output
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.table")
    shell:
     """
        gatk BaseRecalibrator \
            --reference {REF_FILE} \
            --input {input.bam} \
            --output {output} \
            --known-sites {input.discovery_vcf}
    """


################################################################################
rule apply_base_recalibration:
    input:
        bam = recursive_bam,
        cal = rules.calculate_base_recalibration.output
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.bam")
    shell: """
        gatk ApplyBQSR \
            --reference {REF_FILE} \
            --input {input.bam} \
            --bqsr-recal-file {input.cal} \
            --output {output}
    """


################################################################################
rule check_base_recalibration:
    input:
        initial_recal = join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_0/{sample}_0Iter.table"),
        final_recal  = join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.table")
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/{sample}_{iteration}Iter_bsqrCovariates.pdf")
    shell: """
        gatk AnalyzeCovariates \
            -before {input.initial_recal} \
            -after {input.final_recal} \
            -plots {output}
    """


################################################################################
# gatk orthogonal approach
################################################################################
# skips indels
rule bcftool_variants:
    input:
        bam      = rules.merge_add_groups.output,
        bam_idx  = rules.merge_bam_index.output
    output:
        vcf = join(PROJECT_DIR, "01_processing/03_variant_calls/bcftools/individual_vcfs/{sample}.vcf.gz"),
        tbi = join(PROJECT_DIR, "01_processing/03_variant_calls/bcftools/individual_vcfs/{sample}.vcf.gz.tbi")
    params:
        ploidy = config['variant_calling']['ploidy'],
        min_qual = config['trim_galore']['quality']
    threads: 4
    shell: """
		 bcftools mpileup --threads {threads} -q {params.min_qual} \
            --fasta-ref {REF_FILE} {input.bam} | \
            bcftools call -m --ploidy {params.ploidy} -Ou -O z --output {output.vcf};
        tabix -p vcf {output.vcf}
    """
