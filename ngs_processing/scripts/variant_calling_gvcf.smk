wildcard_constraints:
    sample="[^/]+",
    # iteration="[1-9][0-9]*",
    iteration="[0-9]+",
    n="[0-9]+",

def get_bam(wildcards):
    n = int(wildcards.iteration)
    # print(wildcards)
    if n == 0:
        return join(PROJECT_DIR, "02_align/recalibrate/%s_%dIter.bam")  % (wildcards.sample, n)
        # sys.exit('shouldnt be called with zero!')
    if n > 0:
        return join(PROJECT_DIR, "02_align/recalibrate/haplocall_%d/%s_%dIter.bam") % (n-1, wildcards.sample, n-1)
    else:
        raise ValueError("Iteration steps must be an integer: received %s" % iteration)


# could do some kind of checkpoint here where we check the iterations of 
# base recalibration and decide on the final number to use....
# def find_best_iter(wildcards):
    # outputs = checkpoints.metabat.get(**wildcards).output[0]
    # return glob_wildcards(join(outputs, "{metabat_bin}.fa")).metabat_bin


################################################################################
rule call_haplotypes:
    input:
        bam = get_bam,
        bam_idx = rules.merge_bam_index.output, # not directly used in rule
        refdict = rules.create_ref_dict.output, # not directly used in rule
        faidx = rules.build_ref_faidx.output    # not directly used in rule
    output: join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz")
    threads: 2
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
    input: 
        lambda wildcards: expand(join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz"), sample=unique_samples, iteration=wildcards.iteration)
    output: 
        join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/combined.g.vcf")
    params:
        gvcf_string = lambda wildcards: expand("--variant " + join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz"), sample=unique_samples, iteration=wildcards.iteration)
    shell: """
        gatk CombineGVCFs \
            --reference {REF_FILE} \
            {params.gvcf_string} \
            --output {output} \
    """

################################################################################
rule joint_genotyping:
    input: rules.combine_gvcf.output
    output: join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/joint_genotype.vcf.gz")
    shell: """
        gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant {input} \
            --output {output} \
    """

################################################################################
rule snp_discovery_vcf:
    input:  rules.joint_genotyping.output
    output: temp(join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/joint_genotypeFiltered_SNPs.vcf.gz"))
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
    output: join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/joint_genotypeFiltered.vcf.gz")
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
rule calculate_apply_base_recalibration:
    input:
        bam = get_bam,
        discovery_vcf = join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/joint_genotypeFiltered.vcf.gz")
    output: 
        table = join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.table"),
        bam = join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.bam")
    shell: """
        gatk BaseRecalibrator \
            --reference {REF_FILE} \
            --input {input.bam} \
            --output {output.table} \
            --known-sites {input.discovery_vcf}
        gatk ApplyBQSR \
            --reference {REF_FILE} \
            --input {input.bam} \
            --bqsr-recal-file {output.table} \
            --output {output.bam}
    """



################################################################################
rule check_base_recalibration:
    input:
        initial_recal = join(PROJECT_DIR, "02_align/recalibrate/haplocall_0/{sample}_0Iter.table"),
        final_recal  = join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.table")
    output: join(PROJECT_DIR, "02_align/recalibrate/{sample}_{iteration}Iter_bsqrCovariates.pdf")
    shell: """
        gatk AnalyzeCovariates \
            -before {input.initial_recal} \
            -after {input.final_recal} \
            -plots {output}
    """

