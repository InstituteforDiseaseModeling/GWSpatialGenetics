max_iter = config['max_iter']
###############################################################################
rule variant_filter:
    input: rules.filter_discovery_vcf.output
    output: join(PROJECT_DIR, "03_variant_calls/jointGenotype_{iteration}Iter_filtered.vcf.gz")
    params:
        dp_min = config['variant_filter']['read_depth_min'],
        limits = config['variant_filter']['rank_sum']
    shell: """
        gatk VariantFiltration \
            --reference {REF_FILE} \
            --output {output} \
            --variant {input} \
            --filter-expression "ReadPosRankSum <= {params.limits} || ReadPosRankSum >= -{params.limits}" \
            --filter-name "ReadPosRankSum" \
            --filter-expression "ReadQRankSum <= {params.limits} || ReadQRankSum >= -{params.limits}" \
            --filter-name "ReadQRankSum" \
            --filter-expression "MQRankSum <= {params.limits} || MQRankSum >= -{params.limits}" \
            --filter-name "MQRankSum" \
            --filter-expression "ClippingRankSum <= {params.limits} || ClippingRankSum >= -{params.limits}" \
            --filter-name "ClippingRankSum" \
            --filter-expression "DP > {params.dp_min}" \
            --filter-name "DPMin"
    """

###############################################################################
rule vcf2txt:
    input: rules.variant_filter.output
    output: join(PROJECT_DIR, "03_variant_calls/jointGenotype_{iteration}Iter_filtered.txt")
        #join(PROJECT_DIR, "03_variant_calls/bcftools/all_merged.txt")
    shell: """
        gatk VariantsToTable \
            --variant {input} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -GF DP  \
            -raw \
            --output {output}
    """


###############################################################################
# remove samples with low depth as defined in the Rmarkdown processing
# and pick variants that have a call in at least a proportion of samples
rule vcf_max_missing: 
    input: 
        vcf =  join(PROJECT_DIR, f"03_variant_calls/jointGenotype_{max_iter}Iter_filtered.vcf.gz"),
        remove_samples = rules.post_qc_report.output
    output:
        filtered = join(PROJECT_DIR, f"04_variant_calls_final/jointGenotype_{max_iter}Iter_filtered.vcf"),
        decomposed = join(PROJECT_DIR, f"04_variant_calls_final/jointGenotype_{max_iter}Iter_filtered_decomposed.vcf")
    params:
        max_missing = config['variant_filter']['max_missing'],
        low_depth_file = rules.post_qc_report.params.output_low_depth
    shell: """  
        vcftools --gzvcf {input.vcf} --remove {input.remove_samples}  --recode -c | \
        vcftools --vcf - --max-missing {params.max_missing} --recode -c >  {output.filtered}
	vt decompose {output.filtered} > {output.decomposed}
    """
