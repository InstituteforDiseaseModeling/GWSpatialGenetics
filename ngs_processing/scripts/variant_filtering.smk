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
