################################################################################
rule annotate_joint_vcf:
    input:  rules.joint_genotyping_all.output
    output: temp(join(KNOWN_DIR, "vcf_files",  f"{BATCH_NAME}_jointGenotypeAnnotated.vcf.gz"))    
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
           --genotype-filter-expression \"DP < {params.depth}\" \
           --genotype-filter-name "depth" \
           --genotype-filter-expression \"isHet == 1\" \
           --genotype-filter-name "isHetFilter" \
           --output {output}  
    """

################################################################################
rule filter_joint_vcf:
    input:  rules.annotate_joint_vcf.output
    output: join(KNOWN_DIR, "vcf_files",  f"{BATCH_NAME}_jointGenotypeFiltered.vcf.gz")
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
    


