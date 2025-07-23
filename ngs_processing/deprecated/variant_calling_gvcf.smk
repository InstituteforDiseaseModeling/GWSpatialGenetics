################################################################################
rule call_haplotypes:
    input:
        bam = get_bam,
        bam_idx = rules.merge_add_groups.output.bai, # not directly used in rule
        refdict = rules.create_ref_dict.output, # not directly used in rule
        faidx = rules.build_ref_faidx.output    # not directly used in rule
    output: join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz")
    threads: 2
    params:
        #ploidy = config['variant_calling']['ploidy']
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
rule GenomicsDBImport:
    input: 
        lambda wildcards: expand(join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz"), sample=unique_samples, iteration=wildcards.iteration)
    output: join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/genomicsDB/callset.json")
    params:
        db_dir = join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/genomicsDB"),
        gvcf_string = lambda wildcards: expand("-V " + join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz"), sample=unique_samples, iteration=wildcards.iteration)
    shell: """
        rm -rf {params.db_dir}
        # set intervals based on the headers of the chromosome in the fasta file
        intervals=$(grep ">" {REF_FILE} | cut -f 1 -d " " | tr -d ">" | tr "\n" "," | sed "s/,$//g")
        gatk GenomicsDBImport \
            --reference {REF_FILE} \
            {params.gvcf_string} \
            --genomicsdb-workspace-path {params.db_dir} \
            -L $intervals
    """

################################################################################
rule joint_genotyping:
    input: rules.GenomicsDBImport.output
    output: join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/joint_genotype.vcf.gz")
    params: 
        db_dir = join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/genomicsDB"),
    shell: """
        gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant gendb://{params.db_dir} \
            --output {output} \
    """

################################################################################
rule filter_discovery_vcf:
    input:  rules.joint_genotyping.output
    output: join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/joint_genotypeAnnotated.vcf.gz")    
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
rule snp_discovery_vcf:
    input:  rules.filter_discovery_vcf.output
    output: join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/joint_genotypeFiltered.vcf.gz")
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

###############################################################################
rule calculate_apply_discovery_base_recalibration:
    input:
        bam = get_bam,
        discovery_vcf = join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/joint_genotypeFiltered.vcf.gz")
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
rule check_base_discovery_recalibration:
    input:
        initial_recal = join(PROJECT_DIR, "02_align/recalibrate/haplocall_0/{sample}_0Iter.table"),
        final_recal  = join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.table")
    output: join(PROJECT_DIR, "02_align/recalibrate/{sample}_{iteration}Iter_bqsrCovariates.pdf")
    shell: """
        gatk AnalyzeCovariates \
            -before {input.initial_recal} \
            -after {input.final_recal} \
            -plots {output}
    """

################################################################################
rule gvcf_recal_file:
    input: expand(join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz"), iteration = MAX_ITER, sample = list(set(metadata['sample'])))
    output: join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}", "genomicsDBimport_gvcf_list.txt")
    params:     
        gvcf_dir = expand(join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}"), iteration = MAX_ITER)[0]
    run:
        gvcf_list_df(str(params.gvcf_dir), str(output))  