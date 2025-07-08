
# Rules that that mirror iterative base recalibration calling with with known variants.

###############################################################################
rule calculate_apply_known_base_recalibration:
    input:
        bam = rules.merge_add_groups.output.bam,
        known_vcf = config['variant_calling']['known_variants']
    output: 
        table = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "known_variants", "{sample}.table"),
        bam = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "known_variants", "{sample}.bam")
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
    output: join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "summary", "{sample}_known_bqsrCovariates.pdf")
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
    output: join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "known_variants", "{sample}.g.vcf.gz")
    threads: 2
    params:
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
rule gvcf_batch_file:
    input:
        gvcfs = expand(join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "known_variants", "{sample}.g.vcf.gz"), sample=UNIQUE_SAMPLES)
    output:
        list_file = join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "known_variants", "genomicsDBimport_gvcf_list.txt")
    run:
        gvcf_files = list(map(str, input.gvcfs))  
        gvcf_prefixes = [
            re.split(r'_[0-9]+Iter|\.g\.vcf(?:\.gz)?$', os.path.basename(f))[0]
            for f in gvcf_files
        ]
        gvcf_df = pd.DataFrame(np.c_[gvcf_prefixes, gvcf_files])
        gvcf_df.to_csv(output[0], index=False, header=False, sep="\t") 
         


################################################################################
# additional rules for optional single batch joint genotype calling
################################################################################
rule batch_joint_genotyping:
    input: rules.gvcf_batch_file.output
    output: join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "known_variants", "jointGenotype.vcf.gz")
    params: 
        db_dir = join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "known_variants", "genomicsDB"),
    shell: """
        gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant gendb://{params.db_dir} \
            --output {output} \
    """

