'''
Snakemake pipeline to discover variants in non-model organisms using GATK Best Practices.

Assumes some sequencing batches have already been aligned and deduplicated (optional) and will be finding alignment files across batches. 

Author - Jessica Ribado 
'''

################################################################################
def map_samples_to_batches_df(df, unique=True):
    batch_rx = re.compile(r'/batch_[^/]+')
    mapping = defaultdict(list)

    for _, row in df.iterrows():
        sample = row['sample']
        fq1    = row.get('fq1', '') or ''
        m = batch_rx.search(fq1)
        if not m:
            continue
        batch = m.group(0).lstrip('/')

        if unique:
            if batch not in mapping[sample]:
                mapping[sample].append(batch)
        else:
            mapping[sample].append(batch)

    return dict(mapping)


def get_iterative_bam(wc):
    base      = wc.sample
    batch     = wc.batch
    i         = int(wc.iteration)

    if i == 1:
        return join(PARENT_DIR, batch, "out", "02_align", "recalibrate", f"{base}_0Iter.bam")
    else:
        return join(PARENT_DIR, batch, "out", "02_align", "recalibrate", f"bsqr_{i-1}", f"{base}.bam")


################################################################################
checkpoint check_multisamples:
    output: 
        merged_dir = directory(join(PARENT_DIR, RECALIBRATION_BATCH, "out", "02_align", "recalibrate", "multibatch_duplicates"))
    run:
        samples_by_batch = map_samples_to_batches_df(metadata)
        samples_with_replicates = {
            sample: batches
            for sample, batches in samples_by_batch.items()
            if len(batches) > 1
        }

        os.makedirs(output.merged_dir, exist_ok=True)

        # 3) if no replicates, bail out
        if not samples_with_replicates:
            print("No samples with replicates found. Skipping merge step.")
            return

        # 4) otherwise, write one TSV per sample
        for sample, batches in samples_with_replicates.items():
            bam_paths = [
                join(
                    PARENT_DIR, batch,
                    "out", "02_align", "recalibrate",
                    f"{sample}_0Iter.bam"
                )
                for batch in batches
            ]
            out_file = join(output.merged_dir, f"{sample}_bamLocations.tsv")
            with open(out_file, "w") as fh:
                for p in bam_paths:
                    fh.write(f"{sample}\t{p}\n")



################################################################################
rule merge_multibatch_bams:
    input: 
        merged_tsv = lambda wc: os.path.join(
            checkpoints.check_multisamples.get().output.merged_dir,
            f"{wc.sample}_bamLocations.tsv"
        )
    output: 
        bam = join(PARENT_DIR, RECALIBRATION_BATCH, "out", "02_align", "recalibrate", "{sample}_0Iter.bam"),
        bai = join(PARENT_DIR, RECALIBRATION_BATCH, "out", "02_align", "recalibrate", "{sample}_0Iter.bam.bai"),
    shell: """
        mapfile -t bam_array < <(awk -F'\t' '{{print $2}}' {input.merged_tsv})
        samtools merge {output.bam} $(echo "${{bam_array[*]}}")
        samtools index {output.bam} > {output.bai}
    """ 


################################################################################
rule call_haplotypes:
    input:
        bam = get_iterative_bam
    output: join(PARENT_DIR, "{batch}", "out", "03_variant_calls", "discovery", "bsqr_{iteration}", "{sample}.g.vcf.gz")
    params:
        ploidy = 2
    threads: 2
    shell:  """
        mkdir -p $(dirname {output})
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
        gvcfs = lambda wc: [
            join(PARENT_DIR, bd, "out", "03_variant_calls", "discovery", f"bsqr_{wc.iteration}", f"{s}.g.vcf.gz")
            for bd in list(SAMPLES_PER_BATCH.keys())
            for s  in SAMPLES_PER_BATCH[bd]
        ]
    output: join(DISCOVERY_DIR, BATCH_NAME, "out", "bsqr_{iteration}", "genomicsDB", "callset.json")
    params:
        db_dir   = lambda wc: join(DISCOVERY_DIR, BATCH_NAME, "out", f"bsqr_{wc.iteration}", "genomicsDB"),
        vcf_args = lambda wc, input: ' '.join(f"-V {g}" for g in input.gvcfs)
    shell: """
        rm -rf {params.db_dir}

        # list all contigs from the FASTA header
        intervals=$(grep '>' {REF_FILE} | cut -f1 -d' ' | tr -d '>' | paste -sd, -)

        gatk GenomicsDBImport \
          --reference {REF_FILE} \
          {params.vcf_args} \
          --genomicsdb-workspace-path {params.db_dir} \
          -L "$intervals"
    """

################################################################################
rule joint_genotyping:
    input: rules.GenomicsDBImport.output
    output: join(DISCOVERY_DIR, BATCH_NAME, "out", "bsqr_{iteration}", "jointGenotype.vcf.gz")
    params: 
        db_dir = join(DISCOVERY_DIR, BATCH_NAME, "out", "bsqr_{iteration}", "genomicsDB"),
    shell: """
        gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant gendb://{params.db_dir} \
            --output {output} \
    """

################################################################################
rule filter_discovery_vcf:
    input:  rules.joint_genotyping.output
    output: join(DISCOVERY_DIR, BATCH_NAME, "out", "bsqr_{iteration}", "jointGenotypeAnnotated.vcf.gz")    
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
    output: join(DISCOVERY_DIR, BATCH_NAME, "out", "bsqr_{iteration}", "jointGenotypeFiltered.vcf.gz")
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
        bam = get_iterative_bam,
        discovery_vcf = lambda wc: join(DISCOVERY_DIR, BATCH_NAME, "out", f"bsqr_{wc.iteration}", "jointGenotypeFiltered.vcf.gz")
    output:
        table = join(PARENT_DIR, "{batch}", "out", "02_align", "recalibrate", "bsqr_{iteration}", "{sample}.table"),
        bam   = join(PARENT_DIR, "{batch}", "out", "02_align", "recalibrate", "bsqr_{iteration}", "{sample}.bam")
        #table = f"{PARENT_DIR}/{{batch}}/out/02_align/recalibrate/bsqr_{{iteration}}/{{sample}}.table",
        #bam   = f"{PARENT_DIR}/{{batch}}/out/02_align/recalibrate/bsqr_{{iteration}}/{{sample}}.bam"
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
        initial_recal = join(PARENT_DIR, "{batch}", "out", "02_align", "recalibrate", "bsqr_1", "{sample}.table"),
        final_recal   = join(PARENT_DIR, "{batch}", "out", "02_align", "recalibrate", "bsqr_{iteration}", "{sample}.table")
    output: join(PARENT_DIR, "{batch}", "out", "02_align", "recalibrate", "summary", "{sample}_{iteration}_bqsrCovariates.pdf")
    params: 
        output_dir = join(PARENT_DIR, "{batch}", "out", "02_align", "recalibrate", "summary")
    shell: """
        mkdir -p {params.output_dir}
        gatk AnalyzeCovariates \
          -before {input.initial_recal} \
          -after  {input.final_recal} \
          -plots  {output}
    """

################################################################################
rule bsqr_path_file:
    input:
        gvcfs = [
            join(PARENT_DIR, f"{batch}", "out", "03_variant_calls", "discovery", "bsqr_{iteration}", f"{s}.g.vcf.gz")
        for batch in SAMPLES_PER_BATCH.keys()
        for s  in SAMPLES_PER_BATCH[batch]
    ]
    output: 
        sample_list = join(DISCOVERY_DIR, BATCH_NAME, "out", "bsqr_{iteration}", "genomicsDBimport_gvcf_list.txt")
    run:
        with open(output.sample_list, 'w') as f:
            for gvcf in input.gvcfs:
                sample = os.path.basename(gvcf).split(".")[0]
                f.write(f"{sample}\t{gvcf}\n")

     
################################################################################
rule archive_preBSQR_files:
    input: join(DISCOVERY_DIR, RECALIBRATION_BATCH, "out", f"bsqr_{MAX_ITER}", "genomicsDBimport_gvcf_list.txt")
    output:
        archive_marker = join(KNOWN_DIR, "gvcf_lists", f"preBSQR_{RECALIBRATION_BATCH}_archive", ".archived")
    params:
        source = join(KNOWN_DIR, "gvcf_lists"),
    run:
        destination = os.path.dirname(output.archive_marker)
        os.makedirs(destination, exist_ok=True)
        if not exists(input[0]):
            print(f"Reference file {input[0]} not found – skipping archive.")
        else:
            shell(
                "find {params.source} -type f ! -newer {input[0]} "
                "-exec mv {{}} {params.source}/preBSQR_{BATCH_NAME}_archive \\;"
            )
        
        shell("touch {output.archive_marker}")


################################################################################
rule gvcf_path_file:
    """
    Copy the gvcf list file to a new location for joint calling. Will make all these files accessible for variant calling with future batches. 
    """
    input: 
        archive_stamp = rules.archive_preBSQR_files.output.archive_marker,
        batch_list = join(DISCOVERY_DIR, RECALIBRATION_BATCH, "out", f"bsqr_{MAX_ITER}", "genomicsDBimport_gvcf_list.txt")
    output: join(KNOWN_DIR, "gvcf_lists", f"{RECALIBRATION_BATCH}_genomicDB_gvcfs.txt")
    run:
        shell("cp {input.batch_list} {output}")    