################################################################################
def singlebatch_replicate_dicts(df):
    ''' Create a dictionary where the sample prefix is the key and corresponding sample duplicates are the samples. '''
    samp_dict = defaultdict(list)
    seq_pairs = zip(df['sample'], df['ena_accession'])
    for sample, ena in seq_pairs:
        samp_dict[sample].append(ena)
    
    return samp_dict


def get_fastq():
    if REMOVE_ADAPTERS is True:
        return rules.trim_galore.output
    else:
        return expand(join(PARENT_DIR, BATCH_NAME, "out", "00_read_symlinks/{{ena_id}}_R{read}.fastq.gz"), read=['1', '2'])


def get_aligned_bams():
    if RUN_DEDUPLICATE is True:
        return lambda wildcards: expand(join(PARENT_DIR, BATCH_NAME, "out", "02_align", "dedup", "{ena_id}_pairAligned_duplMarked.bam"), ena_id = SAMPLE_DICT[wildcards.sample])
    else:
        return lambda wildcards: expand(join(PARENT_DIR, BATCH_NAME, "out", "02_align", "align", "{ena_id}_pairAligned.bam"), ena_id = SAMPLE_DICT[wildcards.sample]) 


################################################################################
SAMPLE_DICT = singlebatch_replicate_dicts(metadata)

################################################################################
rule align_to_ref:
    ''' Aligns reads to specified reference genome and. Currently only retains and sorts reads where both read mates align to the reference. '''
    input:
        bwa_build = rules.build_ref_index.output, 
        bwa_index = REF_FILE, 
        fastq = get_fastq()
    output:
        join(PARENT_DIR, BATCH_NAME, "out", "02_align", "align", "{ena_id}_pairAligned.bam")
    params:
        rg=r"@RG\tID:{ena_id}\tSM:{ena_id}"
    threads: 12
    shell: """
        bwa mem -R '{params.rg}' -t {threads} {input.bwa_index} {input.fastq} | samtools view -b -F12 - | samtools sort -o {output} -
    """

################################################################################
rule count_aligned_reads:
    input:
        expand(join(PARENT_DIR, BATCH_NAME, "out", "02_align", "align", "{ena_id}_pairAligned.bam"), ena_id=ENA_IDS)
    output:
        join(PARENT_DIR, BATCH_NAME, "out", '00_qc_reports', 'aligned_counts.txt')
    shell: """
        for i in {input}; do
            samtools index "$i"
            bn=$(basename "$i" | sed 's/_pairAligned.bam//g')
            echo "$bn\t$(samtools view -c -F 260 $i)" >> {output}
        done
    """

################################################################################
rule mark_duplicates:
    input: rules.align_to_ref.output
    output:
        dupl_bam = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "derep", "{ena_id}_pairAligned_duplMarked.bam"),
        dupl_log = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "derep", "{ena_id}_pairAligned_duplMarked.log")
    threads: 4
    shell: """
        picard MarkDuplicates \
            INPUT={input} \
            OUTPUT={output.dupl_bam} \
            REMOVE_DUPLICATES=true \
            METRICS_FILE={output.dupl_log}
    """

################################################################################
# this now counts after the MarkDuplicates stage
rule count_aligned_reads_dedup:
    input:
        expand(join(PARENT_DIR, BATCH_NAME, "out", "02_align", "derep", "{ena_id}_pairAligned_duplMarked.bam"), ena_id=ENA_IDS)
    output:
        join(PARENT_DIR, BATCH_NAME, "out", "00_qc_reports", 'aligned_counts_dedup.txt')
    shell: """
        for i in {input}; do
            samtools index "$i"
            bn=$(basename "$i" | sed 's/_pairAligned_duplMarked.bam//g')
            echo "$bn\t$(samtools view -c -F 260 $i)" >> {output}
        done
    """

################################################################################
rule merge_sample_bams:
    ''' Merges aligned deduplicated reads for each replicate. '''
    input: get_aligned_bams()
    output:
        join(PARENT_DIR, BATCH_NAME, "out", "02_align", "merged", "{sample}_merged.bam")
    run:
        if len(input) > 1:
            shell("samtools merge {output} {input}")
        else:
            shell("echo {input}; cp {input} {output} && touch -h {output}")

################################################################################
rule merge_add_groups:
    ''' Adds a new variable to specify combined file name to bam file. Neccessary for calling variants with GATK. '''
    input: rules.merge_sample_bams.output
    output: 
        bam = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "{sample}_0Iter.bam"),
        bai = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "{sample}_0Iter.bam.bai")
    shell: """
        picard AddOrReplaceReadGroups \
            I={input} \
            O={output.bam} \
            RGLB=library1 \
            RGPL=ILLUMINA \
            RGSM={wildcards.sample} \
            RGPU={wildcards.sample}
        samtools index {output.bam} > {output.bai}    
    """
