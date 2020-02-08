sample_dict=assign_sample_replicates(META_FILE)

################################################################################
rule align_to_ref:
    ''' Aligns reads to specified reference genome and. Currently only retains and sorts reads where both read mates align to the reference. '''
    input: bwa_build = rules.build_ref_index.output, \
        bwa_index = REF_FILE, \
        fwd = rules.trim_galore.output.fwd, \
		rev = rules.trim_galore.output.rev
    output:
        join(PROJECT_DIR, "01_processing/02_align/{ena_id}_pairAligned.bam")
    params:
        rg=r"@RG\tID:{ena_id}\tSM:{ena_id}"
    threads: 12
    shell: """
        bwa mem -R '{params.rg}' -t {threads} {input.bwa_index} {input.fwd} {input.rev} | samtools view -b -F12 - | samtools sort -o {output} -
    """

################################################################################
rule mark_duplicates:
    input: rules.align_to_ref.output
    output:
        dupl_bam = join(PROJECT_DIR, "01_processing/02_align/derep/{ena_id}_pairAligned_duplMarked.bam"),
        dupl_log = join(PROJECT_DIR, "01_processing/02_align/derep/{ena_id}_pairAligned_duplMarked.log")
    threads: 8
    shell: """
        picard MarkDuplicates \
            INPUT={input} \
            OUTPUT={output.dupl_bam} \
            REMOVE_DUPLICATES=true \
            METRICS_FILE={output.dupl_log}
    """

################################################################################
rule merge_sample_bams:
    ''' Merges aligned reads for each replicate. '''
    input:
        lambda wildcards: expand(join(PROJECT_DIR, "01_processing/02_align/derep/{ena_id}_pairAligned_duplMarked.bam"), ena_id=sample_dict[wildcards.sample])
    output:
        join(PROJECT_DIR, "01_processing/02_align/merged/{sample}_merged.bam")
    run:
        if len(input) > 1:
            shell("samtools merge {output} {input}")
        else:
            shell("cp {input} {output} && touch -h {output}")


################################################################################
rule merge_add_groups:
    ''' Adds a new variable to specify combined file name to bam file. Neccessary for calling variants with GATK. '''
    input:  rules.merge_sample_bams.output
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/{sample}_0Iter.bam")
    shell: """
        picard AddOrReplaceReadGroups \
            I={input} \
            O={output} \
            RGLB=library1 \
            RGPL=ILLUMINA \
            RGSM={wildcards.sample} \
            RGPU={wildcards.sample}
    """

################################################################################
rule merge_bam_index:
    input: rules.merge_add_groups.output
    output: join(PROJECT_DIR, "01_processing/02_align/recalibrate/{sample}_0Iter.bam.bai")
    shell: """
        samtools index {input} > {output}
    """
