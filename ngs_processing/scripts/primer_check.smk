################################################################################
rule genome_coverage:
    input:  rules.align_to_ref.output
    output: join(PARENT_DIR, BATCH_NAME, "out", "02_align", "coverage", "{ena_id}_pairAligned_coverage.txt"),
    threads: 1
    shell: """
        samtools depth -aa -Q 20 {input} > {output}
    """
################################################################################
rule genome_coverage_dedup:
    input:  rules.mark_duplicates.output.dupl_bam
    output: join(PARENT_DIR, BATCH_NAME, "out", "02_align", "coverage_dedup", "{ena_id}_pairAligned_coverage.txt"),
    threads: 1
    shell: """
        samtools depth -aa -Q 20 {input} > {output}
    """

# trash?
# ################################################################################
# rule primer_counts:
#     input:
#         primers = config['primer_check']['primer_file'],
#         align = rules.align_to_ref.output
#     output: join(PARENT_DIR, BATCH_NAME, "out", "02_align", "primer_counts", "{ena_id}_primerCounts_mismatch{mismatch}.txt")
#     params:
#         map_qual = config['trim_galore']['quality'],
#         mismatch  = config['primer_check']['mismatch']
#     threads: 4
#     shell: """
#          mapfile -t primers < {input.primers}
#          rm -f {output}
#          for i in ${{primers[@]}}; do
#             count=($(samtools view -F 260 -q {params.map_qual} {input.align} | tre-agrep -{params.mismatch} $i | cut -f4 | paste -sd',')); echo "$i    $count" >> {output};
#         done;
#     """

################################################################################
rule filter_fastq:
    input: rules.align_to_ref.output
    output: temp(join(PARENT_DIR, BATCH_NAME, "out", "02_align", "{ena_id}_pairAligned_filt.fq"))
    params:
        map_qual = config['trim_galore']['quality']
    threads: 4
    shell: """
         samtools view -F 260 -q {params.map_qual} -S -b {input} | samtools bam2fq - > {output}
    """

# should be included
# this is the combined version that goes into the Rmkdwn script
# don't need positions where they align?
################################################################################
rule primer_counts:
    input:
        primers = config['primer_file'],
        align = rules.filter_fastq.output
    output: join(PARENT_DIR, BATCH_NAME, "out", "02_align", "primer_counts", "{ena_id}_primerCounts.txt")
    params:
        map_qual = config['trim_galore']['quality']
    threads: 2
    shell: """
        seqkit locate --degenerate --pattern-file {input.primers} {input.align} > {output}
    """

################################################################################
rule primer_summary:
    input:
        primers = config['primer_file'],
        counts = rules.primer_counts.output
    output: join(PARENT_DIR, BATCH_NAME, "out", "02_align", "primer_counts", "{ena_id}_primerCounts_summary.txt")
    threads: 1
    shell: """
        seqkit fx2tab {input.primers} | cut -f 2 |  while read pattern <&0; do echo -e $pattern"\t"$(grep -c $pattern {input.counts}); done > {output}
    """

################################################################################
rule filter_fastq_dedup:
    input: rules.mark_duplicates.output.dupl_bam
    output: temp(join(PARENT_DIR, BATCH_NAME, "out", "02_align", "{ena_id}_pairAligned_dedup_filt.fq"))
    params:
        map_qual = config['trim_galore']['quality']
    threads: 4
    shell: """
         samtools view -F 260 -q {params.map_qual} -S -b {input} | samtools bam2fq - > {output}
    """

################################################################################
rule primer_counts_dedup:
    input:
        primers = config['primer_file'],
        align = rules.filter_fastq_dedup.output
    output: join(PARENT_DIR, BATCH_NAME, "out", "02_align", "primer_counts_dedup", "{ena_id}_primerCounts.txt")
    params:
        map_qual = config['trim_galore']['quality']
    threads: 2
    shell: """
        seqkit locate --degenerate --pattern-file {input.primers} {input.align} > {output}
    """

################################################################################
rule primer_summary_dedup:
    input:
        primers = config['primer_file'],
        counts = rules.primer_counts_dedup.output
    output: join(PARENT_DIR, BATCH_NAME, "out", "02_align" ,"primer_counts_dedup", "{ena_id}_primerCounts_summary.txt")
    threads: 1
    shell: """
        seqkit fx2tab {input.primers} | cut -f 2 |  while read pattern <&0; do echo -e $pattern"\t"$(grep -c $pattern {input.counts}); done > {output}
    """

################################################################################
# final quality check after all other files have been generated
rule mt_qc_report:
    input: get_summary_files()
    output: 
        pdf = join(PARENT_DIR, BATCH_NAME, "out", "00_qc_reports", "mtDNA_quality_report.pdf"),
    params: 
        project_dir = join(PARENT_DIR, BATCH_NAME, "out"),
        primer_fasta_f = config['primer_file'],
        output_low_depth = join(PARENT_DIR, BATCH_NAME, "out", "00_qc_reports", "mtDNA_poor_breadth_samples.txt")
    # conda: join(workflow.basedir, "ngs_qc_env.yaml")
    script: 
        "mtDNA_qualityCheck_auto.Rmd"