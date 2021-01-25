################################################################################
rule genome_coverage:
    input:  rules.align_to_ref.output
    output: join(PROJECT_DIR, "02_align/coverage/{ena_id}_pairAligned_coverage.txt"),
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
#     output: join(PROJECT_DIR, "02_align/primer_counts/{ena_id}_primerCounts_mismatch{mismatch}.txt")
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
    output: temp(join(PROJECT_DIR, "02_align/{ena_id}_pairAligned_filt.fq"))
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
    output: join(PROJECT_DIR, "02_align/primer_counts/{ena_id}_primerCounts.txt")
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
    output: join(PROJECT_DIR, "02_align/primer_counts/{ena_id}_primerCounts_summary.txt")
    threads: 1
    shell: """
        seqkit fx2tab {input.primers} | cut -f 2 |  while read pattern <&0; do echo -e $pattern"\t"$(grep -c $pattern {input.counts}); done > {output}
    """
