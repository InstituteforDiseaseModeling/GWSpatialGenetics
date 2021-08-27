# Rules for fastqc and multiqc

################################################################################
rule read_symlinks:
    input: 
    output: expand(join(PROJECT_DIR, "00_read_symlinks/{ena_id}_R{read}.fastq.gz", ena_id = ena_ids, read = ['1', '2'])
    shell: """
        ln -s {input.fwd} {output.fwd}
        ln -s {input.rev} {output.rev}
    """

rule pre_fastqc:
    input:
        fwd = join(PROJECT_DIR, "00_read_symlinks/{ena_id}_R1.fastq.gz"),
        rev = join(PROJECT_DIR, "00_read_symlinks/{ena_id}_R2.fastq.gz")
    output:
        join(PROJECT_DIR,  "00_qc_reports/pre_fastqc/{ena_id}_R1_fastqc.html"),
        join(PROJECT_DIR,  "00_qc_reports/pre_fastqc/{ena_id}_R2_fastqc.html")
    params:
        outdir = join(PROJECT_DIR, "00_qc_reports/pre_fastqc/")
    shell: """
        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
    """

rule pre_multiqc:
    input: expand(join(PROJECT_DIR, "00_qc_reports/pre_fastqc/{ena_id}_{read}_fastqc.html"), ena_id=ena_ids, read=['R1', 'R2'])
    output: join(PROJECT_DIR,  "00_qc_reports/pre_multiqc/multiqc_report.html")
    params:
        indir  = join(PROJECT_DIR, "00_qc_reports/pre_fastqc"),
        outdir = join(PROJECT_DIR, "00_qc_reports/pre_multiqc/")
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """

################################################################################
rule trim_galore:
    input:
        fwd = join(PROJECT_DIR, "00_read_symlinks/{ena_id}_R1.fastq.gz")
        rev = join(PROJECT_DIR, "00_read_symlinks/{ena_id}_R2.fastq.gz")
    output:
        fwd = join(PROJECT_DIR, "01_trimmed/{ena_id}_R1_val_1.fq" + gz_ext),
        rev = join(PROJECT_DIR, "01_trimmed/{ena_id}_R2_val_2.fq" + gz_ext)
    threads: 2
    params:
        q_min   = config['trim_galore']['quality'],
        min_len = config['trim_galore']['min_read_length'],
        outdir  = join(PROJECT_DIR, "01_trimmed/")
    log:
        join(PROJECT_DIR, "logs/{ena_id}_trim.log")
    shell: """
        mkdir -p {params.outdir}
        trim_galore --quality {params.q_min} \
            --length {params.min_len} \
            --cores {threads} \
            --output_dir {params.outdir} \
            --paired {input.fwd} {input.rev}
    """

################################################################################
rule post_fastqc:
    input:
        fwd = rules.trim_galore.output.fwd,
        rev = rules.trim_galore.output.rev
    output: 
        fwd = join(PROJECT_DIR,  "00_qc_reports/post_fastqc/{ena_id}_R1_val_1_fastqc.html"),
        rev = join(PROJECT_DIR,  "00_qc_reports/post_fastqc/{ena_id}_R2_val_2_fastqc.html")
    params:
        outdir = join(PROJECT_DIR, "00_qc_reports/post_fastqc/")
    shell: """
        fastqc {input} --outdir {params.outdir}
    """

rule post_multiqc:
    input: 
        expand(join(PROJECT_DIR,  "00_qc_reports/post_fastqc/{ena_id}_{read}_fastqc.html"), ena_id=ena_ids, read=['R1_val_1', 'R2_val_2'])
    output: join(PROJECT_DIR, "00_qc_reports/post_multiqc/multiqc_report.html")
    params:
        indir  = join(PROJECT_DIR,  "00_qc_reports/post_fastqc"),
        outdir = join(PROJECT_DIR,  "00_qc_reports/post_multiqc/")
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """


################################################################################
# final quality check after all other files have been generated
# rule post_qc_report:
#     input: expand(join(PROJECT_DIR, "03_variant_calls/jointGenotype_{iteration}Iter_filtered.{extension}"), iteration=list(range(1,config['max_iter']+1)), extension=['vcf.gz', 'txt']),
#     output: 
#         pdf = join(PROJECT_DIR, "mtDNA_quality_report.pdf"),
#     params: 
#         project_dir = PROJECT_DIR,
#         primer_fasta_f = config['primer_file'],
#         output_low_depth = join(PROJECT_DIR, "mtDNA_poor_breadth_samples.txt")
#     conda: "../ngs_qc_env.yaml"
#     script: 
#         "mtDNA_qualityCheck_auto.Rmd"
