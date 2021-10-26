# Rules for fastqc and multiqc

################################################################################
rule read_symlinks:
    output: expand(join(PROJECT_DIR, "00_read_symlinks/{ena_id}_R{read}.fastq.gz"), ena_id = ena_ids, read = ['1', '2'])
    run: symlink_creation(metadata, join(PROJECT_DIR, "00_read_symlinks"))

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
        fwd = join(PROJECT_DIR, "00_read_symlinks/{ena_id}_R1.fastq.gz"),
        rev = join(PROJECT_DIR, "00_read_symlinks/{ena_id}_R2.fastq.gz")
    output:
        fwd = join(PROJECT_DIR, "01_trimmed/{ena_id}_R1_val_1.fq.gz"),
        rev = join(PROJECT_DIR, "01_trimmed/{ena_id}_R2_val_2.fq.gz")
    threads: 2
    params:
        q_min   = config['trim_galore']['quality'],
        min_len = config['trim_galore']['min_read_length'],
        outdir  = join(PROJECT_DIR, "01_trimmed/")
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

