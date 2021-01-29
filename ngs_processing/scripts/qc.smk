# Rules for fastqc and multiqc

################################################################################
rule pre_fastqc:
    input:
        fwd = lambda wildcards: read_map[wildcards.ena_id][0],
        rev = lambda wildcards: read_map[wildcards.ena_id][1],
    output:
        join(PROJECT_DIR,  "00_qc_reports/pre_fastqc/{ena_id}_" + READ_SUFFIX[0] + "_fastqc.html"),
        join(PROJECT_DIR,  "00_qc_reports/pre_fastqc/{ena_id}_" + READ_SUFFIX[1] + "_fastqc.html")
    params:
        outdir = join(PROJECT_DIR, "00_qc_reports/pre_fastqc/")
    shell: """
        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
    """

rule pre_multiqc:
    input: expand(join(PROJECT_DIR, "00_qc_reports/pre_fastqc/{ena_id}_{read}_fastqc.html"), ena_id=ena_ids, read=READ_SUFFIX)
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
        fwd = lambda wildcards: read_map[wildcards.ena_id][0],
        rev = lambda wildcards: read_map[wildcards.ena_id][1],
    output:
        fwd = join(PROJECT_DIR, "01_trimmed/{ena_id}_" + READ_SUFFIX[0] + "_val_1.fq" + gz_ext),
        rev = join(PROJECT_DIR, "01_trimmed/{ena_id}_" + READ_SUFFIX[1] + "_val_2.fq" + gz_ext)
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
        fwd = join(PROJECT_DIR, "01_trimmed/{ena_id}_" + READ_SUFFIX[0] + "_val_1.fq" + gz_ext),
        rev = join(PROJECT_DIR, "01_trimmed/{ena_id}_" + READ_SUFFIX[1] + "_val_2.fq" + gz_ext)
    output: 
        fwd = join(PROJECT_DIR,  "00_qc_reports/post_fastqc/{ena_id}_" + READ_SUFFIX[0] + "_val_1_fastqc.html"),
        rev = join(PROJECT_DIR,  "00_qc_reports/post_fastqc/{ena_id}_" + READ_SUFFIX[1] + "_val_2_fastqc.html")
    params:
        outdir = join(PROJECT_DIR, "00_qc_reports/post_fastqc/")
    shell: """
        fastqc {input} --outdir {params.outdir}
    """

rule post_multiqc:
    input: 
        fwd = expand(join(PROJECT_DIR,  "00_qc_reports/post_fastqc/{ena_id}_" + READ_SUFFIX[0] + "_val_1_fastqc.html"), ena_id=ena_ids),
        rev = expand(join(PROJECT_DIR,  "00_qc_reports/post_fastqc/{ena_id}_" + READ_SUFFIX[1] + "_val_2_fastqc.html"), ena_id=ena_ids)
    output: join(PROJECT_DIR, "00_qc_reports/post_multiqc/multiqc_report.html")
    params:
        indir  = join(PROJECT_DIR,  "00_qc_reports/post_fastqc"),
        outdir = join(PROJECT_DIR,  "00_qc_reports/post_multiqc/")
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """
