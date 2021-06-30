# Rules for fastqc and multiqc

################################################################################
rule read_symlinks:
    input: 
        fwd = lambda wildcards: read_map[wildcards.ena_id][0],
        rev = lambda wildcards: read_map[wildcards.ena_id][1],
    output: 
        fwd = join(PROJECT_DIR, "00_read_symlinks/{ena_id}_" + READ_SUFFIX[0] + ".fastq" + gz_ext),
        rev = join(PROJECT_DIR, "00_read_symlinks/{ena_id}_" + READ_SUFFIX[1] + ".fastq" + gz_ext)
    shell: """
        ln -s {input.fwd} {output.fwd}
        ln -s {input.rev} {output.rev}
    """
rule pre_fastqc:
    input:
        fwd = rules.read_symlinks.output.fwd,
        rev = rules.read_symlinks.output.rev
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
        fwd = rules.read_symlinks.output.fwd,
        rev = rules.read_symlinks.output.rev,
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


################################################################################
# final quality check after all other files have been generated
rule post_qc_report:
    input: expand(join(PROJECT_DIR, "03_variant_calls/jointGenotype_{iteration}Iter_filtered.{extension}"), iteration=list(range(1,config['max_iter']+1)), extension=['vcf.gz', 'txt']),
    output: 
        pdf = join(PROJECT_DIR, "mtDNA_quality_report.pdf"),
    params: 
        project_dir = PROJECT_DIR,
        primer_fasta_f = config['primer_file'],
        output_low_depth = join(PROJECT_DIR, "mtDNA_poor_breadth_samples.txt")
    conda: "../ngs_qc_env.yaml"
    script: 
        "mtDNA_qualityCheck_auto.Rmd"
