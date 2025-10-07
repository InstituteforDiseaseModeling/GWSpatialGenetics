# Rules for fastqc and multiqc

# Note, these rules have not been updated to include the batch wildcard for joint batch recalibration. The assumption of downstream files for joint batch recalibration is that the bam alignment files already exist from the individual batch processing. 

################################################################################
BATCH_DIR = join(PARENT_DIR, BATCH_NAME, "out")

################################################################################
def symlink_creation(metadata, symlink_dir):
    # validate metadata column names
    # only a subset of them are necessary
    req_columns = ['ena_accession', 'sample', 'unit', 'fq1', 'fq2']
    if not(all([r in list(metadata.columns) for r in req_columns])):
        sys.exit("Metadata must contain the following columns " + str(req_columns))

    ena_ids = list(metadata['ena_accession'])
    fq1_list = list(metadata['fq1'])
    fq2_list = list(metadata['fq2'])
    # ensure no duplicates in ena_ids, fq1, fq2
    if (any(ena_ids.count(x) > 1  for x in ena_ids)):
        sys.exit("Check metadata, duplictes in ena_accession column")
    if (any(fq1_list.count(x) > 1  for x in fq1_list)):
        sys.exit("Check metadata, duplictes in fq1 column")
    if (any(fq2_list.count(x) > 1  for x in fq2_list)):
        sys.exit("Check metadata, duplictes in fq2 column")
    # ensure all files specified actually exist
    check_list = [not os.path.exists(a) for a in fq1_list+fq2_list]
    if(any(check_list)):
        missing_files = [b for a,b in zip(check_list, fq1_list + fq2_list) if a]
        print("Missing files: ")
        for a in missing_files: print(a) 
        sys.exit("Some specified fastq files do not exist")

    # Ensure all files are gzipped
    if not all([a.endswith(".gz") for a in fq1_list]) and all([a.endswith(".gz") for a in fq2_list]):
        sys.exit('All input read files must be gzipped!') 

    # get read suffix from the sample names
    # read_suffix_1 = list(set([os.path.basename(a).split(".")[0].split("_R")[-1] for a in fq1_list]))
    # read_suffix_2 = list(set([os.path.basename(a).split(".")[0].split("_R")[-1] for a in fq2_list]))
    # if len(read_suffix_1) != 1 and len(read_suffix_2) != 1:
    #     sys.exit('File read suffixes must be consistent')   

    # map from sample prefix to reads (list of two)
    read_map = {a: [b,c] for a,b,c in zip(ena_ids, fq1_list, fq2_list)}

    # create symlink dir
    if not os.path.exists(symlink_dir):
        os.makedirs(symlink_dir)
    # create symlinks for forward and reverse reads 
    for key in read_map:
        os.symlink(read_map[key][0], join(symlink_dir, "_".join([key, "R1.fastq.gz"])))
        os.symlink(read_map[key][0], join(symlink_dir, "_".join([key, "R2.fastq.gz"])))  


def get_summary_files():
    ''' Returns a list of all summary files that will be created by the single batch pipeline. '''
    QC_STEP = ['pre']
    if REMOVE_ADAPTERS is True:
        QC_STEP = QC_STEP + ['post']
    DEDUP = ['', '_dedup'] if RUN_DEDUPLICATE is True else ['']    
    
    multiqc   = expand(join(BATCH_DIR, "00_qc_reports", "{step}_multiqc","multiqc_report.html"), step=QC_STEP)
    alignment = expand(join(BATCH_DIR, "00_qc_reports", "aligned_counts{dedup}.txt"), dedup=DEDUP)
    coverage  = expand(join(BATCH_DIR, "02_align", "coverage{dedup}", "{ena_id}_pairAligned_coverage.txt"), dedup=DEDUP, ena_id = ENA_IDS)

    summary_files = multiqc + alignment + coverage

    if config['protocol'] == 'amplicon':
        summary_files = summary_files + \
            expand(join(BATCH_DIR, "02_align", "primer_counts{dedup}", "{ena_id}_primerCounts_summary.txt"), dedup=DEDUP, ena_id = ENA_IDS)

    return summary_files


################################################################################
rule read_symlinks:
    output: expand(join(BATCH_DIR, "00_read_symlinks","{ena_id}_R{read}.fastq.gz"), ena_id = ENA_IDS, read = ['1', '2'])
    run: symlink_creation(metadata, join(BATCH_DIR, "00_read_symlinks"))

################################################################################
rule pre_fastqc:
    input:
        fwd = join(BATCH_DIR, "00_read_symlinks/{ena_id}_R1.fastq.gz"),
        rev = join(BATCH_DIR, "00_read_symlinks/{ena_id}_R2.fastq.gz")
    output:
        join(BATCH_DIR,  "00_qc_reports/pre_fastqc/{ena_id}_R1_fastqc.html"),
        join(BATCH_DIR,  "00_qc_reports/pre_fastqc/{ena_id}_R2_fastqc.html")
    params:
        outdir = join(BATCH_DIR, "00_qc_reports/pre_fastqc/")
    shell: """
        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
    """

rule pre_multiqc:
    input: expand(join(BATCH_DIR, "00_qc_reports/pre_fastqc/{ena_id}_{read}_fastqc.html"), ena_id = ENA_IDS, read = ['R1', 'R2'])
    output: join(BATCH_DIR,  "00_qc_reports/pre_multiqc/multiqc_report.html")
    params:
        indir  = join(BATCH_DIR, "00_qc_reports/pre_fastqc"),
        outdir = join(BATCH_DIR, "00_qc_reports/pre_multiqc/")
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """

################################################################################
rule trim_galore:
    input:
        fwd = join(BATCH_DIR, "00_read_symlinks/{ena_id}_R1.fastq.gz"),
        rev = join(BATCH_DIR, "00_read_symlinks/{ena_id}_R2.fastq.gz")
    output:
        fwd = join(BATCH_DIR, "01_trimmed/{ena_id}_R1_val_1.fq.gz"),
        rev = join(BATCH_DIR, "01_trimmed/{ena_id}_R2_val_2.fq.gz")
    threads: 2
    params:
        q_min   = config['trim_galore']['quality'],
        min_len = config['trim_galore']['min_read_length'],
        outdir  = join(BATCH_DIR, "01_trimmed/")
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
        fwd = join(BATCH_DIR,  "00_qc_reports/post_fastqc/{ena_id}_R1_val_1_fastqc.html"),
        rev = join(BATCH_DIR,  "00_qc_reports/post_fastqc/{ena_id}_R2_val_2_fastqc.html")
    params:
        outdir = join(BATCH_DIR, "00_qc_reports/post_fastqc/")
    shell: """
        fastqc {input} --outdir {params.outdir}
    """

rule post_multiqc:
    input: 
        expand(join(BATCH_DIR,  "00_qc_reports/post_fastqc/{ena_id}_{read}_fastqc.html"), ena_id=ENA_IDS, read=['R1_val_1', 'R2_val_2'])
    output: join(BATCH_DIR, "00_qc_reports/post_multiqc/multiqc_report.html")
    params:
        indir  = join(BATCH_DIR,  "00_qc_reports/post_fastqc"),
        outdir = join(BATCH_DIR,  "00_qc_reports/post_multiqc/")
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """

