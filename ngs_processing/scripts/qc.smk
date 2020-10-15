################################################################################
# extract naming convention from file names
files = fn.filter(os.listdir(DATA_DIR), "*.fastq*")
#READ_SUFFIX = list(set([x.split(".")[0] for x in files]))
READ_SUFFIX = ["R1", "R2"]
EXTENSION = ".fastq.gz"
gz_ext = '.gz' if EXTENSION.endswith('.gz') else ''


################################################################################
rule pre_fastqc:
	input:
		join(DATA_DIR, "{ena_id}_" + READ_SUFFIX[0] + EXTENSION),
		join(DATA_DIR, "{ena_id}_" + READ_SUFFIX[1] + EXTENSION)
	output:
		join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{ena_id}_" + READ_SUFFIX[0] + "_fastqc.html"),
		join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_fastqc/{ena_id}_" + READ_SUFFIX[1] + "_fastqc.html")
	params:
		outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc/")
	shell: """
		mkdir -p {params.outdir}
		fastqc {input} --outdir {params.outdir}
	"""

rule pre_multiqc:
	input: expand(join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc/{ena_id}_{read}_fastqc.html"), ena_id=ENA_ID, read=READ_SUFFIX)
	output: join(PROJECT_DIR,  "01_processing/00_qc_reports/pre_multiqc/multiqc_report.html")
	params:
		indir  = join(PROJECT_DIR, "01_processing/00_qc_reports/pre_fastqc"),
		outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/pre_multiqc/")
	shell: """
		multiqc --force {params.indir} -o {params.outdir}
	"""

################################################################################
rule trim_galore:
	input:
		fwd = join(DATA_DIR, "{ena_id}_" + READ_SUFFIX[0] + EXTENSION),
		rev = join(DATA_DIR, "{ena_id}_" + READ_SUFFIX[1] + EXTENSION)
	output:
		fwd = join(PROJECT_DIR, "01_processing/01_trimmed/{ena_id}_" + READ_SUFFIX[0] + "_val_1.fq" + gz_ext),
		rev = join(PROJECT_DIR, "01_processing/01_trimmed/{ena_id}_" + READ_SUFFIX[1] + "_val_2.fq" + gz_ext)
	threads: 4
	params:
		q_min   = config['trim_galore']['quality'],
		# left    = config['trim_galore']['start_trim'],
        # right   = config['trim_galore']['end_trim'],
		min_len = config['trim_galore']['min_read_length'],
		outdir  = join(PROJECT_DIR, "01_processing/01_trimmed/")
	shell: """
		mkdir -p {params.outdir}
		trim_galore --quality {params.q_min} \
			--length {params.min_len} \
			--output_dir {params.outdir} \
			--paired {input.fwd} {input.rev}
	"""
#--clip_R1 {params.left} \
#--clip_R2 {params.left} \
#--three_prime_clip_R1 {params.right} \
#--three_prime_clip_R2 {params.right} \

################################################################################
rule post_fastqc:
	input:  rules.trim_galore.output
	output: join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc/{ena_id}_val_{read}_fastqc.html")
	params:
		outdir = join(PROJECT_DIR, "01_processing/00_qc_reports/post_fastqc/")
	shell: """
		fastqc {input} --outdir {params.outdir}
	"""

rule post_multiqc:
	input: expand(join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc/{ena_id}_{read}_val_{read}_fastqc.html"), ena_id=ENA_ID, read=READ_SUFFIX),
	output: join(PROJECT_DIR, "01_processing/00_qc_reports/post_multiqc/multiqc_report.html")
	params:
		indir  = join(PROJECT_DIR,  "01_processing/00_qc_reports/post_fastqc"),
		outdir = join(PROJECT_DIR,  "01_processing/00_qc_reports/post_multiqc/")
	shell: """
		multiqc --force {params.indir} -o {params.outdir}
	"""
