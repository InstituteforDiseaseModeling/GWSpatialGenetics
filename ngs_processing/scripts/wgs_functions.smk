###############################################################################
# load python modules
import fnmatch as fn
import pandas as pd
import re
from pandas.io.json import json_normalize
from os.path import join, splitext, exists, basename, dirname, getsize
from collections import defaultdict
from snakemake.io import Wildcards, expand
import snakemake.rules
import sys
import json
import glob

################################################################################
# Specify project directories and files
################################################################################
PROJECT_DIR = config["output_directory"]
REF_FILE = config['reference_file']
VARIANT_FILE = config['variant_calling']['known_variants']
FILE_STRUCTURE = PROJECT_DIR.rsplit('/', 1)
PARENT_DIR = FILE_STRUCTURE[0]
BATCH_NAME = FILE_STRUCTURE[1]
MAX_ITER = 0 if os.path.isfile(VARIANT_FILE) else config['max_iter']

# Specific to primers - update for NGS but not worth 
if config['protocol'] == 'amplicon':
    REMOVE_ADAPTERS = True
    RUN_DEDUPLICATE = False

metadata = pd.read_csv(config['metadata_file'], delimiter = '\t')#.head(2)
ena_ids = list(metadata['ena_accession'])
unique_samples = list(metadata['sample'].unique())

wildcard_constraints:
    sample="[^/]+",
    # iteration="[1-9][0-9]*",
    iteration="[0-9]+",
    n="[0-9]+",

################################################################################
# Specify final files based on user defined inputs
# Updates 2021-09 include: 
# - adapter trimming (amplicons)
# - specify inclusion of deplicates (amplicons)
# - skip iterative recalibration with a set of known variants
################################################################################
def user_outputs():

    final_files = [join(JOINT_DIR, "vcf_files",  f"{BATCH_NAME}_jointGenotypeFiltered.vcf.gz")] + \
        [join(PROJECT_DIR, "00_qc_reports", "mtDNA_quality_report.pdf")]

    if os.path.isfile(VARIANT_FILE):
        variant_files = expand(join(PROJECT_DIR, "03_variant_calls/known_variants/{sample}.g.vcf.gz"), sample = list(set(metadata['sample']))) + \
            expand(join(PROJECT_DIR, "02_align/recalibrate/{sample}_known_bqsrCovariates.pdf"), sample = list(set(metadata['sample']))) + \
            [join(PROJECT_DIR, "03_variant_calls", "known_variants", "genomicsDBimport_gvcf_list.txt")] \
            # + [join(PROJECT_DIR, "03_variant_calls", "known_variants", "joint_genotypeFiltered.vcf.gz")]
    else:
        variant_files = expand(join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/{sample}_{iteration}Iter.bam"), \
                iteration = list(range(0, MAX_ITER+1)), sample = list(set(metadata['sample']))) + \
            expand(join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}/{sample}_{iteration}Iter.g.vcf.gz"),  \
                iteration = list(range(0, MAX_ITER+1)), sample = list(set(metadata['sample']))) + \
            expand(join(PROJECT_DIR, "02_align/recalibrate/{sample}_{iteration}Iter_bqsrCovariates.pdf"), \
                iteration = list(range(0, MAX_ITER+1)), sample = list(set(metadata['sample']))) + \
            expand(join(PROJECT_DIR, "03_variant_calls/haplocall_{iteration}", "genomicsDBimport_gvcf_list.txt"), iteration = np.arange(0,MAX_ITER))

    all_files = final_files + variant_files
    #all_files = [join(PROJECT_DIR, "00_qc_reports", "mtDNA_quality_report.pdf")]

    return all_files


################################################################################
# Functions called throughout single-batch pipeline
################################################################################
def assign_sample_replicates(metafile):
    ''' Create a dictionary where the sample prefix is the key and corresponding sample duplicates are the samples. '''
    samp_dict = defaultdict(list)
    seq_pairs = zip(metafile['sample'], metafile['ena_accession'])
    for sample, ena in seq_pairs:
        samp_dict[sample].append(ena)
    return(samp_dict)

################################################################################
def get_fastq():
    if REMOVE_ADAPTERS is True:
        return rules.trim_galore.output
    else:
        return expand(join(PROJECT_DIR, "00_read_symlinks/{{ena_id}}_R{read}.fastq.gz"), read=['1', '2'])

################################################################################
def get_aligned_bams():
    if RUN_DEDUPLICATE is True:
        return lambda wildcards: expand(join(PROJECT_DIR, "02_align/dedup/{ena_id}_pairAligned_duplMarked.bam"), ena_id=sample_dict[wildcards.sample])
    else:
        return lambda wildcards: expand(join(PROJECT_DIR, "02_align/align/{ena_id}_pairAligned.bam"), ena_id=sample_dict[wildcards.sample]) 

################################################################################    
def get_bam(wildcards):
    n = int(wildcards.iteration)
    # print(wildcards)
    if n == 0:
        return join(PROJECT_DIR, "02_align/recalibrate/%s_%dIter.bam")  % (wildcards.sample, n)
        # sys.exit('shouldnt be called with zero!')
    if n > 0:
        return join(PROJECT_DIR, "02_align/recalibrate/haplocall_%d/%s_%dIter.bam") % (n-1, wildcards.sample, n-1)
    else:
        raise ValueError("Iteration steps must be an integer: received %s" % iteration)        

################################################################################
def get_discovery_vcf():
    if os.path.isfile(VARIANT_FILE):
        return VARIANT_FILE
    else:
        return join(PROJECT_DIR, "02_align/recalibrate/haplocall_{iteration}/joint_genotypeFiltered.vcf.gz")    

################################################################################
def get_summary_files():
    QC_STEP = ['pre']
    if REMOVE_ADAPTERS is True:
        QC_STEP = QC_STEP + ['post']
    DEDUP = ['', '_dedup'] if RUN_DEDUPLICATE is True else ['']    
    
    multiqc   = expand(join(PROJECT_DIR, "00_qc_reports/{step}_multiqc/multiqc_report.html"), step=QC_STEP)
    alignment = expand(join(PROJECT_DIR, "00_qc_reports/aligned_counts{dedup}.txt"), dedup=DEDUP)
    coverage  = expand(join(PROJECT_DIR, "02_align/coverage{dedup}/{ena_id}_pairAligned_coverage.txt"), dedup=DEDUP, ena_id = ena_ids)

    summary_files = multiqc + alignment + coverage

    if config['protocol'] == 'amplicon':
        summary_files = summary_files + \
            expand(join(PROJECT_DIR, "02_align/primer_counts{dedup}/{ena_id}_primerCounts_summary.txt"), dedup=DEDUP, ena_id = ena_ids)

    return summary_files

       
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

################################################################################
# Functions called throughout joint-batch pipeline
################################################################################
def gvcf_list_df(gvcf_dir, output_df) :
    gvcf_files = glob.glob(join(gvcf_dir, '*.vcf.gz'))
    gvcf_prefix = [re.split(r'(?<=\w)_[0-9]Iter|.g.vcf', basename(f))[0] for f in gvcf_files]
    gvcf_df = pd.DataFrame(np.c_[gvcf_prefix, gvcf_files])
    gvcf_df.to_csv(output_df, index=False, header=False, sep="\t") 

################################################################################
def get_gvcf_list():
    if os.path.isfile(VARIANT_FILE):
        return(join(PROJECT_DIR, "03_variant_calls", "known_variants", "genomicsDBimport_gvcf_list.txt"))
    else:    
        return(join(PROJECT_DIR, "03_variant_calls", f"haplocall_{MAX_ITER}", "genomicsDBimport_gvcf_list.txt"))

################################################################################
def get_genomicDB_names(name_json):
    with open(name_json, 'r') as j:
        contents = json.loads(j.read())
        df = pd.json_normalize(contents, record_path = ['callsets'])
    return df['sample_name'].tolist()    

################################################################################
def check_sample_names(gvcf_list, db_dir):
    print("Checking names in gvcf list", gvcf_list)
    gvcf_df = pd.read_csv(gvcf_list, sep="\t", names=['sample', 'gvcf_path'])
    print("First few lines of batch gvcf file:")
    print(gvcf_df.head())
    
    gvcf_files = glob.glob(join(db_dir, "gvcf_lists", "*_gvcfs.txt"))
    
    # check the same batch name hasn't been run before - avoids duplicates if rerun without deleting GVCF sample 
    gvcf_rm_batch = [x for x in gvcf_files if BATCH_NAME not in x]
    if len(gvcf_rm_batch) < len(gvcf_files):
        print("Current batch has been previously run.")
        print("Ignoring samples with the same batch name", BATCH_NAME, "for checking duplicates.")

    if len(gvcf_files) > 0:
        joint_tmp = pd.DataFrame()
        for file in gvcf_rm_batch:
            frame = pd.read_csv(file, delimiter="\t", names=['sample', 'gvcf_path'])
            joint_tmp = joint_tmp.append(frame)

        # check for duplicate files in current database
        dup_samples = [x for x in gvcf_df['sample'].tolist() if x in joint_tmp['sample'].tolist()]
        if len(dup_samples) > 0:
            print(len(dup_samples), "duplicate sample(s) detected in newest batch. Identifying batches for duplicate sample files.")
            dup_filt  = gvcf_df[gvcf_df['sample'].isin(dup_samples)]
            joint_dup = joint_tmp[joint_tmp['sample'].isin(dup_samples)]
            dup_df = pd.concat([dup_filt, joint_dup])    

            # get the location of the bam files
            dup_df = dup_df.replace("03_variant_calls", "02_align/recalibrate", regex=True)
            dup_df = dup_df.replace(".g.vcf.gz", ".bam", regex=True)

            # create directory if it does not exist
            if not exists(join(db_dir, BATCH_NAME)):
                os.makedirs(join(db_dir, BATCH_NAME))   
            # create mapping files for each duplicated sample
            for dup_samp in dup_df['sample'].unique():
                dup_tmp = dup_df[dup_df['sample'] == dup_samp].drop_duplicates()
                dup_tmp['gvcf_path'].to_csv(join(db_dir, BATCH_NAME, f"{dup_samp}.txt"), index=False, header=False, sep="\t")

            gvcf_df = gvcf_df[~gvcf_df['sample'].isin(dup_samples)]
    else:
        print("GenomicsDB currently empty. Seeding with", gvcf_list, "sample file.")   
    
    gvcf_df.to_csv(join(db_dir, BATCH_NAME, "gvcf_tmp.tsv"), index=False, header=False, sep="\t") 

################################################################################
def gvcf_list(wildcards):
    checkpoint_output = checkpoints.check_samples.get(**wildcards).output[0]
    print(checkpoint_output)
    putative_files = os.listdir(dirname(checkpoint_output))
    print(putative_files)
    if len(putative_files) > 1:
        merged = expand(join(JOINT_DIR, BATCH_NAME, "{merged_sample}.g.vcf.gz"),
                        merged_sample=glob_wildcards(join(dirname(checkpoint_output), "{merged_sample}.txt")).merged_sample)
        return [checkpoint_output] + merged
    else: 
        return checkpoint_output

################################################################################
def update_gvcf_list(files, output):
    # read in files that did not have a duplicate
    files = files.split(" ")
    if getsize(files[0]) != 0:
        tmp_files = pd.read_csv(files[0], header=None, sep="\t")
    else:    
        tmp_files = pd.DataFrame()    
    
    # update file paths for new .g.vcf files for merged samples
    merged_files  = files[1:]
    merged_prefix = [re.split(r'(?<=\w)_[0-9]Iter|.g.vcf', basename(f))[0] for f in merged_files]
    merged_prefix = [f"{x}-{BATCH_NAME}" for x in merged_prefix]
    merged_df = pd.DataFrame(np.column_stack((merged_prefix, merged_files)))    
    # concatenate and save
    gvcf_final = pd.concat([tmp_files, merged_df])
    gvcf_final.to_csv(output, index=False, header=False, sep="\t")
    



