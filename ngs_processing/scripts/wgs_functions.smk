###############################################################################
# load python modules
import fnmatch as fn
import pandas as pd
import numpy as np
import re
import csv
#from pandas.io.json import json_normalize
from os.path import join, splitext, exists, basename, dirname, getsize, isfile, relpath
from collections import defaultdict
from snakemake.io import Wildcards, expand
import snakemake.rules
import sys
import json
import glob

################################################################################
# Specify project directories and files
################################################################################
PARENT_DIR = config["output_directory"].split("/batch", 1)[0]
BATCH_NAME = re.search(r'(batch[^/]+)', config["output_directory"]).group(1)

REF_FILE = config['reference_file']
VARIANT_FILE = config['variant_calling']['known_variants']

# read in metadata file that contains sample information
metadata = pd.read_csv(config['metadata_file'], delimiter = '\t')#.head(2)
UNIQUE_SAMPLES = list(metadata['sample'].unique())
ENA_IDS = list(metadata['ena_accession'])
    
# Specific to primers - update for NGS but not worth 
if config['protocol'] == 'amplicon':
    REMOVE_ADAPTERS = True
    RUN_DEDUPLICATE = False

wildcard_constraints:
    batch = r"[^/]+",
    sample = "[^/]+",
    merged_sample = "[^/]+",
    iteration = "[0-9]+",
    n = "[0-9]+",

################################################################################
# Specify final files based on user defined inputs
# Updates 2025-06:
# Paths for discovery variant calling output names to harmonize with known variant calling output names
################################################################################
def map_samples_to_batches_df(df, unique=True):
    """
    Same mapping, but takes a DataFrame instead of a file path.
    """
    batch_rx = re.compile(r'/batch_[^/]+')
    mapping = defaultdict(list)

    for _, row in df.iterrows():
        sample = row['sample']
        fq1    = row.get('fq1', '') or ''
        m = batch_rx.search(fq1)
        if not m:
            continue
        batch = m.group(0).lstrip('/')

        if unique:
            if batch not in mapping[sample]:
                mapping[sample].append(batch)
        else:
            mapping[sample].append(batch)

    return dict(mapping)


def invert_with_special_key(samples_per_batch, special_key):
    """
    Given sample→[batches], return a dict mapping:
      - each batch that has exactly one sample → [that sample]
      - special_key → [all samples that had >1 batch]
    """
    inverted = {}
    multi = []
    for sample, batches in samples_per_batch.items():
        if len(batches) > 1:
            multi.append(sample)
        else:
            # exactly one batch → map batch→[sample]
            batch = batches[0]
            inverted.setdefault(batch, []).append(sample)
    # finally, if we found any multi-batch samples, group them under special_key
    if multi:
        inverted[special_key] = multi
    return inverted    


###########################################################################
KNOWN_DIR = join(PARENT_DIR, "joint_calling_known")
if isfile(VARIANT_FILE):
    MAX_ITER = 0    
else:
    DISCOVERY_DIR = join(PARENT_DIR, "joint_calling_discovery")
    RECALIBRATION_BATCH = BATCH_NAME
    MAX_ITER = config['max_iter']
    ITERATIONS = list(range(1, MAX_ITER+1)) 
    SAMPLES_BY_BATCH = map_samples_to_batches_df(metadata)
    SAMPLES_PER_BATCH = invert_with_special_key(SAMPLES_BY_BATCH, BATCH_NAME)

    print(SAMPLES_PER_BATCH)
    

###########################################################################
def user_outputs():
    ''' Returns a list of all files that will be created by the pipeline. '''
    if isfile(VARIANT_FILE):
        summary_files = [
            join(PARENT_DIR, BATCH_NAME, "out", "00_qc_reports", "mtDNA_quality_report.pdf")] + \
            [join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "summary", f"{s}_known_bqsrCovariates.pdf") for s in UNIQUE_SAMPLES] 
        sample_vcf_files = [join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "known_variants", f"{s}.g.vcf.gz") for s in UNIQUE_SAMPLES] 
        batch_vcf_lists = [join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "known_variants", "genomicsDBimport_gvcf_list.txt")]
    else:
        summary_files = [
            join(PARENT_DIR, batch, "out", "02_align", "recalibrate", "summary", f"{s}_{it}_bqsrCovariates.pdf") 
            for batch in SAMPLES_PER_BATCH.keys()
            for it in ITERATIONS if it > 1
            for s  in SAMPLES_PER_BATCH[batch]
        ]    

        sample_vcf_files = [
            join(PARENT_DIR, batch, "out", "03_variant_calls", "discovery", f"bsqr_{it}", f"{s}.g.vcf.gz")
            for batch in SAMPLES_PER_BATCH.keys()
            for it in ITERATIONS
            for s  in SAMPLES_PER_BATCH[batch]
        ]

        batch_vcf_lists = [
            join(DISCOVERY_DIR, BATCH_NAME, "out", f"bsqr_{it}", f"{ext}")
            for it in ITERATIONS
            for ext in ["genomicsDBimport_gvcf_list.txt", "jointGenotypeFiltered.vcf.gz"]
        ]

    joint_vcf = [join(KNOWN_DIR, "vcf_files",  f"{BATCH_NAME}_jointGenotypeFiltered.vcf.gz")]


    all_files = (
        summary_files + 
        sample_vcf_files + 
        batch_vcf_lists + 
        joint_vcf       
    )
        

    return all_files





