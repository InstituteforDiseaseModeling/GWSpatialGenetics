###############################################################################
# load python modules
import fnmatch as fn
import pandas as pd
from os.path import join, splitext
from collections import defaultdict
from snakemake.io import Wildcards, expand
import snakemake.rules
import sys

################################################################################
# specify project directories and files
PROJECT_DIR = config["output_directory"]
REF_FILE = config['reference_file']
metadata = pd.read_csv(config["metadata_file"], delimiter = '\t')
# validate metadata column names
# if (list(metadata.columns) != ['ena_accession', 'sample', 'unit', 'fq1', 'fq2', 'host', 'country', 'variant.data']):
#    sys.exit("Metadata column names must be exactly 'ena_accession', 'sample', 'unit', 'fq1', 'fq2', 'host', 'country', 'variant.data'")

################################################################################
# read in file names and paths from the config files used in multiple rules
ena_ids = list(metadata['ena_accession'])
sample_list = list(metadata['sample'])
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

# map from sample prefix to reads (list of two)
read_map = {a: [b,c] for a,b,c in zip(ena_ids, fq1_list, fq2_list)}

# check all are .gz or not
if all([a.endswith(".gz") for a in fq1_list]) and all([a.endswith(".gz") for a in fq2_list]):
    gz_ext = '.gz'
elif all([not a.endswith(".gz") for a in fq1_list]) and all([not a.endswith(".gz") for a in fq2_list]): 
    gz_ext = ''
else: 
    sys.exit('All filenames must be gzipped, or not gzipped. No mixing.')

# get read suffix from the sample names
read_suffix_1 = list(set([os.path.basename(a).split(".")[0].split("_")[-1] for a in fq1_list]))
read_suffix_2 = list(set([os.path.basename(a).split(".")[0].split("_")[-1] for a in fq2_list]))
if len(read_suffix_1) == 1 and len(read_suffix_2) == 1:
    READ_SUFFIX = [read_suffix_1[0], read_suffix_2[0]]
else: 
    sys.exit('File read suffixes must be consistent')

################################################################################
def assign_sample_replicates(metafile):
    ''' Create a dictionary where the sample prefix is the key and corresponding sample replicates are the samples. '''
    samp_dict = defaultdict(list)
    seq_pairs = zip(metafile['sample'], metafile['ena_accession'])
    for sample, ena in seq_pairs:
        samp_dict[sample].append(ena)
    return(samp_dict)

################################################################################
def recursive_bam(iteration):
    n = int(iteration.iteration)
    if n == 0:
        return join(PROJECT_DIR, "02_align/recalibrate/%s_%dIter.bam")  % (iteration.sample, n)
    if n > 0:
        return join(PROJECT_DIR, "02_align/recalibrate/haplocall_%d/%s_%dIter.bam") % (n-1, iteration.sample, n-1)
    else:
        raise ValueError("Iteration steps must be an integer: received %s" % iteration)

