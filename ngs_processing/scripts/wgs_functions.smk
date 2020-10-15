###############################################################################
# load python modules
import fnmatch as fn
import pandas as pd
from os.path import join, splitext
from collections import defaultdict
from snakemake.io import Wildcards, expand
import snakemake.rules

################################################################################
# specify project directories and files
DATA_DIR    = config["raw_reads_directory"]
PROJECT_DIR = config["output_directory"]
REF_DIR     = config['reference_alignment']['genome_path']
REF_NAME    = config['reference_alignment']['genome_name']
REF_FILE    = join(REF_DIR, REF_NAME)
META_FILE   = pd.read_csv(config["metadata_file"], delimiter = '\t')

################################################################################
# read in file names and paths from the config files used in multiple rules
SAMPLE_PREFIX = list(set(META_FILE['sample']))
ENA_ID = list(set(META_FILE['ena_accession']))

################################################################################
def assign_sample_replicates(metafile):
    ''' Create a dictionary where the sample prefix is the key and corresponding sample replicates are the samples. '''
    samp_dict = defaultdict(list)
    #seq_pairs = metafile['sample']
    seq_pairs = zip(metafile['sample'], metafile['ena_accession'])
    for sample, ena in seq_pairs:
        samp_dict[sample].append(ena)
    return(samp_dict)
    #return(seq_pairs)


################################################################################
def recursive_bam(iteration):
    n = int(iteration.iteration)
    if n == 0:
        return join(PROJECT_DIR, "01_processing/02_align/recalibrate/%s_%dIter.bam")  % (iteration.sample, n)
    if n > 0:
        return join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_%d/%s_%dIter.bam") % (n-1, iteration.sample, n-1)
    else:
        raise ValueError("Iteration steps must be an integer: received %s" % iteration)

# ################################################################################
# def recursive_vcf(iteration):
#     n = int(iteration.iteration)
#     if n == 0:
#         return rules.call_haplotypes.output.vcf
#     if n > 0:
#         return join(PROJECT_DIR, "01_processing/02_align/recalibrate/haplocall_%d/%s_%dIter_filtQD.vcf") % (n-1, iteration.sample, n-1)
#     else:
#         raise ValueError("Iteration steps must be an integer: received %s" % iteration)
