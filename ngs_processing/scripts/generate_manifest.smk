import os, re, sys
import numpy as np
import pandas as pd
from collections import Counter

mapping_file = config['sample_key'] if os.path.isfile(config['sample_key']) else None

#######################################################################
def manifest_generate(fastq_path, output, key_file):
    r1_reads = [f for f in os.listdir(fastq_path) if re.search(r'^(?!Undetermined)[\w\.-]*R1', f)]
    print("First few files in fastq directory:", r1_reads[0:7])
    reads_df = pd.DataFrame({
        'fq1':[os.path.join(fastq_path, y) for y in r1_reads], 
        'fq2':[os.path.join(fastq_path, re.sub("R1", "R2", y)) for y in r1_reads]
    })

    if key_file is not None:
        key_file = pd.read_csv(key_file, sep=",")
        print("Header of key file:")
        print(key_file.head())
        if 'sample' and 'sample_number' not in key_file.columns.tolist():
            sys.exit("Sample key file does not contain sample_number (Qiagen match) and sample (rename) column(s). Update sample key file and try again.")
        
        reads_df['sample_number'] = [int(x.split("-")[1].lstrip('0')) for x in r1_reads]
        reads_df = key_file.merge(reads_df, on='sample_number')
    else:    
        print("Sample key file not provided. Proceeding.")
        # Parsing for Cornell - assumes the metadata and worm number are in the 5th and 6th position
        sample = ["_".join(y.split("_")[4:6]) for y in r1_reads]
        # check that sample name contains year; if not, rename samples by name order for matching in in analysis step
        print('Checking putative sample names follow standard naming convention...')
        sample_ranges = list(map(str, list(range(2010, 2025)) + ["MISC"]))
        sample_matches = [s for s in sample if any(xs in s for xs in sample_ranges)]
        if len(sample) > len(sample_matches):
            print("Sample names do not follow standard naming convention. Updating sample names to numeric order.")
            sample = ["s"+str(x) for x in list(range(1,len(sample)+1))]
        reads_df['sample'] = sample

        
    duplicates = [item for item, count in Counter(reads_df['sample']).items() if count > 1] 
    reads_df['unit'] = reads_df.groupby(['sample']).cumcount() + 1
    reads_df['ena_accession'] = reads_df["sample"] + "-" + reads_df["unit"].astype(str)

    print("Last few lines of the manifest file. Outputting to", output)
    print(reads_df.tail())
    reads_df.to_csv(output, sep="\t", index=False)

#######################################################################
rule all:
    input:
        os.path.join(config['output_file'])

#######################################################################
rule create_manifest:
    input:
        fastq_dir = config['fastq_dir']
    output: os.path.join(config['output_file'])
    run:
            manifest_generate(str(input.fastq_dir), str(output), mapping_file) 
           
