'''
Snakemake pipeline rules to identify technical replicates across batches, merge, call variants, and output the gvcf list file with the correct paths.

Authors - Jessica Ribado 
'''

###############################################################################
def check_db_sample_names(batch_gvcf_list, db_dir, output_dir):
    # Determine current batch sample from input files
    os.makedirs(output_dir, exist_ok=True)
    
    batch_name = re.search(r"batch_[^/]+", batch_gvcf_list[0]).group(0)
    batch_samples = [basename(x).split('.', 1)[0] for x in batch_gvcf_list]
    
    # read & concat all mapping files
    frames = []
    files = glob.glob(join(db_dir, '*_gvcfs.txt'))
    # confirm that current sequencing batch is excluded, if re-running
    files = [fn for fn in files if batch_name not in fn]  
    for fpath in files:
        df = pd.read_csv(
            fpath,
            sep=r"\s+",
            header=None,
            names=['sample', 'gvcf_path'],
            engine='python'
        )
        frames.append(df)
    if not frames:
        print(f"No GenomicDB mapping files in {db_dir}")  
        return  
    else:
        mapping_df = pd.concat(frames, ignore_index=True)
    
    # check if any sample in the batch is already in the database
    overlaps = mapping_df[mapping_df['sample'].isin(batch_samples)] 
    if overlaps.empty:
        print("No overlaps found with existing samples in the database.")
        return
    else:
        print("Samples in the sequencing batch already exist in the database:")
        print(overlaps)

        # Group by sample and write only those with matches
        for sample, group in overlaps.groupby('sample'):
            batch_gvcf = dirname(batch_gvcf_list[0])
            batch_rows = pd.DataFrame({
                'sample': [sample],
                'gvcf_path': [f"{batch_gvcf}/{sample}.g.vcf.gz"] 
            })
            sample_df = pd.concat([group, batch_rows], ignore_index=True)
            
            regex_map = {
                r"03_variant_calls/known_variants": "02_align/recalibrate", 
                r"\.g\.vcf\.gz$": "_0Iter.bam"}

            sample_df['gvcf_path'] = sample_df['gvcf_path'].replace(regex_map, regex=True)
            out_file = join(output_dir, f"{sample}_merged.txt")
            sample_df.to_csv(out_file, sep='\t', header=False, index=False)
        return


def gvcf_list(wildcards):
    # 1) wait for the checkpoint to finish and fetch the directory of your TXT manifests
    merged_txt_dir = checkpoints.check_db_samples.get(**wildcards).output[0]

    # 2) collect the already-existing GVCFs
    known_dir = os.path.join(
        PARENT_DIR, BATCH_NAME,
        "out", "03_variant_calls", "known_variants"
    )
    originals = glob.glob(join(known_dir, "*.g.vcf.gz"))
    files = originals

    # 3) look for the merged‐sample TXT files
    txt_paths = glob.glob(join(merged_txt_dir, "*_merged.txt"))
    merged_gvcfs = []
    if txt_paths:
        # strip “.txt” to get your sample names
        merged_samples = [basename(p).rsplit("_merged.txt", 1)[0]
            for p in txt_paths]
        # build the paths where those g.vcf.gz’s will live
        dup_dir = os.path.join(
            PARENT_DIR, BATCH_NAME,
            "out", "03_variant_calls", "multibatch_duplicates"
        )
        merged_gvcfs = [
            os.path.join(dup_dir, f"{s}.g.vcf.gz")
            for s in merged_samples
        ]
        files = originals + merged_gvcfs

    # 4) return one flat list of everything
    return files


################################################################################
checkpoint check_db_samples:
    input: 
        gvcfs = expand(join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "known_variants", "{sample}.g.vcf.gz"), sample=UNIQUE_SAMPLES)
    output: 
        merged_dir = directory(join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "multibatch_duplicates"))
    params:
        gvcf_dir = join(KNOWN_DIR, "gvcf_lists")   
    run: 
        os.makedirs(output.merged_dir, exist_ok=True)
        if os.path.isdir(params.gvcf_dir): 
            if os.listdir(params.gvcf_dir):
                check_db_sample_names(input.gvcfs, params.gvcf_dir, output.merged_dir)
 

################################################################################
rule merge_multibatch_bams:
    ''' Merges aligned deduplicated reads for each replicate. '''
    input:  
        merged_txt = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "multibatch_duplicates", "{sample}_merged.txt") 
    output: 
        bam = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "multibatch_duplicates", "{sample}_merged.bam"),
        bai = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "multibatch_duplicates", "{sample}_merged.bam.bai"),
    shell: """
        mapfile -t bam_array < <(awk -F'\t' '{{print $2}}' {input.merged_txt})
        samtools merge {output.bam} $(echo "${{bam_array[*]}}")
        samtools index {output.bam} > {output.bai}
    """  

################################################################################
rule merge_multibatch_calls:
    ''' Merges aligned deduplicated reads for each replicate. '''
    input:  
        bam = join(PARENT_DIR, BATCH_NAME, "out", "02_align", "recalibrate", "multibatch_duplicates", "{sample}_merged.bam")
    output: join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "multibatch_duplicates", "{sample}.g.vcf.gz")
    threads: 2
    params:
        # ploidy = config['variant_calling']['ploidy']
        ploidy = 2
    shell: """
        gatk HaplotypeCaller \
            --reference {REF_FILE} \
            --input {input.bam} \
            --output {output} \
            -ploidy {params.ploidy} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads}
    """

################################################################################
rule gvcf_path_file:
    input: 
        batch_file = rules.gvcf_batch_file.output.list_file,
        sample_vcfs = gvcf_list
    output: join(KNOWN_DIR, "gvcf_lists", f"{BATCH_NAME}_genomicDB_gvcfs.txt")
    run:
        # Find all merged gVCFs
        merged_dir = os.path.join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "multibatch_duplicates")
        merged_files = glob.glob(os.path.join(merged_dir, "*.g.vcf.gz"))

        if merged_files:
            print(f"[INFO] Found merged duplicates: {merged_files}")

            batch_gvcf_df = pd.read_csv(input[0], sep="\t", header=None, names=["sample", "gvcf"])

            # Get just the merged sample names
            merged_sample_names = [os.path.basename(f).split(".")[0] for f in merged_files]

            # Remove rows from batch_gvcf_df where the gvcf path ends in a duplicate sample
            batch_no_duplicates = batch_gvcf_df[
                ~batch_gvcf_df["sample"].isin(merged_sample_names)
            ]

            # Build new DataFrame for merged samples
            merged_df = pd.DataFrame({
                "sample": merged_sample_names,
                "gvcf": merged_files
            })

            # Combine and write
            updated = pd.concat([batch_no_duplicates, merged_df])
            updated.to_csv(output[0], sep="\t", index=False, header=False)
        else:
            print(f"[INFO] No merged duplicates found. Copying original list.")
            shell("cp {input[0]} {output[0]}")
