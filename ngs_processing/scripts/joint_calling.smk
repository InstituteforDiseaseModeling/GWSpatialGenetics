'''
Snakemake pipeline to combine batches of sequencing and do joint variant calling
Combined to call single batch processing (gw_processing.snakefile)

Authors - Jessica Ribado and Ben Siranosian
'''

################################################################################
def remove_cumulative_duplicates(gvcf_list_directory, batch_gvcf_list, \
    cumulative_gvcf_path, extension="gvcfs.txt"):
    """
    1) Finds all files under gvcf_list_directory ending with the given extension.
    2) Reads batch_gvcf_list and collects its first-column ('sample') values as a set.
    3) For each existing gvcf list file prior to the current batch:
         - Drops any row where the first-column value ('sample') is in the batch set
         - If a duplicate sample is found in the current batch, removes the first instance of that sample from the previous batch list since the latest batch will include a merged version of the sample across batches. (See merge_batch_replicates.smk for more details)
    """

    # read only the first column of batch_file
    batch_df = pd.read_csv(batch_gvcf_list, sep = "\t", header=None, dtype=str)
    batch_samples = set(batch_df.iloc[:,0])

    pattern = glob.glob(join(gvcf_list_directory, f"*{extension}"))
    files = [f for f in pattern if f not in batch_gvcf_list]
    for fn in files:
        df = pd.read_csv(fn, dtype=str, sep="\t", header=None)
        
        # drop any rows in the old file where it matches the current batch 
        mask = df.iloc[:, 0].isin(batch_samples)
        if mask.any():
            df = df.loc[~mask]
            df.to_csv(join(gvcf_list_directory, basename(fn)), index=False, header=False, sep="\t")

        # save output to a cumulative gvcf file 
        df.to_csv(cumulative_gvcf_path[0], mode="a", index=False, header=False, sep="\t")

    # add the current batch gvcf list to the cumulative gvcf file 
    batch_df.to_csv(cumulative_gvcf_path[0], mode="a", index=False, header=False, sep="\t")

    return 


################################################################################
rule cumulative_gvcf_list:
    input: 
        gvcf_file = rules.gvcf_path_file.output,
    output: temp(join(KNOWN_DIR, "gvcf_lists", f"{BATCH_NAME}_cumulative.txt"))
    params:
        gvcf_list_dir = join(KNOWN_DIR, "gvcf_lists") 
    run:
        # check if this is single batch run, where across batch duplicates were identified
        batch_merge = glob.glob(join(PARENT_DIR, BATCH_NAME, "out", "03_variant_calls", "multibatch_duplicates", f"*.g.vcf.gz"))
        print(f"Batch merge file: {batch_merge}")
        if len(batch_merge) > 0:
            print("Batch merge file exists, removing cumulative duplicates")
            remove_cumulative_duplicates(params.gvcf_list_dir, input[0], output)
        else:
            print("No batch duplicates found, creating cumulative gVCF list")
            shell("cd {params.gvcf_list_dir}; cat *_gvcfs.txt >> {output}")

################################################################################
rule GenomicsDBImport_all:
    input: rules.cumulative_gvcf_list.output
    output: join(KNOWN_DIR, "genomicsDB", "callset.json")
    threads: 10    
    params:
        db_dir = join(KNOWN_DIR, "genomicsDB")
    shell: """
        rm -rf {params.db_dir}
        # set intervals based on the headers of the chromosome in the fasta file
        intervals=$(grep ">" {REF_FILE} | cut -f 1 -d " " | tr -d ">" | tr "\n" "," | sed "s/,$//g")
        gatk GenomicsDBImport \
            --reference {REF_FILE} \
            --genomicsdb-workspace-path {params.db_dir} \
            --batch-size 100 \
            -L $intervals \
            --sample-name-map {input} \
            --reader-threads {threads} 
    """

################################################################################
rule joint_genotyping_all:
    input: rules.GenomicsDBImport_all.output
    output: join(KNOWN_DIR, "vcf_files", f"{BATCH_NAME}_jointGenotype.vcf.gz")
    params:
        db_dir = join(KNOWN_DIR, "genomicsDB")
    shell: """
        gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant gendb://{params.db_dir} \
            --output {output} \
    """

