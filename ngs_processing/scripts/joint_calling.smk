'''
Snakemake pipeline to combine batches of sequencing and do joint variant calling
Combined to call single batch processing (gw_processing.snakefile)

Authors - Jessica Ribado and Ben Siranosian
'''

################################################################################
import sys
from os.path import exists, join

JOINT_DIR =  join(PARENT_DIR, "joint_calling_known") if os.path.isfile(VARIANT_FILE) else join(PARENT_DIR, "joint_calling_recal")

################################################################################
checkpoint check_samples:
    input: get_gvcf_list()
    output: 
        sample_list = join(JOINT_DIR, BATCH_NAME, "gvcf_tmp.tsv")
    run: 
        check_sample_names(str(input), JOINT_DIR)

################################################################################
rule merge_batch_bams:
    ''' Merges aligned deduplicated reads for each replicate. '''
    input:  join(JOINT_DIR, BATCH_NAME, "{merged_sample}.txt") 
    output: 
        bam = join(JOINT_DIR, BATCH_NAME, "{merged_sample}.bam"),
        bai = join(JOINT_DIR, BATCH_NAME, "{merged_sample}.bam.bai")
    shell: """
        mapfile -t bam_array < {input}
        samtools merge {output.bam} $(echo "${{bam_array[*]}}")
        samtools index {output.bam} > {output.bai}
    """    

################################################################################
rule merge_call_haplotypes:
    ''' Merges aligned deduplicated reads for each replicate. '''
    input:  rules.merge_batch_bams.output.bam
    output: join(JOINT_DIR, BATCH_NAME, "{merged_sample}.g.vcf.gz")
    threads: 2
    params:
        # ploidy = config['variant_calling']['ploidy']
        ploidy = 2
    shell: """
        gatk HaplotypeCaller \
            --reference {REF_FILE} \
            --input {input} \
            --output {output} \
            -ploidy {params.ploidy} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads}
    """

################################################################################
rule  update_gvcf_list:
    input:  gvcf_list
    output: join(JOINT_DIR, "gvcf_lists", f"{BATCH_NAME}_genomicDB_gvcfs.txt")
    run:
        if len(input) > 1:
            update_gvcf_list(str(input), str(output))
        else:
            shell("echo {input}; mv {input} {output} && touch -h {output}")   

################################################################################
rule cumulative_gvcf_list:
    input:  rules.update_gvcf_list.output
    output: join(JOINT_DIR, "gvcf_lists", f"{BATCH_NAME}_cumulative.txt")
    params:
        gvcf_dir = join(JOINT_DIR, "gvcf_lists")
    shell: """
        cd {params.gvcf_dir}
        cat *_gvcfs.txt >> {output}
    """
                     
################################################################################
rule GenomicsDBImport_all:
    input: rules.cumulative_gvcf_list.output
    output: join(JOINT_DIR, "genomicsDB/callset.json")
    threads: 10    
    params:
        db_dir = join(JOINT_DIR, "genomicsDB")
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
    output: join(JOINT_DIR, "vcf_files", f"{BATCH_NAME}_jointGenotype.vcf.gz")
    params:
        db_dir = join(JOINT_DIR, "genomicsDB")
    shell: """
        gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant gendb://{params.db_dir} \
            --output {output} \
    """

