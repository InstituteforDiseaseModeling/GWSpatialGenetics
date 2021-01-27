# Snakemake pipeline to combine batches of sequencing and do joint variant calling
# uses gatk CombineGVCFs for now, eventually (maybe) move this to GenomicsDBImport
import sys
from os.path import exists, join

# read configfile 
PROJECT_DIR = config["output_directory"]
REF_FILE = config['reference_file']
file_list_f = config["gvcf_list"]

# read the file list
with open(file_list_f, 'r') as f:
    file_list = [a.strip() for a in f.readlines()]

# ensure all files specified actually exist
check_list = [not exists(a) for a in file_list]
if(any(check_list)):
    missing_files = [b for a,b in zip(check_list, file_list) if a]
    print("Missing files: ")
    for a in missing_files: print(a) 
    sys.exit("Some specified files do not exist")




################################################################################
rule all:
    input:
        join(PROJECT_DIR, "combined.g.vcf"),
        join(PROJECT_DIR, "03_variant_calls/jointGenotype_filtered.txt"),


################################################################################
rule combine_gvcf:
    input: 
        file_list
    output: 
        join(PROJECT_DIR, "01_combined/combined.g.vcf")
    params:
        gvcf_string = lambda wildcards: expand("--variant " + "{file}", file=file_list)
    shell: """
        gatk CombineGVCFs \
            --reference {REF_FILE} \
            {params.gvcf_string} \
            --output {output} \
    """


################################################################################
rule joint_genotyping:
    input: rules.combine_gvcf.output
    output: join(PROJECT_DIR, "02_haplocall/joint_genotype.vcf.gz")
    shell: """
        gatk GenotypeGVCFs \
            --reference {REF_FILE} \
            --variant {input} \
            --output {output} \
    """

################################################################################
rule snp_discovery_vcf:
    input:  rules.joint_genotyping.output
    output: temp(join(PROJECT_DIR, "02_haplocall/joint_genotypeFiltered_SNPs.vcf.gz"))
    shell: """
        gatk SelectVariants \
            --reference {REF_FILE} \
            --variant {input} \
            --select-type-to-include SNP \
            --output {output} \
    """

################################################################################
rule filter_discovery_vcf:
    input:  rules.snp_discovery_vcf.output
    output: join(PROJECT_DIR, "02_haplocall/joint_genotypeFiltered.vcf.gz")
    params:
        qd_score = config['variant_filter']['quality_depth'],
        fisher   = config['variant_filter']['fisher_strand'],
        map_qual = config['variant_filter']['mapping_quality'],
        map_root = config['variant_filter']['mapping_rootsq'],
        map_rank = config['variant_filter']['mapping_rank']
    shell: """
        gatk VariantFiltration  \
            --reference {REF_FILE} \
            --variant {input} \
            --output {output} \
            --filter-expression \"QD < {params.qd_score} || FS > {params.fisher} || MQ < {params.map_qual} || MQRankSum < {params.map_root} || ReadPosRankSum < {params.map_rank}\" \
            --filter-name "gatk_hardFilt" \
    """

###############################################################################
rule variant_filter:
    input: rules.filter_discovery_vcf.output
    output: join(PROJECT_DIR, "03_variant_calls/jointGenotype_filtered.vcf.gz")
    params:
        dp_min = config['variant_filter']['read_depth_min'],
        limits = config['variant_filter']['rank_sum']
    shell: """
        gatk VariantFiltration \
            --reference {REF_FILE} \
            --output {output} \
            --variant {input} \
            --filter-expression "ReadPosRankSum <= {params.limits} || ReadPosRankSum >= -{params.limits}" \
            --filter-name "ReadPosRankSum" \
            --filter-expression "ReadQRankSum <= {params.limits} || ReadQRankSum >= -{params.limits}" \
            --filter-name "ReadQRankSum" \
            --filter-expression "MQRankSum <= {params.limits} || MQRankSum >= -{params.limits}" \
            --filter-name "MQRankSum" \
            --filter-expression "ClippingRankSum <= {params.limits} || ClippingRankSum >= -{params.limits}" \
            --filter-name "ClippingRankSum" \
            --filter-expression "DP > {params.dp_min}" \
            --filter-name "DPMin"
    """

###############################################################################
rule vcf2txt:
    input: rules.variant_filter.output
    output: join(PROJECT_DIR, "03_variant_calls/jointGenotype_filtered.txt")
        #join(PROJECT_DIR, "03_variant_calls/bcftools/all_merged.txt")
    shell: """
        gatk VariantsToTable \
            --variant {input} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -GF DP  \
            -raw \
            --output {output}
    """
