#######################################################################
# Nextstrain pipeline for GW NGS data
# desigend to take in the input of the ngs variant calling pipeline
# and output a json file that can be visualized with nextstrain
# This is mostly just the default nextstrain pipeline, similar to 
#   https://docs.nextstrain.org/en/latest/tutorials/tb_tutorial.html
#   With some modifications and many simplifications. 
#   Many rules, such as translating to amino acid mutations and 
#   reconstruction of ancestral sites have been removed. 
#######################################################################
# Ben Siranosian - Institute for Disease Modeling - June 2021 #########
# Jessica Ribado - Institute for Disease Modeling - Nov  2021 #########
#######################################################################

import sys
from os.path import join, splitext, isfile
from glob import glob

# define functions
def input_vcf():
    if isfile(INPUT_VCF):
        return INPUT_VCF
    else:
        return join(config['diploid_dir'], "{BATCH_NAME}_jointGenotypeFiltered.vcf.gz")

# load config parameters
OUTPUT_DIR = config['output_dir']
INPUT_VCF  = config['diploid_vcf']

if isfile(INPUT_VCF):
    BATCH_NAME = basename(INPUT_VCF).split("_joint")[0]
else:
    BATCH_NAME, = glob_wildcards(join(config['diploid_dir'], "{BATCH_NAME}_jointGenotypeFiltered.vcf.gz")) 
    print("VCF prefixes identified:", BATCH_NAME)
    if len(BATCH_NAME) < 1:
        sys.exit("No VCF files identified with correct suffix in provided directory. Check again.")


# validate input
if config['tree_method'] not in ["iqtree", "raxml", "fasttree"]:
    sys.exit("Check config. tree_method must be in [iqtree, raxml, fasttree]")


################################################################################
rule all:
    input:
        expand(join(OUTPUT_DIR, "auspice", "GW_{BATCH_NAME}.json"), BATCH_NAME=BATCH_NAME[1])


################################################################################
rule haploidize_vcf:
    ''' Converts diploid VCF to haploid and filter low quality samples and sites.'''
    input: input_vcf()
    output: join(OUTPUT_DIR, "filtered_vcf", "{BATCH_NAME}_jointHaploidFilter.vcf.gz") 
    params:
        site_missing  = config['vcf_filter']['site_max_missing'],
        samp_missing  = config['vcf_filter']['samp_max_missing'],
        min_read_depth = config['vcf_filter']['genotype_min_read'],
        het_threshold  = config['vcf_filter']['het_proportion']
    script: "scripts/vcf_diploid2haploid.R"


################################################################################
rule update_meta:
    '''Creates new metadata file with clusters and accompanying color schemes for all metadata columns. '''
    input:
        metadata = config['metadata'],
        haploid_vcf = rules.haploidize_vcf.output
    output:
        meta_clust = join(OUTPUT_DIR, "{BATCH_NAME}", "metadata_clusters.tsv"),
        meta_color = join(OUTPUT_DIR, "{BATCH_NAME}", "metadata_colors.tsv"),
        geo_info = join(OUTPUT_DIR, "{BATCH_NAME}", "geo_info.tsv")
    script: "scripts/vcf_haploid2cluster.R"    


################################################################################
rule tree:
    input:
        aln = rules.haploidize_vcf.output,
        ref = config['reference_file'],
    output:
        aln = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "filtered_recode.vcf.gz"),
        tree  = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "tree.nwk"),
        sites = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "informative_sites.fasta")
    params:
        method = config['tree_method']
    shell: """
        # change deletion asterisk characters to dash
        zcat {input.aln} | sed "s/\\*/\\-/g" | bgzip > {output.aln}
        augur tree --alignment {output.aln} \
            --vcf-reference {input.ref} \
            --method {params.method} \
            --output {output.tree}
       """

################################################################################
rule refine:
    input:
        tree = rules.tree.output.tree,
        aln = rules.tree.output.aln,
        ref = config['reference_file'],
        meta = rules.update_meta.output.meta_clust
    output:
        tree = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "tree_refine.nwk"),
        node_data = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "branch_lengths.json")
    params:
        root = 'min_dev',
        coal = 'opt'
    shell: """
        augur refine --tree {input.tree} \
            --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --metadata {input.meta} \
            --timetree \
            --root {params.root} \
            --coalescent {params.coal} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """
        
################################################################################
rule ancestral:
    input:
        tree = rules.refine.output.tree,
        aln = rules.tree.output.aln,
        ref = config['reference_file'],
    output:
        nt_data = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "nt_muts.json"),
        vcf_out = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "nt_muts.vcf"),
    params:
        inference = "joint"
    shell: """
        augur ancestral --tree {input.tree} \
            --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --inference {params.inference} \
            --output-node-data {output.nt_data} \
            --output-vcf {output.vcf_out}
        """

################################################################################
rule traits:
    input:
        tree = rules.refine.output.tree,
        meta = rules.update_meta.output.meta_clust
    output:
        join(OUTPUT_DIR, "{BATCH_NAME}", "results", "traits.json")
    params:
        traits = 'Country'
    shell: """
        augur traits --tree {input.tree} \
            --metadata {input.meta} \
            --columns {params.traits} \
            --output-node-data {output}
        """

################################################################################
rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.update_meta.output.meta_clust,
        colors = rules.update_meta.output.meta_color,
        geo_info = rules.update_meta.output.geo_info,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output,
        nt_muts = rules.ancestral.output.nt_data,
        auspice_config = config["auspice_config_file"]
    output:
        auspice_json = join(OUTPUT_DIR, "auspice", "GW_{BATCH_NAME}.json"),
    shell: """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --node-data {input.branch_lengths} {input.traits}  {input.nt_muts} \
            --lat-longs {input.geo_info} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
        """
