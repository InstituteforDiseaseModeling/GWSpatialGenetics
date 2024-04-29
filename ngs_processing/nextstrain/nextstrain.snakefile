#######################################################################
# Nextstrain pipeline for GW NGS data
# Designed to take in the input of the ngs variant calling pipeline
# and output a json file that can be visualized with Nextstrain
# This is mostly just the default Nextstrain pipeline, similar to 
#   https://docs.nextstrain.org/en/latest/tutorials/tb_tutorial.html
#   With some modifications and many simplifications. 
#   Many rules, such as translating to amino acid mutations and 
#   reconstruction of ancestral sites have been removed. 
#######################################################################
# Ben Siranosian - Institute for Disease Modeling - June 2021 #
# Jessica Ribado - Institute for Disease Modeling - Nov  2021 and April 2024 
#######################################################################

import sys
from os.path import join, splitext, isfile
from glob import glob
import os.path as path
from os.path import basename


# load config parameters
PROJECT_DIR = config['project_dir']
OUTPUT_DIR = config['output_dir']
INPUT_VCF  = config['diploid_vcf']

if isfile(INPUT_VCF):
    BATCH_NAME = basename(INPUT_VCF).split("_joint")[0]

# validate input
if config['tree_method'] not in ["iqtree", "raxml", "fasttree"]:
    sys.exit("Check config. tree_method must be in [iqtree, raxml, fasttree]")


################################################################################
rule all:
    input:
        expand(join(OUTPUT_DIR, "auspice", "{BATCH_NAME}.json"), BATCH_NAME=BATCH_NAME)


################################################################################
rule haploidize_vcf:
    input:
       diploid_vcf = config['diploid_vcf']
    output:
        haploid_vcf = join(PROJECT_DIR, "{BATCH_NAME}", "{BATCH_NAME}_jointHaploidFilter.vcf.gz"),
        filter_summary = join(PROJECT_DIR, "{BATCH_NAME}", "{BATCH_NAME}_filterSummary.tsv")
    params:
        site_missing  = config['vcf_filter']['site_max_missing'],
        samp_missing  = config['vcf_filter']['samp_max_missing'],
        min_read_depth = config['vcf_filter']['genotype_min_read'],
        het_threshold  = config['vcf_filter']['het_proportion']    
    script: "scripts/vcf_diploid2haploid.R"   


################################################################################
rule create_barcodes:
    input:
        metadata = config['metadata'],
        haploid_vcf = rules.haploidize_vcf.output.haploid_vcf
    output:
        barcode_metadata = join(PROJECT_DIR, "{BATCH_NAME}", "GenomicSamplesMetadata_Database_{BATCH_NAME}.tsv")
    params:
        missing_barcode_max = config['vcf_filter']['missing_barcode_max']    
    script: "scripts/vcf_haploid2cluster.R"


################################################################################
rule colors_to_metadata:
    '''Creates new metadata file with clusters and accompanying color schemes for all metadata columns. '''
    input:
        barcode_metadata = rules.create_barcodes.output.barcode_metadata
    output:
        meta_clust = join(OUTPUT_DIR, "{BATCH_NAME}", "metadata_clusters.tsv"),
        meta_color = join(OUTPUT_DIR, "{BATCH_NAME}", "metadata_colors.tsv"),
        geo_info = join(OUTPUT_DIR, "{BATCH_NAME}", "geo_info.tsv")
    params:
        old_barcodes = join(workflow.basedir, "scripts/original_protocol_barcodes.txt")
    script: "scripts/nextstrain_metadata_check.R"   


################################################################################
rule update_description:
    '''Updates description file to contain filtering parameters. '''
    input:
        readme = config['base_readme'],
        summary = rules.haploidize_vcf.output.filter_summary
    output: join(OUTPUT_DIR, "{BATCH_NAME}", "updated_description.md")
    shell: """
        temp_file=$(mktemp)
        tr "\\t" "\\\|" < {input.summary} > temp_file
        cat {input.readme} temp_file > {output}
        rm temp_file
    """    

################################################################################
rule tree:
    input:
        aln = rules.haploidize_vcf.output.haploid_vcf,
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
        meta = rules.colors_to_metadata.output.meta_clust
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
        ref = config['reference_file']
    output:
        nt_data = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "nt_muts.json"),
        vcf_out = join(OUTPUT_DIR, "{BATCH_NAME}", "results", "nt_muts.vcf")
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
        meta = rules.colors_to_metadata.output.meta_clust
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
        metadata = rules.colors_to_metadata.output.meta_clust,
        colors = rules.colors_to_metadata.output.meta_color,
        geo_info = rules.colors_to_metadata.output.geo_info,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output,
        nt_muts = rules.ancestral.output.nt_data,
        auspice_config = config["auspice_config_file"],
        description = rules.update_description.output
    output:
        auspice_json = join(OUTPUT_DIR, "auspice", "{BATCH_NAME}.json"),
    shell: """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --maintainers "Institute for Disease Modeling <www.idmod.org>" \
            --node-data {input.branch_lengths} {input.traits}  {input.nt_muts} \
            --lat-longs {input.geo_info} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --output {output.auspice_json} \
        """
