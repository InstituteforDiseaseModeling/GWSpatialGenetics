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
#######################################################################
import sys
from os.path import join, splitext

# load config parameters
tree_method = config["tree_method"]
outdir = config["outdir"]
seq_file = config["seq_file"]
ref_file = config["ref_file"]
meta_file = config["meta_file"]
exclude_file = config["exclude_file"]
auspice_config_file = config["auspice_config_file"]
geo_info_file = config["geo_info_file"]

# validate input
if tree_method not in ["iqtree", "raxml", "fasttree"]:
    sys.exit("Check config. tree_method must be in [iqtree, raxml, fasttree]")

rule all:
    input:
        auspice_json = join(outdir, "auspice/GW.json"),
        new_metadata = splitext(meta_file)[0] + "_new_clusters.tsv"

rule filter:
    input:
        seq = seq_file,
        meta = meta_file,
        exclude = exclude_file
    output:
        join(outdir, "results/filtered.vcf")
    shell: """
        augur filter --sequences {input.seq} \
            --metadata {input.meta} \
            --exclude {input.exclude} \
            --output {output}
        # change deletion asterisk characters to dash
        sed -i "s/\\*/\\-/g" {output}
        """

rule tree:
    input:
        aln = rules.filter.output,
        ref = ref_file,
    output:
        tree = join(outdir, "results/tree_raw.nwk"),
        sites = join(outdir, "results/informative_sites.fasta")
    params:
        method = config['tree_method']
    shell: """
        augur tree --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --method {params.method} \
            --output {output.tree}
        """

# decompose the input VCF (changes biallelic variants)
rule vcf_decompose:
    input: rules.filter.output
    output: join(outdir, "results/filtered_decomposed.vcf")
    shell: """
        vt decompose {input} > {output}
    """

# extract genotype matrix from the decomposed vcf
rule extract_genotype:
    input: rules.vcf_decompose.output
    output: join(outdir, "results/filtered_decomposed_genotype.GT.FORMAT")
    params: 
        outstring = join(outdir, "results/filtered_decomposed_genotype")
    shell: """
        vcftools --vcf {input} --extract-FORMAT-info GT --out {params.outstring}
    """

rule add_new_clusters:
    input:
        genotype = rules.extract_genotype.output,
        metadata = meta_file
    output:
        new_metadata = splitext(meta_file)[0] + "_new_clusters.tsv"
    params:
        hclust_height = config['hclust_height']
    script: "scripts/add_new_clusters.R"


# alternative version that uses cd_hit clustering on the informative sites
# no longer used
# rule cd_hit:
#     input:
#         join(outdir, "results/informative_sites.fasta")
#     output:
#         join(outdir, "results/cd_hit_clustering.clstr")
#     params:
#         out_arg = join(outdir, "results/cd_hit_clustering")
#     shell: """
#         cd-hit-est -i {input} -o  {params.out_arg} -g 1 -sc 1
#     """

# rule add_new_clusters_cd_hit:
#     input:
#         clusters = join(outdir, "results/cd_hit_clustering.clstr"),
#         metadata = meta_file
#     output:
#         new_metadata = splitext(meta_file)[0] + "_new_clusters.tsv"
#     script: "scripts/add_new_clusters.R"


rule refine:
    input:
        tree = rules.tree.output.tree,
        aln = rules.filter.output,
        metadata = meta_file,
        ref = ref_file
    output:
        tree = join(outdir, "results/tree.nwk"),
        node_data = join(outdir, "results/branch_lengths.json"),
    params:
        root = 'min_dev',
        coal = 'opt'
    shell: """
        augur refine --tree {input.tree} \
            --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --metadata {input.metadata} \
            --timetree \
            --root {params.root} \
            --coalescent {params.coal} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.filter.output,
        ref = ref_file
    output:
        nt_data = join(outdir, "results/nt_muts.json"),
        vcf_out = join(outdir, "results/nt_muts.vcf"),
    params:
        inference = "joint"
    shell: """
        augur ancestral --tree {input.tree} \
            --alignment {input.alignment} \
            --vcf-reference {input.ref} \
            --inference {params.inference} \
            --output-node-data {output.nt_data} \
            --output-vcf {output.vcf_out}
        """

rule traits:
    input:
        tree = rules.refine.output.tree,
        meta = meta_file,
    output:
        join(outdir, "results/traits.json")
    params:
        traits = 'location'
    shell: """
        augur traits --tree {input.tree} \
            --metadata {input.meta} \
            --columns {params.traits} \
            --output-node-data {output}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_new_clusters.output.new_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output,
        nt_muts = rules.ancestral.output.nt_data,
        auspice_config = auspice_config_file,
        geo_info = geo_info_file,
    output:
        auspice_json = rules.all.input.auspice_json,
    shell: """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits}  {input.nt_muts} \
            --lat-longs {input.geo_info} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
        """
