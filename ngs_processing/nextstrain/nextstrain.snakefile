rule all:
    input:
        auspice_json = "auspice/GW.json",

# names of files used in the analysis
seq_file = "/home/internal.idm.ctr/bsiranosian/batch1_processing/03_variant_calls/jointGenotype_2Iter_filtered.vcf.gz"
ref_file = "/home/internal.idm.ctr/bsiranosian/references/d_medinensis_mitochondrion.fasta"
meta_file = "/home/internal.idm.ctr/bsiranosian/batch1_processing/nextstrain/nextstrain_metadata.tsv"
exclude_file = "/home/internal.idm.ctr/bsiranosian/batch1_processing/exclude_samples.txt"
mask_file = ""
drms_file = ""
sites_file = ""
generef_file = ""
genes_file = ""
clades_file = ""
colors_file = ""
auspice_config_file = "/home/internal.idm.ctr/bsiranosian/projects/ben_fork/GWSpatialGenetics/ngs_processing/nextstrain/auspice_config.json"
geo_info_file = "/home/internal.idm.ctr/bsiranosian/batch1_processing/nextstrain/geo_info.tsv"



rule filter:
    input:
        seq = seq_file,
        meta = meta_file,
        exclude = exclude_file
    output:
        "results/filtered.vcf"
    shell: """
        augur filter --sequences {input.seq} \
            --metadata {input.meta} \
            --exclude {input.exclude} \
            --output {output}
        # change deletion asterisk characters to dash? 
        sed -i "s/\\*/\\-/g" {output}
        """
            # --exclude {input.exclude} \

# NOT USED
'''
rule mask:
    input:
        seq = rules.filter.output,
        mask = mask_file
    output:
       "results/masked.vcf.gz"
    shell: """
        augur mask --sequences {input.seq} \
            --mask {input.mask} \
            --output {output}
        """
'''

rule tree:
    input:
        # aln = rules.mask.output,
        aln = rules.filter.output,
        ref = ref_file,
        # sites = sites_file
    output:
        "results/tree_raw.nwk"
    params:
        method = 'iqtree'
    shell: """
        augur tree --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --method {params.method} \
            --output {output}
        """
            # --exclude-sites {input.sites} \

rule refine:
    input:
        tree = rules.tree.output,
        aln = rules.filter.output,
        metadata = meta_file,
        ref = ref_file
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json",
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
        nt_data = "results/nt_muts.json",
        vcf_out = "results/nt_muts.vcf"
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

'''
rule translate:
    input:
        tree = rules.refine.output.tree,
        ref = ref_file,
        gene_ref = generef_file,
        vcf = rules.ancestral.output.vcf_out,
        genes = genes_file
    output:
        aa_data = "results/aa_muts.json",
        vcf_out = "results/translations.vcf",
        vcf_ref = "results/translations_reference.fasta"
    shell: """
        augur translate --tree {input.tree} \
            --vcf-reference {input.ref} \
            --ancestral-sequences {input.vcf} \
            --genes {input.genes} \
            --reference-sequence {input.gene_ref} \
            --output-node-data {output.aa_data} \
            --alignment-output {output.vcf_out} \
            --vcf-reference-output {output.vcf_ref}
        """
'''

rule clades:
    input:
        tree = rules.refine.output.tree,
        # aa_muts = rules.translate.output.aa_data,
        nuc_muts = rules.ancestral.output.nt_data,
        # clades = clades_file
    output:
        clade_data = "results/clades.json"
    shell: """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} \
            --output-node-data {output.clade_data}
        """
            # --mutations {input.nuc_muts} {input.aa_muts} \
            # --clades {input.clades} \

rule traits:
    input:
        tree = rules.refine.output.tree,
        meta = meta_file
    output:
        "results/traits.json"
    params:
        traits = 'location'
    shell: """
        augur traits --tree {input.tree} \
            --metadata {input.meta} \
            --columns {params.traits} \
            --output-node-data {output}
        """

'''
rule seqtraits:
    input:
        align = rules.ancestral.output.vcf_out,
        ref = ref_file,
        # trans_align = rules.translate.output.vcf_out,
        # trans_ref = rules.translate.output.vcf_ref,
        # drms = drms_file
    output:
        drm_data = "results/drms.json"
    params:
        field_to_count = "traits",
        label = "Drug_Resistance"
    shell: """
        augur sequence-traits \
            --ancestral-sequences {input.align} \
            --vcf-reference {input.ref} \
            --count {params.field_to_count} \
            --label {params.label} \
            --output-node-data {output.drm_data}
        """
            # --translations {input.trans_align} \
            # --vcf-translate-reference {input.trans_ref} \
            # --features {input.drms} \
'''

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = meta_file,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output,
        nt_muts = rules.ancestral.output.nt_data,
        # aa_muts = rules.translate.output.aa_data,
        # drms = rules.seqtraits.output.drm_data,
        # color_defs = colors_file,
        auspice_config = auspice_config_file,
        geo_info = geo_info_file,
        # clades = rules.clades.output.clade_data
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
            # --node-data {input.branch_lengths} {input.traits} {input.drms} {input.aa_muts} {input.nt_muts} {input.clades} \
            # --auspice-config {input.auspice_config} \
            # --colors {input.color_defs} \
            # --lat-longs {input.geo_info} \

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

