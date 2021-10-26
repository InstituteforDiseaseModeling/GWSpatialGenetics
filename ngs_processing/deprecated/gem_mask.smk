#  Masking discovery to match James's Cotton guinea worm processing pipeline.
# Unfortuntely, GEM mappability is not in Conda for easy install. Followed the following tutorial to download and run.
# https://evodify.com/gem-mappability/
# Got new errors that had to do with the program itself. 
# Instead, I am using another mappability tool available in conda. 
# Only keeping the most unique regions of the genome to use as known variants in the "edited" mappability file.

################################################################################
rule ref_index:
    input:
        ref = REF_FILE,
    output: join(REF_DIR, "ref_mappability/index.rev.lf.drp.sbl")
    params:
        index_dir = join(REF_DIR, "ref_mappability")
    threads: 10
    shell:  """
        rm -rf {params.index_dir}
        genmap index -F {input.ref} -I {params.index_dir}
    """

################################################################################
rule ref_mappability:
    input:
        index = rules.ref_index.output
    output: join(REF_DIR, "ref_mappability",  REF_NAME.split(".")[0] + ".genmap.bed")
    params:
        index_dir = join(REF_DIR, "ref_mappability"),
        kmer_size = config['gem_masking']['kmer_size'],
        mismatch = config['gem_masking']['mistmatch']
    threads: 10
    shell:  """
        genmap map -T {threads} \
            -I {params.index_dir} -O {params.index_dir} \
            -K {params.kmer_size} --bed
    """

################################################################################
rule ref_mappability_edit:
    ''' Subset the high mappability (i.e. most unique) regions of the genome to use for recalibration. 
    Addtionally removes additional information from the contig identifier to match GATK produced VCFs. '''
     input:  join(REF_DIR, "ref_mappability",  REF_NAME.split(".")[0] + ".genmap.bed")
    output: join(REF_DIR, "ref_mappability",  REF_NAME.split(".")[0] + ".genmapEdit.bed")
    run:
        with open(output[0], "w") as fout:
            with open(input[0], "r") as bed:
                for line in bed:
                    chr, start, end, name, score = line.strip().split('\t')
                    if float(score) > 0.9:
                        chr_match = chr.split(" ")[0]
                        new_line = '\t'.join([chr_match, start, end, name, score])
                        fout.write(new_line+'\n')

################################################################################
rule mappability_to_gatk:
    input: rules.ref_mappability_edit.output
    output: join(REF_DIR, "ref_mappability",  REF_NAME.split(".")[0] + ".genmapEdit.bed.idx")
    shell: """
        gatk IndexFeatureFile -I {input}
    """
