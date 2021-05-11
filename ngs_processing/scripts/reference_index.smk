from os.path import splitext

################################################################################
rule build_ref_index:
    ''' Create the index neccessary for BWA alignment. '''
    input: REF_FILE
    output: 
        expand(REF_FILE + ".{bwa_ext}",  bwa_ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
    shell: """
        bwa index {input}
    """

################################################################################
rule build_ref_faidx:
    ''' Create the fasta index for alignment. '''
    input: REF_FILE
    output: REF_FILE + ".fai"
    shell: """
        samtools faidx {input}
    """

################################################################################
rule create_ref_dict:
    ''' Create sequence dictionary for GATK variant calling. '''
    input:  REF_FILE
    output: splitext(REF_FILE)[0] + ".dict"
    shell: """
        gatk CreateSequenceDictionary --REFERENCE {input}
    """
