
################################################################################
rule build_ref_index:
    ''' Create the index neccessary for BWA alignment. '''
    input: REF_FILE
    output: expand(join(REF_DIR, REF_NAME + ".{bwa_ext}"),  bwa_ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
	shell: """
		bwa index {input}
	"""

################################################################################
rule create_ref_dict:
    ''' Create sequence dictionary for GATK variant calling. '''
    input:  REF_FILE
    output: join(REF_DIR, REF_NAME.split(".")[0] + ".dict")
    shell: """
        gatk CreateSequenceDictionary --REFERENCE {input}
    """
