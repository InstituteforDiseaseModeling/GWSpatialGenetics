# EICC Version specific to EICC server hardware, July 2025.

The Emory Integrated Computational Core at Emory University has updated certain scripts to natively run on their server cluster. Additional scripts, such as wrapper scripts for Qiagen-sequenced data which contain its own metadata generation line, wrapper scripts for Cornell-sequenced data, Cornell old (pre-2024) and current metadata generation scripts, and the like are tailored for EICC's servers and are provided with zero support on an as-is basis per request from Jessica Ribado. These are written by RAA to run on his server cluster and should only be used as a model for one to attempt to generate their own local cluster files for a SLURM batch job submission script UNIX based server. 

The files that may be of importance to all are the sites_only.vcf* files, which are also copied in under input_files to gw_known_amplicon.vcf.gz*. These are from the most recent independent joint calling attempt with novel discovery that did not rely on old known variants that was run on July 14, 2025.

**07/2025: EICC_Version branch contains files updated by Robert A. Arthur at the Emory Integrated Computational Core. They are intended only as reference and to be used on EICC hardware only.**

