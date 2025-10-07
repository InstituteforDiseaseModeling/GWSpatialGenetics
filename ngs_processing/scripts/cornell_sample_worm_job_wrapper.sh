#!/bin/bash

#SBATCH -p nodes
#SBATCH -D .
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 60
#SBATCH -J GWEP_worm_analysis_batch_Apr032025
#SBATCH -o worm_batch_Apr032025.out
#SBATCH -e worm_batch_Apr032025.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rarthu3@emory.edu

pwd
echo "Worm analysis pipeline; `date`"

source activate worm

#New version for Skynet implementation, 7/29/24 RAA
#added switches around main and uploading job submission portion to make it easier
#to use what is needed at the time.

#set data dir to correct since we have to work around conda install from a working dir
#that is not the same

datadir="/home/rarthu3/projects/worm/batch_Apr032025"
workdir=$(pwd)
#workdir is "worm/test_install" for skynet
echo -e "Data directory where .fastqs and metadata manifest are located is:\n"$datadir"\nWorking directory where conda install is pre-done is:\n"$workdir"\n"

#make sure directories are unlocked...
#snakemake -s ~/git/worm_software/gw_processing.snakefile --configfile config.txt --cores 10 --jobs 10 --unlock

#echo "Generating manifest via new snakemake tool"
#If Cornell - use my own script to make the manifest/metadata instead
#use Jessica's new tool to make manifest
#snakemake -s /home/rarthu3/git/worm/scripts/generate_manifest.smk --config fastq_dir='/home/rarthu3/projects/worm/batch_Apr032025/' output_file='/home/rarthu3/projects/worm/batch_Apr032025/metadata.tsv' sample_key='sample-key.csv'

#initiate actual analytic pipeline
callmain=true
if [ "$callmain" = true ]; then
	echo -e "beginning main pipeline\nNote that the downloading and installation of requisite conda env may take up to 60-90 minutes!"
	snakemake -s ~/git/worm/gw_processing.snakefile --configfile $datadir/config.txt --cores 60 --jobs 60 --nolock --use-conda --keep-going --latency-wait 30 -R check_samples
	echo "Main analysis done at `date`"
else
	echo "Not running main analysis"
fi

conda deactivate
cd $datadir
pwd

uploading=true
if [ "$uploading" = true ]; then
	echo "Submitting tar compression/uploading script"
	sbatch outcopy.sh
else
	echo "Not uploading data"
fi

echo "Script complete at `date`"
