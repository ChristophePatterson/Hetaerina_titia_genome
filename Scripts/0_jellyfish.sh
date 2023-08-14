#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=200G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:200G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 12:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

module load bioinformatics
module load jellyfish/2.3.0

# Commands to be run:

dir_path=(/nobackup/tmjj24/H_titia_genome/Master/jellyfish)
mkdir -p $dir_path

# run jellyfish count command to calculate number of k-mers across sequence data

pb_fastq=(/nobackup/tmjj24/H_titia_genome/Master/HiFiAdapterFilt/multiple_movies.hifi_reads.filt.fastq.gz)

zcat $pb_fastq > /nobackup/tmjj24/H_titia_genome/jellyfish/multiple_movies.hifi_reads.fastq

jellyfish count -t 24 -C -m 21 -s 5G -o $dir_path/k21_reads.jf /nobackup/tmjj24/H_titia_genome/HiFiAdapterFilt/multiple_movies.hifi_reads.fastq

jellyfish histo -t 24 $dir_path/k21_reads.jf > /nobackup/tmjj24/H_titia_genome/Master/jellyfish/k21_reads.histo
