#!/bin/bash

#SBATCH -c 4 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 12:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out


# Output directory
dir_location=(/nobackup/tmjj24/H_titia_genome/lastz/)

HetTit=(/nobackup/tmjj24/H_titia_genome/Master/Draft_genome_fasta/HetTit1.0.p.genome.fasta)
IolEle=(/nobackup/tmjj24/Odonata_genomes/ncbi-genomes-2023-06-15/GCA_921293095.2_ioIscEleg1.2_genomic.fna)

mkdir -p $dir_location
cd $dir_location

lastz_32 $IolEle[multiple] $HetTit \
      --notransition --step=20 --nogapped \
      --progress=1 --gfextend --chain \
      --rdotplot=IolEle_vs_titia.r --format=blastn --ambiguous=iupac > IolEle_vs_HetTit.dat 
