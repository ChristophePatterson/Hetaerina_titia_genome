#!/bin/bash

#SBATCH -c 4 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 12:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out


# Output directory
dir_location=(/nobackup/tmjj24/H_titia_genome/lastz/)

HetTit=(/nobackup/tmjj24/H_titia_genome/Master/Draft_genome_fasta/HetTit1.0.p.genome.fasta)
HetAme=(/nobackup/tmjj24/H_americana_draft_genome/GCA_022747635.1_ioHetAmer1.0.p_genomic.fna)

mkdir -p $dir_location
cd $dir_location

lastz_32 $HetAme[multiple] $HetTit \
      --notransition --step=20 --nogapped \
      --progress=1 --gfextend --chain \
      --rdotplot=amer_vs_titia.r --format=blastn --ambiguous=iupac > HetAme_vs_HetTit.dat 
