#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 08:00:00         # time limit in format dd-hh:mm:ss


# Output directory
dir_location=(/nobackup/tmjj24/H_titia_genome/jupiter/)

mkdir -p $dir_location
cd $dir_location

module load bioinformatics
module load circos/0.69.8
module load samtools/1.15
module load minimap2/2.24

Prim_asm_dp=(/nobackup/tmjj24/H_americana_draft_genome/GCA_022747635.1_ioHetAmer1.0.p_genomic.fna)
Alt_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/Draft_genome_fasta/HetTit1.0.p.genome.fasta)


# /home/tmjj24/apps/JupiterPlot/jupiter name=$prefix ref=$reference fa=$scaffolds

/home/tmjj24/apps/JupiterPlot/jupiter name=HetAmer1.0_vs_HetTit1.0 \
	ref=$Prim_asm_dp \
	fa=$Alt_asm_dp \
	m=100000 ng=85 minBundleSize=50000 maxGap=100000 

cd ~