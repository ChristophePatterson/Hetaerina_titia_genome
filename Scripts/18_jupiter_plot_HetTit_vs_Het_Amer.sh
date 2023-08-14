#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=40G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:40G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 02:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Output directory
dir_location=(/nobackup/tmjj24/H_titia_genome/jupiter/)

mkdir -p $dir_location
cd $dir_location


module load bioinformatics
module load circos/0.69.8
module load samtools/1.17
module load minimap2/2.24

HetTit=(/nobackup/tmjj24/H_titia_genome/Master/Draft_genome_fasta/HetTit1.0.p.genome.fasta)
HetAmer=(/nobackup/tmjj24/H_americana_draft_genome/GCA_022747635.1_ioHetAmer1.0.p_genomic.fna)
CalSpe=(/nobackup/tmjj24/Odonata_genomes/Calopteryx_splendens/ncbi_dataset/data/GCA_002093875.1/GCA_002093875.1_Calsple1.0_genomic.fna)
IolEle=(/nobackup/tmjj24/Odonata_genomes/ncbi-genomes-2023-06-15/GCA_921293095.2_ioIscEleg1.2_genomic.fna)

output_name=(HetTit1.0_vs_HetAmer1.0)

# /home/tmjj24/apps/JupiterPlot/jupiter name=$prefix ref=$reference fa=$scaffolds

/home/tmjj24/apps/JupiterPlot/jupiter name=$output_name.v8\
	ref=$HetTit \
	fa=$HetAmer \
	m=4000000 ng=95 t=24 # v8
#	t=24 m=100000 ng=80 minBundleSize=2500 maxGap=1000000 labels=both # v7
#	t=24 m=100000 ng=85 minBundleSize=1000 maxGap=1000000 # v4
#	t=24 m=100000 ng=85 minBundleSize=5000 maxGap=1000000 # v3
#	t=24 m=100000 ng=85 minBundleSize=50000 maxGap=100000 
# sam=/nobackup/tmjj24/H_titia_genome/jupiter/IolEle1.2_vs_HetTit1.0_v3-agp.sam \

bioawk -c fastx '{ print $name, length($seq), "HetTit" }' $HetTit > $output_name.Chromosome_length.txt
bioawk -c fastx '{ print $name, length($seq), "HetAmer" }' $HetAmer >> $output_name.Chromosome_length.txt

cd ~