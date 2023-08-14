#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 06:00:00         # time limit in format dd-hh:mm:ss

# Commands to execute start here

module load bioinformatics
module load quast/5.2.0

unset OMP_PROC_BIND
unset OMP_PLACES

## VARIABLES

dir_path=(/nobackup/tmjj24/H_titia_genome/Master/QUAST/HetTit1.0)
Prim_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/Draft_genome_fasta/HetTit1.0.p.genome.fasta)
Alt_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/Draft_genome_fasta/HetTit1.0.a.genome.fasta)

REF=(/nobackup/tmjj24/H_americana_draft_genome/GCA_022747635.1_ioHetAmer1.0.p_genomic.fna)
# FEATURES=(/nobackup/tmjj24/H_titia_genome/Master/eggnog-mapper/eggnog-runs/Hap1_run_1/H_titia_Hap1_eggnog_genome_mmseqs.emapper.genepred.gff)
PACBIO=(/nobackup/tmjj24/H_titia_genome/Master/HiFiAdapterFilt/multiple_movies.hifi_reads.filt.fastq.gz)

mkdir -p $dir_path

quast.py $Prim_asm_dp $Alt_asm_dp -o $dir_path -t 24 --split-scaffolds --eukaryote --large # -r $REF --pacbio $PACBIO # -g $FEATURES
