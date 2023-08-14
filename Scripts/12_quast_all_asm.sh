#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=200G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:200G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 72:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

module load bioinformatics
module load quast/5.2.0

unset OMP_PROC_BIND
unset OMP_PLACES

## VARIABLES

dir_path=(/nobackup/tmjj24/H_titia_genome/Master/QUAST/H_titia_All_asm_HiRise)

HIFIASM_hap1=(/nobackup/tmjj24/H_titia_genome/Master/Hifiasm/H_titia.asm.bp.hap1.p_ctg.fa)
HIFIASM_hap2=(/nobackup/tmjj24/H_titia_genome/Master/Hifiasm/H_titia.asm.bp.hap2.p_ctg.fa)

PURGE_DUPS_hap1=(/nobackup/tmjj24/H_titia_genome/Master/purge_dups/H_titia.asm.bp.hap1.p_ctg.fa.purged.fa)
PURGE_DUPS_hap2=(/nobackup/tmjj24/H_titia_genome/Master/purge_dups/H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa)

PYLUM_UNCONT_hap1=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit/H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.fa)
PYLUM_UNCONT_hap2=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit/H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.fa)

MITO_UNCONT_hap1=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit/H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.filtered.fa)
MITO_UNCONT_hap2=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit/H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.filtered.fa)

SCAFF_Hap1=(/nobackup/tmjj24/H_titia_genome/Master/RagTag/Hap1/H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.filtered.scaff.fa)
SCAFF_Hap2=(/nobackup/tmjj24/H_titia_genome/Master/RagTag/Hap2/H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.filtered.scaff.fa)

HiRise=(/nobackup/tmjj24/H_titia_genome/Master/HiRise/HetTit1.0.p.genome.fasta)

# REF=(/nobackup/tmjj24/H_americana_draft_genome/GCA_022747635.1_ioHetAmer1.0.p_genomic.fna)
# PACBIO=(/nobackup/tmjj24/H_titia_genome/H_titia_Hifi_reads/multiple_movies.hifi_reads.fastq.gz)

mkdir -p $dir_path

quast.py $Prim_asm_dp \
	$HIFIASM_hap1 \
	$HIFIASM_hap2 \
	$PURGE_DUPS_hap1 \
	$PURGE_DUPS_hap2 \
	$PYLUM_UNCONT_hap1 \
	$PYLUM_UNCONT_hap2 \
	$MITO_UNCONT_hap1 \
	$MITO_UNCONT_hap2 \
	$SCAFF_Hap1 \
	$SCAFF_Hap2 \
	$HiRise \
	-o $dir_path -t 24  \
	--eukaryote --large 
	# --pacbio $PACBIO
