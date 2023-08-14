#!/bin/bash

#SBATCH -c 64 
#SBATCH --mem=100G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:100G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

module load bioinformatics
module load busco/5.3.2

unset OMP_PROC_BIND
unset OMP_PLACES

# Defaults
dir_path=(/nobackup/tmjj24/H_titia_genome/Master/)
mkdir -p $dir_path/blobtools
mkdir -p $dir_path/busco

Prim_asm_dp=($dir_path/purge_dups/H_titia.asm.bp.hap1.p_ctg.fa.purged.fa)
Alt_asm_dp=($dir_path/purge_dups/H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa)

BUSCO_DIR=(/nobackup/tmjj24/H_titia_genome/BUSCO/Busco_calculations/busco_downloads/)

# SPLITING FASTA FILE INTO SMALLER CONTIG FILES RUN BEFORE THE ARRAY

arrary_dir=($dir_path/blobtools/blastn_array)

mkdir -p $arrary_dir
cd $arrary_dir

cd $arrary_dir
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%50==0){file=sprintf("Hap1_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $Prim_asm_dp
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%50==0){file=sprintf("Hap2_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $Alt_asm_dp

# HOW MANY Hap1 files there are. Add -1 to the task array below

ls -lR Hap1_* | wc -l > Hap_1_split_num
ls -lR Hap2_* | wc -l > Hap_2_split_num


cd $dir_path/busco


# data must be downloaded from https://busco-data.ezlab.org/v5/data/ and placed in busco download folder
# subfolders of information, lingeages, are needed and files must be unziped using tar -zxvf <file_name> 

echo "Running!"

busco -i $Prim_asm_dp \
	-l arthropoda_odb10 \
	-o H_titia_purge_Hifiasm_arthropoda_odb10_Purge_Hap1 \
	-m genome \
	-f \
	--offline \
	--cpu 64 \
	--download_path $BUSCO_DIR

echo "Arthropoda complete Hap1"

busco -i $Prim_asm_dp \
	-l insecta_odb10 \
	-o H_titia_purge_Hifiasm_insecta_odb10_Purge_Hap1 \
	-m genome \
	-f \
	--offline \
	--cpu 64 \
	--download_path $BUSCO_DIR

echo "Insecta complete Hap1"


busco -i $Alt_asm_dp \
	-l arthropoda_odb10 \
	-o H_titia_purge_Hifiasm_arthropoda_odb10_Purge_Hap2 \
	-m genome \
	-f \
	--offline \
	--cpu 102 \
	--download_path $BUSCO_DIR

echo "Arthropoda complete Hap2"

busco -i $Alt_asm_dp \
	-l insecta_odb10 \
	-o H_titia_purge_Hifiasm_insecta_odb10_Purge_Hap2 \
	-m genome \
	-f \
	--offline \
	--cpu 102 \
	--download_path $BUSCO_DIR


echo "Insecta complete Hap2"

echo "finished."

cd ~