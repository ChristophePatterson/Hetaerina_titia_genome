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

dir_path=(/nobackup/tmjj24/H_titia_genome/Master/busco)
mkdir -p $dir_path
cd $dir_path


HiRise_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/Draft_genome_fasta/HetTit1.0.a.genome.fasta)
BUSCO_DIR=(/nobackup/tmjj24/H_titia_genome/BUSCO/Busco_calculations/busco_downloads/)

# busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]

# data must be downloaded from https://busco-data.ezlab.org/v5/data/ and placed in busco download folder
# subfolders of information, lingeages, are needed and files must be unziped using tar -zxvf <file_name> 

echo "Running!"

busco -i $HiRise_asm_dp \
	-l arthropoda_odb10 \
	-o HetTit1.0.p_arthropoda_odb10 \
	-m genome \
	-f \
	--offline \
	--cpu 64 \
	--download_path $BUSCO_DIR

echo "Arthropoda complete Hap1"

busco -i $HiRise_asm_dp \
	-l insecta_odb10 \
	-o HetTit1.0.p_insecta_odb10 \
	-m genome \
	-f \
	--offline \
	--cpu 64 \
	--download_path $BUSCO_DIR

echo "Insecta complete Hap1"

echo "finished."

cd ~