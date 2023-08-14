#!/bin/bash

#SBATCH -c 124 
#SBATCH --mem=200G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:200G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 72:00:00         # time limit in format dd-hh:mm:ss

# Commands to execute start here

module load bioinformatics
module load blobtoolkit2/3.1.6
module load samtools
module load minimap2
module load blast
module load busco/5.3.2

unset OMP_PROC_BIND
unset OMP_PLACES

## VARIABLES

dir_path=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit/H_titia_Scaffold_Hap1)
Prim_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/RagTag/Hap1/H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.filtered.scaff.fa)
Alt_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/RagTag/Hap2/H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.filtered.scaff.fa)
Prim_asm=(H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.filtered.scaff.fa)
Alt_asm=(H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.filtered.scaff.fa)

pb_list=(/nobackup/tmjj24/H_titia_genome/HiFiAdapterFilt/multiple_movies.hifi_reads_adaptfilt/multiple_movies.hifi_reads.filt.fastq.gz)

mkdir -p $dir_path

cp -u $Prim_asm_dp $dir_path
cp -u $Alt_asm_dp $dir_path
cp -u /home/tmjj24/scripts/job_scripts/Blob_ToolKit/meta_file_Hetaerina_titia.yaml $dir_path

echo "All files copied"

####################################
## RUNS BUSCO ######################
####################################

dir_path_busco=(/nobackup/tmjj24/H_titia_genome/Master/busco/)
mkdir -p $dir_path_busco
cd $dir_path_busco

BUSCO_DIR=(/nobackup/tmjj24/H_titia_genome/BUSCO/Busco_calculations/busco_downloads/)

# busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]

# data must be downloaded from https://busco-data.ezlab.org/v5/data/ and placed in busco download folder
# subfolders of information, lingeages, are needed and files must be unziped using tar -zxvf <file_name> 

echo "Running!"

busco -i $Prim_asm_dp \
	-l arthropoda_odb10 \
	-o H_titia_purge_Hifiasm_arthropoda_odb10_Purge_Hap1_scaff \
	-m genome \
	--offline \
	--cpu 120 \
	--download_path $BUSCO_DIR

echo "Arthropoda complete Hap1"

busco -i $Prim_asm_dp \
	-l insecta_odb10 \
	-o H_titia_purge_Hifiasm_insecta_odb10_Purge_Hap1_scaff \
	-m genome \
	--offline \
	--cpu 120 \
	--download_path $BUSCO_DIR

echo "Insecta complete Hap1"


busco -i $Alt_asm_dp \
	-l arthropoda_odb10 \
	-o H_titia_purge_Hifiasm_arthropoda_odb10_Purge_Hap2_scaff \
	-m genome \
	--offline \
	--cpu 120 \
	--download_path $BUSCO_DIR

echo "Arthropoda complete Hap2"

busco -i $Alt_asm_dp \
	-l insecta_odb10 \
	-o H_titia_purge_Hifiasm_insecta_odb10_Purge_Hap2_scaff \
	-m genome \
	--offline \
	--cpu 120 \
	--download_path $BUSCO_DIR


echo "Insecta complete Hap2"

echo "finished."

