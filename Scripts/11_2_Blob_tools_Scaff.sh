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

dir_path=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit/H_titia_Scaffold)
Prim_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/RagTag/Hap1/H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.filtered.scaff.fa)
Alt_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/RagTag/Hap2/H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.filtered.scaff.fa)
Prim_asm=(H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.filtered.scaff.fa)
Alt_asm=(H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.filtered.scaff.fa)

Prim_blob=(H_titia_Hap1_scaff)
Alt_blob=(H_titia_Hap2_scaff)
pb_list=(/nobackup/tmjj24/H_titia_genome/Master/HiFiAdapterFilt/multiple_movies.hifi_reads.filt.fastq.gz)

mkdir -p $dir_path

cp -u $Prim_asm_dp $dir_path
cp -u $Alt_asm_dp $dir_path
cp -u /home/tmjj24/scripts/job_scripts/Blob_ToolKit/meta_file_Hetaerina_titia.yaml $dir_path

echo "All files copied"

####################################
## CREATE BLOBTOOLS DIRECTORY ######
####################################

cd $dir_path

blobtools create \
    --fasta $Prim_asm \
    ./$Prim_blob

blobtools create \
    --fasta $Alt_asm \
    ./$Alt_blob

## ADD BUSCO

blobtools add \
	--busco /nobackup/tmjj24/H_titia_genome/Master/busco/H_titia_purge_Hifiasm_arthropoda_odb10_Purge_Hap1_scaff/run_arthropoda_odb10/full_table.tsv \
	./$Prim_blob

blobtools add \
	--busco /nobackup/tmjj24/H_titia_genome/Master/busco/H_titia_purge_Hifiasm_insecta_odb10_Purge_Hap1_scaff/run_insecta_odb10/full_table.tsv \
	./$Prim_blob

# blobtools add \
# 	--busco /nobackup/tmjj24/H_titia_genome/Master/busco/H_titia_purge_Hifiasm_arthropoda_odb10_Purge_Hap2_scaff/run_arthropoda_odb10/full_table.tsv \
# 	./$Alt_blob

# blobtools add \
# 	--busco /nobackup/tmjj24/H_titia_genome/Master/busco/H_titia_purge_Hifiasm_insecta_odb10_Purge_Hap2_scaff/run_insecta_odb10/full_table.tsv \
# 	./$Alt_blob


##################
### COVERAGE #####
##################

echo "############### Coverage ##############"

minimap2 -ax map-pb \
         -t 124 $Prim_asm \
         $pb_list \
| samtools sort -@16 -O BAM -o $Prim_asm.reads.bam -

blobtools add \
    --cov $Prim_asm.reads.bam \
    ./$Prim_blob

minimap2 -ax map-pb \
         -t 124 $Alt_asm \
         $pb_list \
| samtools sort -@16 -O BAM -o $Alt_asm.reads.bam -

blobtools add \
    --cov $Alt_asm.reads.bam \
    ./$Alt_blob

