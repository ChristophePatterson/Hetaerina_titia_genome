#!/bin/bash

#SBATCH -c 48 
#SBATCH --mem=100G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:100G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 72:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

module load bioinformatics
module load blobtoolkit2/3.1.6
module load samtools
module load minimap2

unset OMP_PROC_BIND
unset OMP_PLACES

## VARIABLES
BUSCO_DIR=(/nobackup/tmjj24/H_titia_genome/Master/busco)
dir_path_1=(/nobackup/tmjj24/H_titia_genome/Master)
dir_path=($dir_path_1/BlobToolKit/HetTit1.0.p.genome)
Prim_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/HiRise/HetTit1.0.p.genome.fasta)
Prim_asm=(HetTit1.0.p.genome.fasta)

Prim_blob=(HetTit1.0.p.genome)
pb_list=(/nobackup/tmjj24/H_titia_genome/Master/HiFiAdapterFilt/multiple_movies.hifi_reads.filt.fastq.gz)

mkdir -p $dir_path

cp -u $Prim_asm_dp $dir_path
# cp -u /nobackup/tmjj24/H_titia_genome/H_titia_Hifi_reads/multiple_movies.hifi_reads.fastq.gz $dir_path
cp -u /home/tmjj24/scripts/job_scripts/Blob_ToolKit/meta_file_Hetaerina_titia.yaml $dir_path

echo "All files copied"

####################################
## CREATE BLOBTOOLS DIRECTORY ######
####################################

echo "############### Create Blobtoolkit ##############"

cd $dir_path

blobtools create \
    --fasta $Prim_asm \
    ./$Prim_blob

## ADD BUSCO

blobtools add \
	--busco $BUSCO_DIR/HetTit1.0.p_arthropoda_odb10/run_arthropoda_odb10/full_table.tsv \
	./$Prim_blob


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

