#!/bin/bash

# Request resources (per task):
#SBATCH -c 2           # 1 CPU core
#SBATCH --mem=5G       # 1 GB RAM
#SBATCH --time=72:00:00   # 6 hours (hours:minutes:seconds)
#SBATCH --output=slurm-%x.%j.%A.%a.out

# Run on the shared queue
#SBATCH -p shared


# Specify the tasks to run:
#SBATCH --array=0-82   # Create 0 to the number of fasta files -1 tasks

# Each separate task can be identified based on the SLURM_ARRAY_TASK_ID
# environment variable:

echo "I am task number $SLURM_ARRAY_TASK_ID"

module load bioinformatics
module load blobtoolkit2/3.1.6
module load blast/2.13.0

dir_path=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit/blastn_array)

cd $dir_path

FILES=(Hap2_*)

FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "This job should use $FILE"

head "${FILE}" --bytes 100
echo "...."
tail "${FILE}" -c 100 

echo ""

echo "Starting Blastn"

export BLASTDB='/nobackup/tmjj24/H_titia_genome/BlobToolKit/nt'


blastn -db nt \
       -query ${FILE} \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 2 \
       -out blastn.${FILE}.out

echo "finished blastn"
   

# Run this code using sbatch -o sbatch_%A_%a.out /home/tmjj24/scripts/job_scripts/Blob_ToolKit/Blob_ToolKit_man_blastn_array.sh

# Merge all files together (untested 10-10-2022) blastn.Hap2_*.out > blastn.all.Hap2.out


