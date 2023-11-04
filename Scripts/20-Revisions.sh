#!/bin/bash

#SBATCH -c 124 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --array=1-10   # Create 1 to the number of fasta files 8 tasks
#SBATCH --output=slurm-%x.%j.out

module load bioinformatics
module load busco/5.3.2

unset OMP_PROC_BIND
unset OMP_PLACES

dir_path=(/nobackup/tmjj24/H_titia_genome/Master/Revisions/busco)
mkdir -p $dir_path
cd $dir_path

# Genome from file
line_num=$(expr $SLURM_ARRAY_TASK_ID)
echo "$line_num"
genome=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master_genome_assembly/Master_HifiAdapt_Filter/genome_files.txt)


# Locations of genomes
input_dir=(/nobackup/tmjj24/Odonata_genomes/All_genomes)

# Busco directory
BUSCO_DIR=(/nobackup/tmjj24/H_titia_genome/BUSCO/Busco_calculations/busco_downloads/)

# busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]

# data must be downloaded from https://busco-data.ezlab.org/v5/data/ and placed in busco download folder
# subfolders of information, lingeages, are needed and files must be unziped using tar -zxvf <file_name> 

# Get file of chromsome lengths
bioawk -c fastx '{ print $name, length($seq), "'"$genome"'" }' $input_dir/$genome.fna > $genome.Chromosome_length.txt

echo "Running BUSCO!"

busco -i $input_dir/$genome.fna \
	-l insecta_odb10 \
	-o ${genome}_insecta_odb10 \
	-m genome \
	-f \
	--offline \
	--cpu 124 \
	--download_path $BUSCO_DIR

busco -i $input_dir/$genome.fna \
	-l arthropoda_odb10 \
	-o ${genome}_arthropoda_odb10 \
	-m genome \
	-f \
	--offline \
	--cpu 124 \
	--download_path $BUSCO_DIR

echo "Arthropoda complete"

echo "Insecta complete"

echo "finished."


# Get full list of busco scores

cp ${genome}_arthropoda_odb10/run_arthropoda_odb10/full_table.tsv ${genome}_arthropoda_odb10_full_table.tsv
cp ${genome}_insecta_odb10/run_insecta_odb10/full_table.tsv ${genome}_insecta_odb10_full_table.tsv

cd ~