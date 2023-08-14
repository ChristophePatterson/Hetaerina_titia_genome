# Run using an interactive slurm job

srun -t 8:00:00 -c 24 --mem=50G --gres=tmp:50G --pty bash

# Clear all files from directory other than your fasta data


# Open a singularity container with bash
singularity run --bind /nobackup/tmjj24/:/nobackup/tmjj24/ --home /nobackup/tmjj24/H_titia_genome/Master_HifiFilt/mitohifi/mitohifi_p50_MK /home/tmjj24/apps/MitoHiFi/mitohifi_2.2_cv1.sif bash


# Open a singularity container with bash
singularity run --bind /nobackup/tmjj24/:/nobackup/tmjj24/ --home /nobackup/tmjj24/H_titia_genome/Master/mitohifi/mitohifi_p50_MK /home/tmjj24/apps/MitoHiFi/mitohifi_2.2_cv1.sif bash

findMitoReference.py --species "Hetaerina titia" --email christophe.patterson@durham.ac.uk --outfolder data --min_length 16000

# Run the mitohifi pipeline
mitohifi.py -r /nobackup/tmjj24/H_titia_genome/Master/HiFiAdapterFilt/multiple_movies.hifi_reads.filt.fastq.gz -f data/MK722304.1.fasta -g data/MK722304.1.gb -t 24 -o 5

