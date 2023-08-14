#!/bin/bash

#SBATCH -c 20 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 12:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out


module load python/3.9.9
module load bioinformatics
module load minimap2/2.24

# Setting parameters
dir_path=(/nobackup/tmjj24/H_titia_genome/Master/HiRise/ragtag_alt)
REF=(/nobackup/tmjj24/H_titia_genome/Master/HiRise/HetTit1.0.p.genome.fasta)
QUR_dir=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit)
# QUR_1=(H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.filtered)
QUR_2=(H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.filtered)
mitogenome=(/nobackup/tmjj24/H_titia_genome/Master/mitohifi/mitohifi_p50_MK/final_mitogenome.fasta)

mkdir -p $dir_path
cd $dir_path

# Corrects by breaking up incorrectly assembled contigs based off the primary HiRise scaffolded assembly
ragtag.py correct -t 20 -u -o $dir_path $REF $QUR_dir/$QUR_2.fa


# Scaffolds the alternative assembly by mapping to the  primary HiRise scaffolded assembly
ragtag.py scaffold -t 20 -u -r -o $dir_path $REF ragtag.correct.fasta

