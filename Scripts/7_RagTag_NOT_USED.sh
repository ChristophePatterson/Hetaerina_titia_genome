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
dir_path=(/nobackup/tmjj24/H_titia_genome/Master/RagTag)
REF=(/nobackup/tmjj24/H_americana_draft_genome/GCA_022747635.1_ioHetAmer1.0.p_genomic.fna)
QUR_dir=(/nobackup/tmjj24/H_titia_genome/Master/BlobToolKit)
QUR_1=(H_titia.asm.bp.hap1.p_ctg.fa.purged.fa.uncont.filtered.filtered)
QUR_2=(H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa.uncont.filtered.filtered)
mitogenome=(/nobackup/tmjj24/H_titia_genome/Master/mitohifi/mitohifi_p50_MK/final_mitogenome.fasta)

mkdir -p $dir_path
cd $dir_path

mkdir -p Hap1
mkdir -p Hap2

# Builds scaffold of haplotype 1
ragtag.py scaffold -t 20 -u -r -o $dir_path/Hap1 $REF $QUR_dir/$QUR_1.fa

chmod u+x $dir_path/Hap1/ragtag.scaffold.fasta

cat $dir_path/Hap1/ragtag.scaffold.fasta > $dir_path/Hap1/$QUR_1.scaff.fa

# Builds scaffold of haplotype 2
ragtag.py scaffold -t 20 -u -r -o $dir_path/Hap2 $REF $QUR_dir/$QUR_2.fa

chmod u+x $dir_path/Hap2/ragtag.scaffold.fasta

cat $dir_path/Hap2/ragtag.scaffold.fasta > $dir_path/Hap2/$QUR_2.scaff.fa

# Adds mitocondrial draft genome to haplotype 1

# cat $mitogenome >> $dir_path/Hap1/$QUR_1.scaff.fa

