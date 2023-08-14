#!/bin/bash

# Request resources (per task):
#SBATCH -c 64           # 1 CPU core
#SBATCH --mem=200G       # 1 GB RAM
#SBATCH --time=72:00:00   # 6 hours (hours:minutes:seconds)
#SBATCH --output=slurm-%x.%j.out

# Run on the shared queue
#SBATCH -p shared

module load bioinformatics
module load bamtools/2.5.1
module load blast/2.13.0

# Directory and file name for output
dir_path=(/nobackup/tmjj24/H_titia_genome/Master)

OUT_FOLDER=($dir_path/HiFiAdapterFilt)

mkdir -p $OUT_FOLDER

############################
######## HIFIASM ###########
############################

# Filter PacBio seq from HifAdaptFilt
pb_list=($dir_path/HiFiAdapterFilt/multiple_movies.hifi_reads.filt.fastq.gz)

module load hifiasm/0.16.1

mkdir -p $dir_path/Hifiasm

cd $dir_path/Hifiasm

hifiasm -o H_titia.asm -t 64 $pb_list

awk '/^S/{print ">"$2;print $3}' H_titia.asm.bp.hap1.p_ctg.gfa > H_titia.asm.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' H_titia.asm.bp.hap2.p_ctg.gfa > H_titia.asm.bp.hap2.p_ctg.fa


############################
######## purge_dups ########
############################

module load purge_dups
module load minimap2

echo "Set cd and assign asm and PacBio reads"

mkdir -p $dir_path/purge_dups
cd $dir_path/purge_dups

dir_purge=($dir_path/purge_dups)
hap1_asm=($dir_path/Hifiasm/H_titia.asm.bp.hap1.p_ctg.fa)
hap1_asm_name=(H_titia.asm.bp.hap1.p_ctg.fa)
hap2_asm=($dir_path/Hifiasm/H_titia.asm.bp.hap2.p_ctg.fa)
hap2_asm_name=(H_titia.asm.bp.hap2.p_ctg.fa)

echo "Varibles made"

### STEP 1 - Hap 1 Purge
echo "Run minimap2 to align pacbio data and generate paf files on Hap1"

minimap2 -map-hifi $hap1_asm $pb_list -t 124 | gzip -c - > $hap1_asm_name.paf.gz

# produces PB.base.cov and PB.stat files

pbcstat *.paf.gz 
calcuts PB.stat > cutoffs 2>calcults.log

echo "Split an assembly and do a self-self alignment"

split_fa $hap1_asm > $hap1_asm_name.split
minimap2 -xasm5 -DP $hap1_asm_name.split $hap1_asm_name.split -t 100 | gzip -c - > $hap1_asm_name.split.self.paf.gz

purge_dups -2 -T cutoffs -c PB.base.cov $hap1_asm_name.split.self.paf.gz > dups.bed 2> purge_dups.log

echo "Get purged primary and haplotig sequences from draft assembly"

# get_seqs dups.bed $pri_asm -p H_titia.asm_purged_all
get_seqs -e dups.bed $hap1_asm -p $hap1_asm_name

echo "Merge hapoloypes and alternative assembly"

cat $hap2_asm > $hap2_asm_name.purged.fa
cat $hap1_asm_name.hap.fa >> $hap2_asm_name.purged.fa

echo "DONE HAP 1 YAY!"

echo "Running duplication on purged hap 2!"

### STEP 2 - Hap 2 Purge
###
echo "Run minimap2 to align pacbio data and generate paf files on Hap2"

minimap2 -map-hifi $hap2_asm_name.purged.fa $pb_list -t 124 | gzip -c - > $hap2_asm_name.purged.paf.gz

# produces PB.base.cov and PB.stat files

pbcstat $hap2_asm_name.purged.paf.gz 
calcuts PB.stat > cutoffs 2>calcults_Hap2.log

echo "Split an assembly and do a self-self alignment"

split_fa $hap2_asm_name.purged.fa > $hap2_asm_name.purged.split
minimap2 -xasm5 -DP $hap2_asm_name.purged.split $hap2_asm_name.purged.split -t 124 | gzip -c - > $hap2_asm_name.purged.split.self.paf.gz

purge_dups -2 -T cutoffs -c PB.base.cov $hap2_asm_name.purged.split.self.paf.gz > dups_Hap2.bed 2> purge_dups_hap2_purged.log

echo "Get purged primary and haplotig sequences from draft assembly"

get_seqs -e dups_Hap2.bed $hap2_asm_name.purged.fa -p $hap2_asm_name.removed.seq.fa

echo "DONE HAP 2 YAY!"
