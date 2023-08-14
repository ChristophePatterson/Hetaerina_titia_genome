#!/bin/bash

#SBATCH -c 124 
#SBATCH --mem=200G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:200G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 72:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

module load bioinformatics
module load blobtoolkit2/3.1.6
module load samtools
module load minimap2
module load blast

unset OMP_PROC_BIND
unset OMP_PLACES

## VARIABLES
BUSCO_DIR=(/nobackup/tmjj24/H_titia_genome/Master/busco)
dir_path_1=(/nobackup/tmjj24/H_titia_genome/Master)
dir_path=($dir_path_1/BlobToolKit)
Prim_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/purge_dups/H_titia.asm.bp.hap1.p_ctg.fa.purged.fa)
Alt_asm_dp=(/nobackup/tmjj24/H_titia_genome/Master/purge_dups/H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa)
Prim_asm=(H_titia.asm.bp.hap1.p_ctg.fa.purged.fa)
Alt_asm=(H_titia.asm.bp.hap2.p_ctg.fa.purged.purged.fa)

Prim_blob=(H_titia_Hifiasm_purged_Hap1)
Alt_blob=(H_titia_Hifiasm_purged_Hap2)
pb_list=(/nobackup/tmjj24/H_titia_genome/HiFiAdapterFilt/multiple_movies.hifi_reads_adaptfilt/multiple_movies.hifi_reads.filt.fastq.gz)

mkdir -p $dir_path

cp -u $Prim_asm_dp $dir_path
cp -u $Alt_asm_dp $dir_path
# cp -u /nobackup/tmjj24/H_titia_genome/H_titia_Hifi_reads/multiple_movies.hifi_reads.fastq.gz $dir_path
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
	--busco $BUSCO_DIR/H_titia_purge_Hifiasm_arthropoda_odb10_Purge_Hap1/run_arthropoda_odb10/full_table.tsv \
	./$Prim_blob

blobtools add \
	--busco $BUSCO_DIR/H_titia_purge_Hifiasm_insecta_odb10_Purge_Hap1/run_insecta_odb10/full_table.tsv \
	./$Prim_blob

# blobtools add \
# 	--busco $BUSCO_DIR/H_titia_purge_Hifiasm_arthropoda_odb10_Purge_Hap2/run_arthropoda_odb10/full_table.tsv \
# 	./$Alt_blob

# blobtools add \
# 	--busco $BUSCO_DIR/H_titia_purge_Hifiasm_insecta_odb10_Purge_Hap2/run_insecta_odb10/full_table.tsv \
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

#Set parameters
Hap1_dir=($Prim_blob)
Hap2_dir=($Alt_blob)
Hap1_fasta=($Prim_asm)
Hap2_fasta=($Alt_asm)
mitogenome=(/nobackup/tmjj24/H_titia_genome/Master/mitohifi/mitohifi_p50_MK/final_mitogenome.fasta)

# calculate length of mitogenome
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' $mitogenome > mitogenome_length.txt
mitogenome_length=`cat mitogenome_length.txt | tail -c 6`
echo $mitogenome_length
cd $dir_path

# Merge all blast files together
cd blastn_array
cat blastn.Hap1*.out > blastn.all.Hap1.out
cat blastn.Hap2*.out > blastn.all.Hap2.out

cd $dir_path

# Add blast hits to blobtools
blobtools add \
    --hits blastn_array/blastn.all.Hap1.out \
    --taxrule bestsum \
    --taxdump /nobackup/tmjj24/H_titia_genome/BlobToolKit/taxdump \
    --replace \
    $Hap1_dir

blobtools add \
    --hits blastn_array/blastn.all.Hap2.out \
    --taxrule bestsum \
    --taxdump /nobackup/tmjj24/H_titia_genome/BlobToolKit/taxdump \
    --replace \
    $Hap2_dir

# Create new fasta files to fill with uncont seq

cat $Hap1_fasta \
 > $Hap1_fasta.uncont.fa

cat $Hap2_fasta \
 > $Hap2_fasta.uncont.fa

# Filter out contaiminates
blobtools filter \
    --param bestsum_phylum--Keys=no-hit,Arthropoda,Chordata \
    --invert \
    --output $Hap1_dir-uncont \
    --fasta $Hap1_fasta.uncont.fa \
    $Hap1_dir

blobtools filter \
    --param bestsum_phylum--Keys=no-hit,Arthropoda,Chordata \
    --invert \
    --output $Hap2_dir-uncont \
    --fasta $Hap2_fasta.uncont.fa \
    $Hap2_dir

# Matching to mitogenome

echo "Starting blast hap1"
blastn -query $mitogenome \
       -subject $Hap1_fasta.uncont.filtered.fa \
       -outfmt "6 slen qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 0.001 \
       -out blastn_mitogenome_hap1.out

echo "Finished blast hap1"

echo "Starting blast hap2"

blastn -query $mitogenome \
       -subject $Hap2_fasta.uncont.filtered.fa \
       -outfmt "6 slen qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 0.001 \
       -out blastn_mitogenome_hap2.out

echo "Finished blast hap2"

## R CODE THAT CALCUALATES WHICH CONTIGS have a greater than 90 percentage match AND are shorter than the mitogenome

R CMD BATCH /home/tmjj24/scripts/job_scripts/Master_genome_assembly/Master_HifiAdapt_Filter/6_R_mitolength_AND_90per_match.R

# Adding mitogenome matches to blobtools

echo "Adding mitogenome matches in blobtools"
blobtools add \
    --text blastn_mitogenome_hap1_AND.out \
    --text-cols 6=identifiers,16=percent_length_mitomatch \
    $Hap1_dir-uncont

blobtools add \
    --text blastn_mitogenome_hap2_AND.out \
    --text-cols 6=identifiers,16=percent_length_mitomatch \
    $Hap2_dir-uncont

echo "Filtering out potential mitogenome contigs"

blobtools filter \
    --param percent_length_mitomatch--Keys=FALSE,NA \
    --invert \
    --output $Hap1_dir-uncont-mito \
    --fasta $Hap1_fasta.uncont.filtered.fa \
    $Hap1_dir-uncont

blobtools filter \
    --param percent_length_mitomatch--Keys=FALSE,NA \
    --invert \
    --output $Hap2_dir-uncont-mito \
    --fasta $Hap2_fasta.uncont.filtered.fa \
    $Hap2_dir-uncont

