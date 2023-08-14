# Hetaerina_titia_genome
A repository for code used to construct the draft genome of Hetaerina titia

Draft genome of a female *Hetaerina titia* (Sample ID - HXRCb13) collected from Santa María Huatulco, Oaxaca, Mexico (15.74, -96.298) in June 2021 and stored in RNAlater. 

![P1010242_HiFi_Female_HXRCb13](https://github.com/ChristophePatterson/Hetaerina_titia_genome/docs/P1010242_HiFi_Female_HXRCb13.jpeg)

Figure 1: HXRCb13 collected from María Huatulco, Oaxaca, Mexico (15.74, -96.298) in June 2021.

All data, scripts, and output files used to create the final draft genome, including the contigs only version and the final scaffolded genome are stored within the folder `Scripts`. 

# Raw PacBio Hifi reads
The raw PacBio reads, sequenced using one 8M SMRT Cell, Sequel II sequencing chemistry 2.0, and 10-hr movies on a PacBio Sequel II sequencer, are achieved on NCBI Sequence read achieve (SRR23023424 - under embargo until publication). Raw reads were filtered for pacbio adapters using `Hifiadaptfilt` removing 581 reads.

Estimation of genome from is [archieved here](http://genomescope.org/analysis.php?code=f6FBTv63MLxk5BCyfcep)

# Draft genome

The final primary and alternative draft genome are upload to NCBI (underembargo until publication)

The following programmes were used during the contruction.

- HifiAdapterFilt - for filtering out of remanant pacbio adapters from the raw sequence data.
- Jellyfish/Genomescope - for estimation of genome size and statistics built from raw reads.
- Hifiasm - for initial contig assembly
- Purge_dups - for removal of haplotig duplications in the primary and alternative assembly
- MitoHifi - mitogenome assembly
- Mitos - mitogenome annotation
- Blob toolkit/Blastn - contamination and figure creation
- HiRise for Hi-C scaffold construction
- BUSCO - genome completion statistics
- QUAST - Genome statistics

# Scripts

All scripts used to construct the draft genome are stored within the directory `Scripts`. Scripts are ordered from 0 to 19 and indicate which programme or stage the each should be used. Each script is separated into batch jobs that can be ran within the 3 day time limit of Hamilton (Durham Univerisity's HPC). Scripts 2 and 3 are contained within script 1 and do not need to be run. All, but one script, are written as slurm scripts and can be submitted to the Hamilton HPC. File `5_srun_mitohifi.sh` requires manually opening an interactive slurm job and then an interactive singularly image -  details of which are contained within the script. Scripts 4_2 can be run in parallel. 

Because the final scaffolded assembly was created using the proprietary software HiRise it is not possible, currently, to construct the Hi-C scaffolded assembly ourselves

From start to finish running all scripts will take approximately 5 days.

- 0_jellyfish.sh
- 1_HifiAdatp-Hifiasm-purge_dups.sh
	- 2_Hifiasm-purge_dups.sh
	- 3_Purge_dups.sh
- 4_File_split_busco.sh
- 4_2_Blastn_uncont_Hap1_array.sh
- 4_2_Blastn_uncont_Hap2_array.sh
- 5_srun_mitohifi.sh
- 6_blobtools_contaimination.sh
	- 6_R_mitolength_AND_90per_match.R
- 7_RagTag.sh
- 11_1_Blob_tools_Busco_Scaff.sh
- 11_2_Blob_tools_Scaff.sh
- 12_quast_all_asm.sh
- 12_quast.sh
- 13_buco_HetTit1.0.sh
- 14_Blobtools_HetTit1.0.sh
- 15_Jupiter_plot.sh
- 16_lastz_HetTit_HetAme.sh
- 16_lastz_IsoEle_vs_HetTit.sh
- 18_jupiter_plot_Cal_vs_Het
- 19_jupiter_plot_HetTit_vs_Het_Amer
	- 18_circlize_synteny_plot_HetTit_vs_HetAmer.R
	- 18_circlize_synteny_plot_Iolele_vs_HetTit.R


# Hi-C scaffolding

Hi-C sequencing was attained from a Hetaerina titia collected from NMSC01 in April 2023. Sequencing was done by [Cantata-bio](https://cantatabio.com/). A new scaffolded assembly was received for the primary assembly. This new assembly is highly contiguous and contains 12 chromosome scale scaffolds (75 scaffolds in total).

# Genome assembly overview - Quast

Overviews of genome assembly can be viewed using the output of quest stored within the directory `Quast`.

Two Quast reports are retained. Both quast reports can be view in either `.html` or `.pdf` format. `html` format is recommended.

### (1) Detailed Quast output for the scaffolded primary and alternative assemblies.

Detailed summary statistics for the primary and alternative scaffolded assemblies are retained within the subdirectory ``HetTit1.0`.`.

### (2) Basic Quast output for all assemblies
Basic summary statistics for all assemblies covering the entire pipeline of genome assembly can be found within the subdirectory `H_titia_All_asm_HetTit1.0`. 

# Genome QC table

|                                               |   |               |       | Primary    |                       | Alternate  |
|-----------------------------------------------|---|---------------|-------|------------|-----------------------|------------|
| No. contigs                                   |   |               |       | 4094       |                       | 4054       |
| Contig N50 (bp)                               |   |               |       | 639        |                       | 664        |
| Contig L50                                    |   |               |       | 638359     |                       | 631942     |
| Longest contigs                               |   |               |       | 4037799    |                       | 3160621    |
| Number of scaffolds                           |   |               |       | 75         |                       | NA         |
| Scaffold N50                                  |   |               |       | 120343728  |                       | NA         |
| Scaffold L50                                  |   |               |       | 6          |                       | NA         |
| Largest scaffold                              |   |               |       | 151479759  |                       | NA         |
| Size of final assembly (bp)                   |   |               |       | 1444697970 |                       | 1429703052 |
| BUSCO completeness (insect_odb10) n=1367      |   | C             | S     | D          | F                     | M          |
|                                               | P | 97.0%         | 96.0% | 1.0%       | 1.6%                  | 1.4%       |
|                                               | A | 95.9%         | 94.9% | 1.0%       | 2.0%                  | 2.1%       |
| BUSCO completeness (arthropoda_odb10) n= 1013 |   | C             | S     | D          | F                     | M          |
|                                               | P | 97.6%         | 97.3% | 0.3%       | 0.8%                  | 1.6%       |
|                                               | A | 96.2%         | 95.7% | 0.5%       | 1.3%                  | 2.5%       |
| Organelle                                     |   | Mitochondrion |       |            | NCBI Accession Number | OQ363879   |


---




