# Custom R code to create an file that states which contigs 
# Have a greater than 90 match to the mitogenome AND 
# are shorter than the mitogenome

# read in mitogenome blast results
blast_data_hap1 <- read.table("blastn_mitogenome_hap1.out", header = F)
blast_data_hap2 <- read.table("blastn_mitogenome_hap2.out", header = F)
# Read in mitogenome length
mitolength <- read.table("mitogenome_length.txt")

blast_data_hap1$V16 <- blast_data_hap1[,1]<=mitolength[2,]&blast_data_hap1[,7]>=90
blast_data_hap2$V16 <- blast_data_hap2[,1]<=mitolength[2,]&blast_data_hap2[,7]>=90

write.table(blast_data_hap1, "blastn_mitogenome_hap1_AND.out", row.names = F, col.names = F)
write.table(blast_data_hap2, "blastn_mitogenome_hap2_AND.out", row.names = F, col.names = F)
