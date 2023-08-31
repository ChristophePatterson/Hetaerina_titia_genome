# Reading lastz output

blastn.names <- c("query id",  "subject_id",  "Per_identity",  "alignment_length",  "mismatches",  "gap_opens",  "q.start",  "q.end",  "s.start",  "s.end",  "evalue",  "bit_score")
Ele_tit <- read.table("data/Lastz/IolEle_vs_HetTit.dat")
Ele_amer <- read.table("data/Lastz/IolEle_vs_HetAmer.dat")
Tit_amer <- read.table("data/Lastz/HetAme_vs_HetTit.dat")
colnames(Ele_tit) <- blastn.names
colnames(Ele_amer) <- blastn.names
colnames(Tit_amer) <- blastn.names

# Sum of difference / sum of alligned length
# Elegans vs titia
(1-sum(Ele_tit$mismatches)/sum(Ele_tit$alignment_length))*100
# Elegans vs amer
(1-sum(Ele_amer$mismatches)/sum(Ele_amer$alignment_length))*100
# Amer vs titia
(1-sum(Tit_amer$mismatches)/sum(Tit_amer$alignment_length))*100
