## Comparing busco scorces
  library(ggplot2)
## Read in busco genes information
insecta_info <- readLines("Data/Revisions/BUSCO/insecta_links_to_ODB10.txt")
insecta_info <- do.call("rbind.data.frame", strsplit(insecta_info, split = "\t"))
colnames(insecta_info) <- c("gene", "info", "web")
arthropoda_info <- readLines("Data/Revisions/BUSCO/arthropoda_links_to_ODB10.txt")
arthropoda_info <- do.call("rbind.data.frame", strsplit(arthropoda_info, split = "\t"))
colnames(arthropoda_info) <- c("gene", "info", "web")

head(insecta_info)
# Which are mitochondrial genes

insecta_info$mito <- grepl(insecta_info$info,pattern = "mitochondrial")
arthropoda_info$mito <- grepl(arthropoda_info$info,pattern = "mitochondrial")
table(insecta_info$mito)


# Get busco run files
busco.runs <- list.files("Data/Revisions/BUSCO/")[grepl(list.files("Data/Revisions/BUSCO/"), pattern = "full_table")]
busco.simple <- list()
busco.missing <- list()
for(i in 1:length(busco.runs)){

busco <- readLines(paste0("Data/Revisions/BUSCO/",busco.runs[i]))
#remove leading lines
busco <- busco[4:length(busco)]

# Create into data frame all the non-missing records
busco.df <- do.call("rbind.data.frame", strsplit(busco[!grepl(busco, pattern = "Missing")], split = "\t"))
colnames(busco.df) <- c("Busco_id", "Status", "Sequence", "Gene_Start", "Gene_End", "Strand", "Score", "Length", "OrthoDB_url", "Description")

# Get all the missing records
busco.m <- do.call("rbind.data.frame", strsplit(busco[grepl(busco, pattern = "Missing")], split = "\t"))
# Add 8 columns of NA on the end so that they are the same size as the completed buscos
busco.m[,3:10] <- NA
colnames(busco.m) <- c("Busco_id", "Status", "Sequence", "Gene_Start", "Gene_End", "Strand", "Score", "Length", "OrthoDB_url", "Description")

# Merging missing data with complete records
busco.df <- rbind(busco.df, busco.m)

# Calculate busco score (need to remove excess duplicate lines)
busco.simple[[i]] <- table(busco.df[!duplicated(paste(busco.df$Busco_id, busco.df$Status)),]$Status)
busco.missing[[i]] <- busco.m$Busco_id
}
busco.simple
busco.df <- as.data.frame(do.call("rbind", busco.simple))
colnames(busco.df) <- c("Single", "Duplicated", "Fragmented", "Missing")

busco.df$n <- rowSums(busco.df)
busco.df$genome <- substr(busco.runs, 1, nchar(busco.runs)-15)
busco.df$Completed <- busco.df$Single+busco.df$Duplicated
busco.df <- busco.df[,c(6,7,1,2,3,4,5)]

busco.df

busco.df$genome.text <- gsub(busco.df$genome, replacement = "", pattern = "^GCA_\\d+\\.\\d+_")

busco.df$genome.text <- gsub(busco.df$genome.text, replacement = "", 
                                pattern = "_genomic_insecta_odb10")
busco.df$genome.text <- gsub(busco.df$genome.text, replacement = "", 
                                pattern = "genomic_arthropoda_odb10")
busco.df$genome.text <- gsub(busco.df$genome.text, replacement = "", 
                                pattern = ".genome_insecta_odb10")



busco.pivot <- tidyr::pivot_longer(busco.df, cols = c("Single", "Duplicated", "Fragmented", "Missing"))
busco.pivot$name <- factor(busco.pivot$name, levels = c("Missing", "Fragmented", "Duplicated", "Single"))


busco.pivot$busco.data <- paste("C:", busco.pivot)
col.bar <-  c("darkred", "yellow","lightblue4","lightblue")

col.bar <-  c("#FF968A", "#FFFFB5","#55CBCD","#ABDEE6")

busco.df.percent <- busco.df
busco.df.percent[,2:6] <- round((busco.df.percent[,2:6]/busco.df.percent$n)*100,digits = 1)
busco.df.percent
busco.df.percent$busco.text <- paste0("C:", busco.df.percent$Completed, " [S:", busco.df.percent$Single, " D:", busco.df.percent$Duplicated, "],",
       " F:", busco.df.percent$Fragmented, ", M:", busco.df.percent$Missing)
busco.pivot$genome.text=="HetTit1.0.p"



# re-order scores based on complete
busco.df.percent$genome <- factor(busco.df.percent$genome, levels = busco.df.percent$genome[order(busco.df.percent$Completed)])
busco.pivot$genome <- factor(busco.pivot$genome, levels = busco.df.percent$genome[order(busco.df.percent$Completed)])

# Order of scores
busco.df.percent.insect <- busco.df.percent[grepl(busco.df.percent$genome, pattern = "insecta"),]
busco.pivot.insect <- busco.pivot[grepl(busco.pivot$genome, pattern = "insecta"),]



q <- ggplot(busco.pivot.insect) +
  geom_bar(aes(y = genome.text, x = (value/n)*100, fill = name), stat = "identity") +
  geom_text(data = busco.df.percent.insect, aes(y = genome.text, x = 5, label = busco.text),size = 5, hjust = "inward") +
  geom_point(aes(y = genome.text, x = 2, col = genome.text=="HetTit1.0.p"), size = 7, shape = 18, show.legend = F) +
  scale_fill_manual(values = col.bar) +
  scale_color_manual(values = c(rgb(0,0,0,0), "orange")) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 30)) +
  xlab("Percentage of busco genes (%)") +
  ylab("Draft genome") +
  labs(fill = element_blank()) +
  scale_y_discrete(limits = busco.df.percent.insect$genome.text[order(as.numeric(busco.df.percent.insect$Single))])


q  

ggsave(filename = "Plots/BUSCO_all_odonata_2.png", q, width = 12, height = 8)

# Which genes are missing across multiple busco runs
missing.genes <- cbind(busco.missing[[1]])
for(i in 2:length(busco.runs)){
  temp <- cbind(busco.missing[[i]], busco.runs[i])
  missing.genes <- rbind(missing.genes, temp)
}

#Which genes are not found in any genome assemblies
missing.genes.all <- names(table(missing.genes[,1])[table(missing.genes[,1])==9])
insecta_info[insecta_info$gene%in%missing.genes.all,]

write.table(busco.df, file = "Data/Revisions/BUSCO/insecta_all_genome.txt",quote = F, row.names = F, sep = "\t")
write.table(busco.df.percent, file = "Data/Revisions/BUSCO/insecta_all_genome_percentage.txt",quote = F, row.names = F, sep = "\t")
