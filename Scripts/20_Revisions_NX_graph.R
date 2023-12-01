#libraries
library(ggplot2)
library(ggrepel)

# Revisions plot NX graph from contig lengths
#Read in first table
df <- read.table(paste0("data/Revisions/NX/",list.files("data/Revisions/NX")[1]))
#order by size of contig largest to smallest
df <- df[order(df$V2, decreasing = T),]
df$cum.sum <- cumsum(df$V2)
#Divide by the total length of the assembly
df$NX <- df$cum.sum/sum(df$V2)
# Create new first row
df <- rbind(df[1,], df)
# Add in 0 value to step plot starts at 0
df[1,"NX"] <- 0


plot(df$V2, df$NX)

for(i in 2:length(list.files("data/Revisions/NX"))){
  df.temp <- read.table(paste0("data/Revisions/NX/",list.files("data/Revisions/NX")[i]))
  #order by size of contig largest to smallest
  df.temp <- df.temp[order(df.temp$V2, decreasing = T),]
  #Calculate the cummulative sum
  df.temp$cum.sum <- cumsum(df.temp$V2)
  #Divide by the total length of the assembly
  df.temp$NX <- df.temp$cum.sum/sum(df.temp$V2)
  # Create new first row
  df.temp <- rbind(df.temp[1,], df.temp)
  # Add in 0 value to step plot starts at 0
  df.temp[1,"NX"] <- 0
  df <- rbind(df, df.temp)
}

# Rename columns
colnames(df) <- c("contig", "length", "genome", "cum.sum", "NX")

unique(df$genome)
df$is.titia <- df$genome=="HetTit1.0.p.genome"


## Read data table on species etc
ncbi <- read.table("Data/Revisions/ncbi_dataset.tsv", header = T, sep = "\t")
ncbi[dim(ncbi)[1]+1,] <- c("", "HetTit1.0.p", "Hetaerina titia",rep(NA, dim(ncbi)[2]-3))
ncbi$genome <- paste(ncbi$Assembly.Accession, ncbi$Assembly.Name, "genomic", sep = "_")
ncbi$genome[ncbi$genome=="_HetTit1.0.p_genomic"] <- "HetTit1.0.p.genome"

ncbi$longest_contig <- sapply(ncbi$genome, function(x) max(df$length[df$genome==x]))
ncbi$longest_contig[is.infinite(ncbi$longest_contig)] <- NA


p <- ggplot(df) +
  geom_step(aes(x = NX*100, y = length/1000000, group = genome, col = genome), linewidth = 2, show.legend = F) + 
  geom_text_repel(data = ncbi, aes(x = 0, y = longest_contig/1000000, label = Assembly.Name), 
                  nudge_x = -15, hjust = "outward", nudge_y = 0, size = 10) +
  scale_color_manual(values = c("grey20", "grey25","grey30", "grey35","grey40","grey45","grey50","grey55","grey60","orange")) +
  ylab("Scaffold length (Mbp)") +
  xlab("NX (X = %)") +
  theme_bw() +
  theme(text = element_text(size = 30))

p

ggsave(filename = "Plots/NX_all_odonata.png", p, width = 12, height = 8)

# install.packages("patchwork")
library(patchwork)


pq <- q + p + plot_annotation(tag_levels = "A")

ggsave(filename = "Plots/NX_busco_all_odonata_v2.png", pq, width = 24, height = 12)
ggsave(filename = "Plots/NX_busco_all_odonata_v2.pdf", pq, width = 24, height = 12)
