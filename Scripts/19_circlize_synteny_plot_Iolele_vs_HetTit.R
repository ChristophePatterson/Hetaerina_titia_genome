# PLotting synteny plots using 
# https://www.royfrancis.com/beautiful-circos-plots-in-r/
# https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#concatenating-two-genomes

# install.packages("circlize")

library(circlize)
# Needed data
#input.file <- "HetTit1.0_vs_HetAmer1.0"
input.file <- "IolEle1.2_vs_HetTit1.0"
spp.1 <- strsplit(input.file,"_")[[1]][1]
spp.1 <- substr(spp.1, 1, nchar(spp.1)-3)
spp.2 <- strsplit(input.file,"_")[[1]][3]
spp.2 <- substr(spp.2, 1, nchar(spp.2)-3)
version.input <- "_v7"
#version.input <- ".v8"

# Simple plot
chrom_length <- read.table(paste0("Data/Jupiter/",input.file,".Chromosome_length.txt"))
names(chrom_length) <- c("chr", "length", "Spp")
chrom_length <- chrom_length[!grepl(chrom_length$chr, pattern = "CAKLCU02")|chrom_length$Spp=="HetTit",]

# Percentage length of chrom
order(chrom_length[chrom_length$Spp==spp.2,]$length)==75:1
(cumsum(chrom_length[chrom_length$Spp==spp.2,]$length)/sum(chrom_length[chrom_length$Spp==spp.2,]$length))[12]


#Extract chrom size spp 1
spp.1.chrom <- chrom_length[chrom_length$Spp==spp.1,]
spp.1.chrom$cum.length <- cumsum(spp.1.chrom$length)
spp.1.chrom$start <- spp.1.chrom$cum.length-spp.1.chrom$length
spp.1.chrom$end <- spp.1.chrom$cum.length
spp.1.chrom$chr.old <- spp.1.chrom$chr
spp.1.chrom <- spp.1.chrom[!grepl(spp.1.chrom$chr, pattern = "CAKLCU02"),]
duplicated(spp.1.chrom$chr)
spp.1.chrom$chr <- paste0(spp.1,"-chr",as.numeric(substr(sapply(strsplit(spp.1.chrom$chr, split = "\\."), "[", 1 ),nchar(spp.1.chrom$chr)-3,nchar(spp.1.chrom$chr))))
#remorder largest to smallest
spp.1.chrom <- spp.1.chrom[order(spp.1.chrom$end-spp.1.chrom$start),]


#Extract chrom spp 1
spp.2.chrom <- chrom_length[chrom_length$Spp==spp.2,]
spp.2.chrom$cum.length <- cumsum(spp.2.chrom$length)
spp.2.chrom$start <- spp.2.chrom$cum.length-spp.2.chrom$length
spp.2.chrom$end <- spp.2.chrom$cum.length
spp.2.chrom$chr.old <- spp.2.chrom$chr
spp.2.chrom$chr <- paste0(spp.2,"-chr",sapply(strsplit(spp.2.chrom$chr, split = "-"), "[", 2 ))
spp.2.chrom <- spp.2.chrom[order(spp.2.chrom$end-spp.2.chrom$start),]
spp.2.chrom <- spp.2.chrom[spp.2.chrom$end-spp.2.chrom$start>10000000,]
#remorder largest to smallest
spp.2.chrom <- spp.2.chrom[order(spp.2.chrom$end-spp.2.chrom$start),]

#combine
chrom <- rbind(spp.1.chrom[,c(1,5,6,3,7)], spp.2.chrom[,c(1,5,6,3,7)])

# Input names of chromosome into links dataframe
Karyotype <- read.table(paste0("Data/Jupiter/",input.file,version.input,".Karyotype"))
Karyotype <- Karyotype[Karyotype$V1=="chr",]

links.final <- read.table(paste0("Data/Jupiter/",input.file,version.input,".links.final"), sep = " ")[,1:6]
names(links.final) <- c("Reference.chromosome.ID","Reference.start.position", "Reference.end.position","Target.chromosome","Target.start.position", "Target.end.position") 

links.final$Reference.chromosome.ID <- Karyotype$V4[do.call("c",lapply(links.final$Reference.chromosome.ID, function(x) which(Karyotype$V3==x)))]
links.final$Target.chromosome <- Karyotype$V4[do.call("c",lapply(links.final$Target.chromosome, function(x) which(Karyotype$V3==x)))]

links.final$Reference.chromosome.ID <- paste0(spp.1,"-chr",as.numeric(substr(sapply(strsplit(links.final$Reference.chromosome.ID , split = "\\."), "[", 1 ),nchar(links.final$Reference.chromosome.ID)-3,nchar(links.final$Reference.chromosome.ID))))
links.final$Target.chromosome <- paste0(spp.2,"-chr",sapply(strsplit(links.final$Target.chromosome, split = "-"), "[", 2 ))

links.final$Reference.start.position <- links.final$Reference.start.position+spp.1.chrom$start[do.call("c",lapply(links.final$Reference.chromosome.ID, function(x) which(spp.1.chrom$chr==x)))]
links.final$Reference.end.position <- links.final$Reference.end.position+spp.1.chrom$start[do.call("c",lapply(links.final$Reference.chromosome.ID, function(x) which(spp.1.chrom$chr==x)))]
links.final$Target.start.position <- links.final$Target.start.position+spp.2.chrom$start[do.call("c",lapply(links.final$Target.chromosome, function(x) which(spp.2.chrom$chr==x)))]
links.final$Target.end.position <- links.final$Target.end.position+spp.2.chrom$start[do.call("c",lapply(links.final$Target.chromosome, function(x) which(spp.2.chrom$chr==x)))]

link1 <- links.final[, c("Reference.chromosome.ID","Reference.start.position", "Reference.end.position")]
link2 <- links.final[, c("Target.chromosome","Target.start.position", "Target.end.position")]
names(link1) <- c("chr", "start","end")
names(link2) <- c("chr", "start","end")


# Get order of spp2 contigs in relation to there position on spp1
library(tidyverse)
ref.target.order <- group_by(links.final, by = Target.chromosome) %>%
  summarise(mn = mean(Reference.start.position))
ref.target.order <- ref.target.order$by[order(ref.target.order$mn)]

ref.target.order <- c(spp.1.chrom$chr, ref.target.order)
chrom$chr %in% ref.target.order

chrom <- chrom[do.call("c",lapply(ref.target.order, function(x) which(chrom$chr==x))),]

chrom$chr[chrom$chr=="IolEle-chr6"] <- "IolEle-chrX"
chrom$chr[chrom$chr=="HetTit-chr12"] <- "HetTit-chrX"
link1$chr[link1$chr=="IolEle-chr6"] <- "IolEle-chrX"
link2$chr[link2$chr=="HetTit-chr12"] <- "HetTit-chrX"

chrom$chr[chrom$chr=="IolEle-chr13"] <- "IolEle-chrM"
link1$chr[link1$chr=="IolEle-chr14"] <- "IolEle-chrM"



# shift X chromosome to the end
X_order <- c("IolEle-chrX",chrom$chr[chrom$Spp==spp.1][!grepl(chrom$chr[chrom$Spp==spp.1], pattern = "X")],
             chrom$chr[chrom$Spp==spp.2][!grepl(chrom$chr[chrom$Spp==spp.2], pattern = "X")],"HetTit-chrX")
chrom <- chrom[do.call("c",lapply(X_order, function(x) which(chrom$chr==x))),]

col.data.spp <- cbind(c(spp.1, spp.2), c("lightblue","orange"))
chrom$spp.col <- do.call("c", lapply(chrom$Spp, function(x) col.data.spp[,2][which(col.data.spp[,1]==x)]))


png(filename = paste0("Plots/Synteny_plot_", input.file,".png"),width = 5000, height = 5000)

circos.initializeWithIdeogram(chrom, plotType = NULL, 
                              chromosome.index = chrom$chr)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  
  circos.rect(xlim[1], 0, xlim[2],1,lwd = 2, col = chrom$spp.col[which(chrom$chr==CELL_META$sector.index)])
  
}, track.height = 0.15, bg.border = NA)

#highlight.chromosome(chrom$chr[chrom$Spp==spp.2], 
#                     col = "darkolivegreen3", track.index = 1, padding = c(-1,0,-2,0))

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chrom.text <- gsub(".*chr", "", CELL_META$sector.index)
  if(chrom.text!="X"&chrom.text!="M"){
    if(as.numeric(chrom.text)>12){
      chrom.text <- ""
    }}
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(50), 
              chrom.text, cex = 8, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

highlight.chromosome(chrom$chr[chrom$chr=="HetTit-chr1"], 
                     col = "orange", track.index = 2, padding = c(5,0,8,0))
highlight.chromosome(c("IolEle-chr12", "IolEle-chr9"), 
                     col = "lightblue", track.index = 2, padding = c(5,0,8,0))

set.seed(1)
colour.pall <- rand_color(length(unique(link1$chr)),luminosity = "bright", transparency = 0)[sapply(link1$chr, function(x) which(unique(link1$chr)==x))]
# plot(1:length(colour.pall), col = colour.pall, cex = 15, pch = 19)

length(link1$chr)
circos.genomicLink(link1, link2, col = colour.pall, lwd = 5)
text(-0.9, -0.8, substitute(italic("I. elegans")), cex = 10)
text(0.9, 0.8, substitute(italic("H. titia")), cex = 10)
text(-0.9, 0.9, "A", cex = 15)

dev.off()

#################################################

pdf(file = paste0("Plots/Synteny_plot_", input.file,".pdf"))

circos.initializeWithIdeogram(chrom, plotType = NULL, 
                              chromosome.index = chrom$chr)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  
  circos.rect(xlim[1], 0, xlim[2],1, col = chrom$spp.col[which(chrom$chr==CELL_META$sector.index)])
  
}, track.height = 0.15, bg.border = NA)

#highlight.chromosome(chrom$chr[chrom$Spp==spp.1], 
#                     col = "orange", track.index = 1)
#highlight.chromosome(chrom$chr[chrom$Spp==spp.2], 
#                     col = "darkolivegreen3", track.index = 1, padding = c(-1,0,-2,0))

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chrom.text <- gsub(".*chr", "", CELL_META$sector.index)
  if(chrom.text!="X"&chrom.text!="M"){
    if(as.numeric(chrom.text)>12){
      chrom.text <- ""
    }}
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(5), 
              chrom.text, cex = 0.6, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

highlight.chromosome(chrom$chr[chrom$chr=="HetTit-chr1"], 
                     col = "orange", track.index = 2, padding = c(1,0,1,0))
highlight.chromosome(c("IolEle-chr12", "IolEle-chr9"), 
                     col = "lightblue", track.index = 2, padding = c(1,0,1,0))

set.seed(1)
colour.pall <- rand_color(length(unique(link1$chr)),luminosity = "bright", transparency = 0)[sapply(link1$chr, function(x) which(unique(link1$chr)==x))]
# plot(1:length(colour.pall), col = colour.pall, cex = 15, pch = 19)

length(link1$chr)
circos.genomicLink(link1, link2, col = colour.pall)
text(-0.9, -0.8, substitute(italic("I. elegans")), cex = 1)
text(0.9, 0.8, substitute(italic("H. titia")), cex = 1)
text(-0.9, 0.9, "A", cex = 2)

dev.off()

