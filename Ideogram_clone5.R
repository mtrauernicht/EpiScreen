#ideogram of insertions on chromosomes
#what you need:
#a file with chromosome lengths (chromsizes) for plotting chromosomes as bar graph (example file attached, downloaded from UCSC)
#a data frame with your integrations  to add lines to the bar graph, format c(â€œchrâ€,â€posâ€)
#I split the integrations into several data frames by expression groups (ins_lo etc.) because I wanted to plot them in discrete colors
#below I also added another version of the graph with insertions color-coded by expression as gradient
library(ggplot2)
chromsizes<-read.table("/home/r.schep/mydata/data/genomes/GRCh38/hg38_chromsizes2.txt",header=F, stringsAsFactors = F)
chromsizes
names(chromsizes)<-c("chr_r","pos_r")

granges.clone5 <- read.csv2("/DATA/usr/m.trauernicht/projects/EpiScreen/Mapping_Clone5_9/mt20190509_clone5_genomicranges.csv")
setnames(granges.clone5, "start", "start_pos")
setnames(granges.clone5, "seqnames", "seqname")

mutationintegrations <- granges.clone5[,c("seqname","start_pos")]
mutationintegrations$start_pos <- mutationintegrations$start_pos+500
colnames(mutationintegrations) <- c("chr_r","pos_r")
dim(mutationintegrations)

# Added by Tom
changeLevels <- function(object, column = "chr_r", levels = c(paste0("chr", 1:22), "chrX")) {
  object[, column] <- as.character(object[, column])
  object[, column] <- factor(object[, column], levels)
  object
}

mutationintegrations <- changeLevels(mutationintegrations)
chromsizes <- changeLevels(chromsizes)

dim(mutationintegrations)

chromsizes <- chromsizes[!is.na(chromsizes$chr_r), ]
mutationintegrations <- mutationintegrations[!is.na(mutationintegrations$chr_r), ]


p <- ggplot(chromsizes, aes(x=chr_r,y=as.numeric(pos_r)))

# All integrations with mutation data
p +  geom_bar(stat='identity',fill="#e2e2e2", color="black")+
  geom_hline(data = mutationintegrations, 
             aes(yintercept= as.numeric(pos_r)),
             color="navy") +
  facet_grid(. ~ chr_r,
             scale = "free",
             space = "free_x")+
  scale_y_continuous(breaks = c(50000000, 
                                100000000, 
                                150000000, 
                                200000000,
                                250000000), 
                     labels=c("50 Mb",
                              "100 Mb",
                              "150 Mb",
                              "200 Mb",
                              "250 Mb"))+
  theme_bw()
theme(axis.title.y = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle=90, 
                                 size = 15),
      axis.title.x = element_blank(), 
      axis.text.x = element_blank(),
      strip.text.x = element_text(angle= 90))

