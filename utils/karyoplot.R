library(tidyverse)
library(plotly)
library(ggplotlyExtra)

args = commandArgs(trailingOnly=TRUE)

setwd("/g/data/te53/hrp561/wga/tiny")

speciesinfo <- read_delim("./metadata/species.txt", delim = "\t")
chrinfo <- read_delim("./metadata/chrmergedassemblyreport.txt", delim="\t")

setwd("/g/data/te53/hrp561/wga/lastzaln")

files <- c(list.files(".", "*.target.chainpairs.tab"), list.files(".", "*.query.chainpairs.tab"))

i <- as.numeric(args[1])
print(files[i])
chainfile <- files[i]
qspecies <- str_split(chainfile, "\\.", simplify = T)[1]
tspecies <- str_split(chainfile, "\\.", simplify = T)[3]
qcname <- filter(speciesinfo, uniprotspeciescode == qspecies) %>% pull("Common name")
tcname <- filter(speciesinfo, uniprotspeciescode == tspecies) %>% pull("Common name")
labelx <- paste(tcname, " (", tspecies, ")", sep = "")
labely <- paste(qcname, " (", qspecies, ")", sep = "")

if (str_split(chainfile, "\\.", simplify = T)[4] == "query"){
  labelx <- paste(qcname, " (", qspecies, ")", sep = "")
  labely <- paste(tcname, " (", tspecies, ")", sep = "")
}
  
b <- read_delim(chainfile, "\t",col_names = F)
colnames(b) <- c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand")
b <- left_join(b, dplyr::select(chrinfo, tacc = seqid, tchr = chrname, trank = krank, toffset=offset), by = "tacc") 
b <- left_join(b, dplyr::select(chrinfo, qacc = seqid, qchr = chrname, qrank = krank, qoffset=offset), by = "qacc")
b <- mutate(b, tlen=tend-tstart, qlen=qend-qstart)
  
  
##with colors for strand and displaying strands correctly
minalnlen <- 150
g <- drop_na(b) %>% dplyr::filter(qlen>minalnlen | tlen>minalnlen) %>% 
  mutate(adjstart = if_else(strand=="+", qstart, qend), adjend=if_else(strand=="+",qend,qstart)) %>% 
  ggplot() + 
  geom_segment(aes(x=tstart+toffset, y=adjstart+qoffset, xend=tend+toffset, yend=adjend+qoffset,color=strand)) + 
  theme_bw() + 
  scale_x_continuous(breaks = drop_na(b) %>% arrange(trank) %>% pull(toffset) %>% unique(), minor_breaks = NULL, labels = drop_na(b) %>% arrange(trank) %>% pull(tchr) %>% unique()) + 
  scale_y_continuous(breaks = drop_na(b) %>% arrange(qrank) %>% pull(qoffset) %>% unique(), minor_breaks = NULL, labels = drop_na(b) %>% arrange(qrank) %>% pull(qchr) %>% unique()) +
  xlab(labelx) +
  ylab(labely)
pdf(paste("../dotplots/", paste(str_split(chainfile, "\\.", simplify = T)[1:5], collapse = "."), ".pdf", sep = ""), width = 15, height = 15)
g
dev.off()
#dotplot <- ggplotly(g)
#htmlwidgets::saveWidget(dotplot, paste("../dotplots/", paste(str_split(chainfile, "\\.", simplify = T)[1:5], collapse = "."), ".html", sep = ""), selfcontained=F,  libdir="/g/data/te53/hrp561/wga/dotplots/lib")

##calculate pairwise alignment stats
##prepare alluvial plots
