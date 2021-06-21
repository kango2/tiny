library(tidyverse)
library(GenomicRanges)
speciesinfo <- read_delim("species.txt", delim="\t")
speciesinfo <- speciesinfo %>% filter(heat_map > 0) %>% arrange(desc(heat_map)) %>% mutate(rphyloorder = row_number())
alnstats <- list()
minalnlen <- 1000 
counter <- 0
for (i in 1:nrow(speciesinfo)){
 for (j in 1:nrow(speciesinfo)){

  if (i == j) {
    next
  }
  qcode <- dplyr::filter(speciesinfo, rphyloorder == i) %>% pull(uniprotspeciescode)
  tcode <- dplyr::filter(speciesinfo, rphyloorder == j) %>% pull(uniprotspeciescode)
  f1 <- paste("lastzaln/",qcode,".vs.",tcode,".target.chainpairs.tab", sep = "")
  f2 <- paste("lastzaln/",tcode,".vs.",qcode,".query.chainpairs.tab", sep = "")
  aln <- NULL
  if (file.exists(f1)){
    aln <- f1
  }
  if (file.exists(f2)){
    aln <- f2
  }
  if (is.null(aln)){
    next
  }
  print(paste(qcode, tcode, aln))
  b <- NULL
  b <- read_delim(aln, "\t",col_names = F)
  colnames(b) <- c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand")
  b <- filter(b, (tend - tstart > minalnlen | qend - qstart > minalnlen))
  ttotal <- GRanges(seqnames = pull(b, tacc), ranges = IRanges(start = pull(b, tstart), end = pull(b, tend)))
  ttotal <- reduce(ttotal)
  qtotal <- GRanges(seqnames = pull(b, qacc), ranges = IRanges(start = pull(b, qstart), end = pull(b, qend)))
  qtotal <- reduce(qtotal)
  counter <- counter + 1
  alnstats[[counter]] <- tibble(qcode = qcode, tcode = tcode, qbases = sum(width(qtotal)), tbases = sum(width(ttotal)))
}
}

x <- bind_rows(alnstats)
x <- left_join(x, select(speciesinfo, tcode = uniprotspeciescode, tgenomelen = genomelength)) 
x <- left_join(x, select(speciesinfo, qcode = uniprotspeciescode, qgenomelen = genomelength))
x <- left_join(x, select(speciesinfo, tcode = uniprotspeciescode, torder = rphyloorder)) 
x <- left_join(x, select(speciesinfo, qcode = uniprotspeciescode, qorder = rphyloorder)) 

pdf("alnproportion.pdf", height = 20, width = 20)
x %>% mutate(qpercent = qbases * 100 / qgenomelen) %>% ggplot(aes(x=torder, y=qorder, fill= qpercent )) + geom_tile() + geom_text(aes(label=sprintf("%0.2f", round(qpercent, digits = 1)))) + scale_fill_viridis_b() + scale_x_continuous(breaks = 1:nrow(speciesinfo), labels = arrange(speciesinfo, rphyloorder) %>%  pull(`Common name`)) + scale_y_continuous(breaks = 1:nrow(speciesinfo), labels = arrange(speciesinfo, rphyloorder) %>%  pull(`Common name`)) + theme_bw() + theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust=1), legend.position="top", legend.title = element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Target species", y = "Query species", fill = "Aln %")
dev.off()
