library(tidyverse)
library(Biostrings)

speciesinfo <- read_delim("./metadata/species.txt", delim = "\t")

tlist <- list()
for (i in 1:nrow(speciesinfo)) {
  print(i)
  print(speciesinfo[i,"uniprotspeciescode"])
  asminfo <- read_delim(paste("./metadata/", str_match(speciesinfo[i,"assemblyreport"], "\\S.*\\/(\\S+$)")[,2], sep = ""), delim="\t", comment = "#", col_names = F)
  colnames(asminfo) <- c("seqname", "seqrole", "assignedmol", "assignedlocationtype", "genbankacc", "relation", "refseqacc", "asmunit", "seqlength", "ucscname")
  asminfo <-
    mutate(asminfo, speciescode = as.character(speciesinfo[i,2]),
           cname = as.character(speciesinfo[i,3]),
           krank = row_number(),
           offset = lag(seqlength),
           offset = replace_na(offset, 0),
           offset = cumsum(offset)
    )
  fasta <- readDNAStringSet(paste("../genomes/", speciesinfo[i,"uniprotspeciescode"], "/", str_match(speciesinfo[i,"fasta"], "\\S.*\\/(\\S+$)")[,2], sep = ""))
  l <- data.frame(letterFrequency(fasta,c("A","C","G","T")))
  l <- bind_cols(l, seqid = str_extract(names(fasta), "(^\\S+)"))
  l <- mutate(l, gc = ((G + C) * 100) / (A + C + G + T)) 
  fasta <- NULL
  r <- sum(asminfo$refseqacc %in% l$seqid)
  g <- sum(asminfo$genbankacc %in% l$seqid)
  if (g > r) {
  	asminfo$seqid <- asminfo$genbankacc
  } else {
  	asminfo$seqid <- asminfo$refseqacc
  }
  asminfo <- left_join(asminfo, l, by = c('seqid'))
  asminfo <- mutate(asminfo, chrname = 
  		case_when(
  			seqrole == "assembled-molecule" & assignedlocationtype == "Chromosome" ~ paste("chr",assignedmol, sep = ""),
  			seqrole == "assembled-molecule" & assignedlocationtype == "Linkage Group" ~ assignedmol,
  			seqrole == "assembled-molecule" & assignedlocationtype == "Mitochondrion" ~ "ChrM",
  			seqrole == "unplaced-scaffold" ~ seqname,
  			seqrole == "unlocalized-scaffold" ~ seqname
  			)
  		)
  tlist[[i]] <- asminfo
}

asminfo <- bind_rows(tlist)
write_delim(asminfo, "./metadata/mergedassemblyreport.txt", delim="\t")
write_delim(filter(asminfo, seqrole == "assembled-molecule" & assignedlocationtype != "Mitochondrion"), "./metadata/chrmergedassemblyreport.txt", delim="\t")

