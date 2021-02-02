library(tidyverse)
library(plotly)
library(ggplotlyExtra)

speciesinfo <- read_delim("metadata/species.txt", delim = "\t")
speciesinfo <- speciesinfo %>% mutate(lareport = str_match(assemblyreport, "\\S.*\\/(\\S+$)")[,2])

##process assembly report files
tlist <- list()
for (i in 1:nrow(speciesinfo)) {
  asminfo <- read_delim(paste("metadata/", speciesinfo[i,"lareport"], sep = ""), delim="\t", comment = "#", col_names = F)
  colnames(asminfo) <- c("seqname",	"seqrole", "assignedmol", "assignedlocationtype", "genbankacc", "relation", "refseqacc", "asmunit", "seqlength", "ucscname")
  asminfo <- 
    filter(asminfo, seqrole == "assembled-molecule" & assignedlocationtype != "Mitochondrion") %>% 
    mutate(speciescode = as.character(speciesinfo[i,2]), 
           cname = as.character(speciesinfo[i,3]),
           krank=row_number(),
           offset = lag(seqlength),
           offset = replace_na(offset, 0),
           offset = cumsum(offset)
    )
  tlist[[i]] <- asminfo
}
asminfo <- bind_rows(tlist)
asminfo <- mutate(asminfo, chrname = case_when(assignedlocationtype == "Chromosome" ~ paste("chr",assignedmol, sep = ""), assignedlocationtype == "Linkage Group" ~ assignedmol))

#> files
#[1] "GCF_000002315.6_GRCg6a.vs.GCA_000241765.4_Chrysemys_picta_BioNano-3.0.4.chainpairs.tab" "GCF_000002315.6_GRCg6a.vs.GCA_002798355.1_Ogye1.0.chainpairs.tab"          
#[3] "GCF_000002315.6_GRCg6a.vs.GCA_003400415.2_UTA_CroVir_3.0.chainpairs.tab"                "GCF_000002315.6_GRCg6a.vs.GCA_009733165.1_Nana_v5.chainpairs.tab"          
#[5] "GCF_000002315.6_GRCg6a.vs.GCA_009764565.2_rDerCor1.pri.v2.chainpairs.tab"               "GCF_000002315.6_GRCg6a.vs.GCA_009769625.1_bCygOlo1.pri.chainpairs.tab"     
#[7] "GCF_000002315.6_GRCg6a.vs.GCA_013407035.1_ASM1340703v1.chainpairs.tab"                  "GCF_000002315.6_GRCg6a.vs.GCA_014706415.1_LU_Pmuni_1.1.chainpairs.tab"     
#[9] "GCF_000002315.6_GRCg6a.vs.GCA_015237465.1_rCheMyd1.pri.chainpairs.tab"                  "GCF_000002315.6_GRCg6a.vs.GCF_000002315.6_GRCg6a.chainpairs.tab"           
#[11] "GCF_000002315.6_GRCg6a.vs.GCF_000090745.1_AnoCar2.0.chainpairs.tab"                     "GCF_000002315.6_GRCg6a.vs.GCF_004115215.1_mOrnAna1.p.v1.chainpairs.tab"   
#[13] "GCF_000002315.6_GRCg6a.vs.GCF_004329235.1_PodMur_1.0.chainpairs.tab"                    "GCF_000002315.6_GRCg6a.vs.GCF_007399415.2_rGopEvg1_v1.p.chainpairs.tab"   
#[15] "GCF_000002315.6_GRCg6a.vs.GCF_009769535.1_rThaEle1.pri.chainpairs.tab"                  "GCF_000002315.6_GRCg6a.vs.GCF_009819535.1_rLacAgi1.pri.chainpairs.tab"        
#[17] "GCF_000002315.6_GRCg6a.vs.GCF_011800845.1_UG_Zviv_1.chainpairs.tab"                     "GCF_000002315.6_GRCg6a.vs.GCF_013100865.1_CAS_Tse_1.0.chainpairs.tab"    

setwd("data/")
files <- c(list.files(".", "*.chainpairs.target.tab"), list.files(".", "*.chainpairs.query.tab"))
seqinfo <- bind_rows(select(asminfo, accession = genbankacc, seqlength, speciescode, krank, offset, chrname), 
                     select(asminfo, accession = refseqacc, seqlength, speciescode, krank, offset, chrname))
for (i in 1:length(files)) {
  print(files[i])
  chainfile <- files[i]
  qspecies <- str_split(chainfile, "\\.", simplify = T)[1]
  tspecies <- str_split(chainfile, "\\.", simplify = T)[3]
  qcname <- filter(speciesinfo, uniprotspeciescode == qspecies) %>% pull("Common name")
  tcname <- filter(speciesinfo, uniprotspeciescode == tspecies) %>% pull("Common name")
  labelx <- paste(tcname, " (", tspecies, ")", sep = "")
  labely <- paste(qcname, " (", qspecies, ")", sep = "")
  
  if (str_split(chainfile, "\\.", simplify = T)[5] == "query"){
    labelx <- paste(qcname, " (", qspecies, ")", sep = "")
    labely <- paste(tcname, " (", tspecies, ")", sep = "")
  }
  
  b <- read_delim(chainfile, "\t",col_names = F)
  colnames(b) <- c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand")

  b <- left_join(b, dplyr::select(seqinfo, tacc = accession, tchr = chrname, trank = krank, toffset=offset), by = "tacc") 
  b <- left_join(b, dplyr::select(seqinfo, qacc = accession, qchr = chrname, qrank = krank, qoffset=offset), by = "qacc")
  b <- mutate(b, tlen=tend-tstart, qlen=qend-qstart)
  
  ##calculate pairwise alignment stats
  
  #g <- b %>% mutate(tmid=toffset + ((tstart+tend)/2), qmid=qoffset + ((qstart+qend)/2)) %>% dplyr::filter(qlen>200) %>% ggplot() + geom_segment(aes(x=tstart+toffset, y=qstart+qoffset, xend=tend+toffset, yend=qend+qoffset,)) + theme_bw() + scale_x_continuous(breaks = arrange(b, trank) %>% pull(toffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, trank) %>% pull(tchr) %>% unique()) + scale_y_continuous(breaks = arrange(b, qrank) %>% pull(qoffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, qrank) %>% pull(qchr) %>% unique())
  #dotplot <- ggplotly(g)
  #dotplot
  #htmlwidgets::saveWidget(dotplot, "dotplot.html")
  
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
  dotplot <- ggplotly(g)
  dotplot
  htmlwidgets::saveWidget(dotplot, paste(paste(str_split(chainfile, "\\.", simplify = T)[1:5], collapse = "."), ".html", sep = ""))
}
