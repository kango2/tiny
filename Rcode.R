library(tidyverse)
library(ggplotlyExtra)

setwd("../micro/")

a <- read_delim("chrnames.txt", "\t")
a <- a %>% group_by(assemblyid) %>% mutate(krank=row_number(), offset = lag(length), offset = replace_na(offset, 0), offset = cumsum(offset)) %>% ungroup()

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

files <- list.files(".", "GCF_000002315.6_GRCg6a.vs.*")

b <- read_delim(files[4], "\t",col_names = F)
colnames(b) <- c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand")

b <- left_join(b, dplyr::select(a, tacc = accession, tchr = chrname, trank = krank, toffset=offset), by = "tacc") 
b <- left_join(b, dplyr::select(a, qacc = accession, qchr = chrname, qrank = krank, qoffset=offset), by = "qacc") 
b <- mutate(b, tlen=tend-tstart, qlen=qend-qstart)

g <- b %>% mutate(tmid=toffset + ((tstart+tend)/2), qmid=qoffset + ((qstart+qend)/2)) %>% dplyr::filter(qlen>500) %>% ggplot() + geom_segment(aes(x=tstart+toffset, y=qstart+qoffset, xend=tend+toffset, yend=qend+qoffset,)) + theme_bw() + scale_x_continuous(breaks = arrange(b, trank) %>% pull(toffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, trank) %>% pull(tchr) %>% unique()) + scale_y_continuous(breaks = arrange(b, qrank) %>% pull(qoffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, qrank) %>% pull(qchr) %>% unique())
dotplot <- ggplotly(g)
dotplot
htmlwidgets::saveWidget(dotplot, "dotplot.html")

##with colors for strand and displaying strands correctly
g <- b %>% dplyr::filter(qlen>500 | tlen>500) %>% mutate(adjstart = if_else(strand=="+", qstart, qend), adjend=if_else(strand=="+",qend,qstart)) %>% ggplot() + geom_segment(aes(x=tstart+toffset, y=adjstart+qoffset, xend=tend+toffset, yend=adjend+qoffset,color=strand)) + theme_bw() + scale_x_continuous(breaks = arrange(b, trank) %>% pull(toffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, trank) %>% pull(tchr) %>% unique()) + scale_y_continuous(breaks = arrange(b, qrank) %>% pull(qoffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, qrank) %>% pull(qchr) %>% unique())
dotplot <- ggplotly(g)
dotplot
htmlwidgets::saveWidget(dotplot, "dotplot.html")

