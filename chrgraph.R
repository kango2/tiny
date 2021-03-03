library(tidyverse)
library(ggplotlyExtra)
library(plotly)

chrinfo <- read_delim("tiny/metadata/chrinfo.2.txt", delim = "\t")
chrinfo <- group_by(chrinfo, speciescode) %>% arrange(krank) %>% mutate(offset = lag(seqlength), offset = replace_na(offset, 0), offset = cumsum(offset), xoffset = offset + ((krank - 1) * 2000000))
dnazoo <- read_delim("tiny/metadata/dnazoo.scfinfo.txt", delim = "\t")
dnazoo <- group_by(dnazoo, speciescode) %>% arrange(krank) %>% mutate(offset = lag(seqlength), offset = replace_na(offset, 0), offset = cumsum(offset), xoffset = offset + ((krank - 1) * 2000000))

chrinfo <- bind_rows(dnazoo, select(chrinfo, colnames(dnazoo)))

speciesinfo <- read_delim("tiny/metadata/species.txt", delim = "\t")
speciesinfo <- speciesinfo %>% filter(figureorder > 0) %>% arrange(desc(figureorder)) %>% mutate(rphyloorder = row_number())
chrinfo <- left_join(chrinfo, select(speciesinfo, uniprotspeciescode, rphyloorder), by = c("speciescode"="uniprotspeciescode"))
chrinfo <- filter(chrinfo, ! is.na(rphyloorder))

files <- list.files("tiny/data/", ".tab")
minalnlen <- 100000
mergedaln <- list()
for (i in 1:nrow(speciesinfo)){
  j <- i + 1
  qcode <- dplyr::filter(speciesinfo, rphyloorder == i) %>% pull(uniprotspeciescode)
  tcode <- dplyr::filter(speciesinfo, rphyloorder == j) %>% pull(uniprotspeciescode)
  f1 <- paste("tiny/data/",qcode,".vs.",tcode,".target.chainpairs.tab", sep = "")
  f2 <- paste("tiny/data/",tcode,".vs.",qcode,".query.chainpairs.tab", sep = "")
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
  b <- left_join(b, dplyr::select(chrinfo, tacc = seqid, tchr = chrname, trank = krank, toffset=xoffset, torder = rphyloorder, tchrtype = chrtype, tcode = speciescode), by = "tacc") 
  b <- left_join(b, dplyr::select(chrinfo, qacc = seqid, qchr = chrname, qrank = krank, qoffset=xoffset, qorder = rphyloorder, qchrtype = chrtype, qcode = speciescode), by = "qacc")
  mergedaln[[i]] <- b 
}
b <- bind_rows(mergedaln)
c <- mutate(b, tlen=tend-tstart, qlen=qend-qstart) %>% 
  filter((qlen >= minalnlen | tlen >= minalnlen) & ! is.na(qchr) & ! is.na(tchr)) %>% 
  mutate(alnid = paste("aln", row_number(), sep = "")) %>%
  mutate( alntype = paste(tchrtype, qchrtype, sep = ":"), x1 = tstart, x2 = tend, x1 = x1 + toffset, x2 = x2 + toffset, x3 = if_else(strand == '-', qstart + qoffset, qend + qoffset), x4 = if_else(strand == '-', qend + qoffset, qstart + qoffset), y1 = torder - 0.01, y2 = torder - 0.01, y3 = qorder + 0.16, y4 = qorder + 0.16) %>%
  select(alnid,x1,x2,x3,x4,y1,y2,y3,y4,alntype)

t1 <- select(c, alnid, alntype, x1:x4) %>% pivot_longer(x1:x4, names_to = "ctype", values_to = "xcoords") %>% select(alnid, alntype, xcoords)
t2 <- select(c, alnid, alntype, y1:y4) %>% pivot_longer(y1:y4, names_to = "ctype", values_to = "ycoords") %>% select(alnid, alntype, ycoords)
t <- bind_cols(t1, t2)
colnames(t) <- c("alnid","alntype", "xcoords","tmp1","tmp2", "ycoords")
t <- t %>% mutate( alntype = case_when(alntype == "macro:macro" ~ "macro", alntype == "micro:micro" ~ "micro", alntype == "macro:micro" ~ "macro X micro", alntype == "micro:macro" ~ "macro X micro"))

p <- ggplot() + geom_rect(data = chrinfo, aes(xmin=xoffset, xmax=xoffset+seqlength, ymin=rphyloorder,ymax= rphyloorder+0.15, fill=chrtype), alpha = 0.6) +
  geom_text(data = chrinfo, aes(x=xoffset, y = rphyloorder + 0.075, label = chrname), hjust = 0, size = 2) + 
  geom_polygon(data = t, aes(x=xcoords,y=ycoords, fill = alntype, group = alnid), alpha = 0.4) +
  scale_y_continuous(breaks = 1:nrow(speciesinfo), labels = arrange(speciesinfo, rphyloorder) %>% pull("Common name")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        text = element_text(size = 16), 
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
ggplotly(p)





###create synteny blocks

b %>% select(1:7, tcode, qcode) %>% 
  group_by(tcode, qcode, tacc) %>% 
  add_tally() %>% 
  dplyr::arrange(tstart, .by_group = TRUE) %>%
  mutate(nt = lead(tacc), qt = lead(qacc), id = if_else((nt == tacc & qt == qacc), 1, 0))
  ##check if current [qt]acc is same as next [qt]acc

