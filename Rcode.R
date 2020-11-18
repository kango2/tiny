
a <- read_delim("../micro/chrnames.txt", "\t")
a <- a %>% group_by(assemblyid) %>% mutate(krank=row_number(), offset = lag(length), offset = replace_na(offset, 0), offset = cumsum(offset)) %>% ungroup()


files <- list.files(".", "GCF_000002315.6_GRCg6a.vs.*")

#b %>% mutate(tmid=toffset + ((tstart+tend)/2), qmid=qoffset + ((qstart+qend)/2)) %>% dplyr::filter(qlen>500) %>% ggplot(aes(x=tmid, y=qmid)) + geom_point( size = 0.3) + theme_bw() + scale_x_continuous(breaks = arrange(b, trank) %>% pull(toffset) %>% unique(), labels = arrange(b, trank) %>% pull(tchr) %>% unique()) + scale_y_continuous(breaks = arrange(b, qrank) %>% pull(qoffset) %>% unique(), labels = arrange(b, qrank) %>% pull(qchr) %>% unique())
#b %>% mutate(tmid=toffset + ((tstart+tend)/2), qmid=qoffset + ((qstart+qend)/2)) %>% dplyr::filter(qlen>500) %>% ggplot() + geom_segment(aes(x=tstart+toffset, y=qstart+qoffset, xend=tend+toffset, yend=qend+qoffset,)) + theme_bw() + scale_x_continuous(breaks = arrange(b, trank) %>% pull(toffset) %>% unique(), labels = arrange(b, trank) %>% pull(tchr) %>% unique()) + scale_y_continuous(breaks = arrange(b, qrank) %>% pull(qoffset) %>% unique(), labels = arrange(b, qrank) %>% pull(qchr) %>% unique())
#b %>% dplyr::filter(trank>20) %>% mutate(tmid=toffset + ((tstart+tend)/2), qmid=qoffset + ((qstart+qend)/2)) %>% dplyr::filter(qlen>500) %>% ggplot() + geom_segment(aes(x=tstart+toffset, y=qstart+qoffset, xend=tend+toffset, yend=qend+qoffset,)) + theme_bw() + scale_x_continuous(breaks = arrange(b, trank) %>% pull(toffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, trank) %>% pull(tchr) %>% unique()) + scale_y_continuous(breaks = arrange(b, qrank) %>% pull(qoffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, qrank) %>% pull(qchr) %>% unique())

b <- read_delim(files[4], "\t",col_names = F)
colnames(b) <- c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand")

b <- left_join(b, dplyr::select(a, tacc = accession, tchr = chrname, trank = krank, toffset=offset), by = "tacc") 
b <- left_join(b, dplyr::select(a, qacc = accession, qchr = chrname, qrank = krank, qoffset=offset), by = "qacc") 
b <- mutate(b, tlen=tend-tstart, qlen=qend-qstart)

b %>% mutate(tmid=toffset + ((tstart+tend)/2), qmid=qoffset + ((qstart+qend)/2)) %>% dplyr::filter(qlen>500) %>% ggplot() + geom_segment(aes(x=tstart+toffset, y=qstart+qoffset, xend=tend+toffset, yend=qend+qoffset,)) + theme_bw() + scale_x_continuous(breaks = arrange(b, trank) %>% pull(toffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, trank) %>% pull(tchr) %>% unique()) + scale_y_continuous(breaks = arrange(b, qrank) %>% pull(qoffset) %>% unique(), minor_breaks = NULL, labels = arrange(b, qrank) %>% pull(qchr) %>% unique())
