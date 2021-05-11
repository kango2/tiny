library(tidyverse)

chrinfo <- read_delim("chrinfo.2.txt", delim = "\t")
chrinfo <- group_by(chrinfo, speciescode) %>% arrange(krank) %>% mutate(offset = lag(seqlength), offset = replace_na(offset, 0), offset = cumsum(offset), xoffset = offset + ((krank - 1) * 2000000))
dnazoo <- read_delim("dnazoo.scfinfo.txt", delim = "\t")
dnazoo <- group_by(dnazoo, speciescode) %>% arrange(krank) %>% mutate(offset = lag(seqlength), offset = replace_na(offset, 0), offset = cumsum(offset), xoffset = offset + ((krank - 1) * 2000000))
chrinfo <- bind_rows(dnazoo, select(chrinfo, colnames(dnazoo)))
speciesinfo <- read_delim("species.txt", delim = "\t")

##get proportions of micro and macro chromosomes
write_delim(left_join(chrinfo %>% group_by(speciescode, chrtype) %>% summarise(l = sum(seqlength)) %>% pivot_wider(names_from = chrtype, values_from = l), filter(speciesinfo, in_paper>0) %>% select(`Common name`, speciescode = uniprotspeciescode, genomelength, in_paper)) %>% filter(!is.na(genomelength)) %>% arrange(in_paper) %>% mutate(microp = micro/genomelength, macrop = macro/genomelength), "chrproportions.txt", delim = "\t")
