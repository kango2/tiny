library(tidyverse)
library(xml2)

##assembly_results.xml
##add <documents> at the begining of the document manually
##add </documents> at the end of the document manually

x <- as_list(read_xml("Downloads/assembly_result-4.xml"))
xmlinfo <- list()
asmreports <- list()
for (i in 1:length(x$documents)) {
  thisdoc <- x$documents[i]$DocumentSummary
  infotags <- c("AssemblyAccession", "LastMajorReleaseAccession", "AssemblyName", "Taxid", "SpeciesName", "Organism", "AssemblyStatus", "WGS", "SubmitterOrganization", "RefSeq_category", "FtpPath_GenBank", "FtpPath_RefSeq", "FtpPath_Assembly_rpt", "FtpPath_Stats_rpt")
  fields <- thisdoc[infotags]
  fields <- unlist(fields)
  fields <- tibble(names=names(fields), values = unlist(fields))
  syn <- unlist(thisdoc$Synonym)
  syn <- tibble(names = names(syn), values = unlist(syn))
  fields <- bind_rows(fields, syn)
  fields <- bind_rows(fields, tibble(names = "properties", values = paste(unlist(thisdoc$PropertyList), collapse = ";")))
  asmreportcols <- c("Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name")
  asmrep <- read_delim(filter(fields, names == "FtpPath_Assembly_rpt") %>% pull(values), delim = "\t", comment = "#", col_names = F)
  colnames(asmrep) <- asmreportcols
  asmrep$docid <- attributes(thisdoc)$uid
  fields <- bind_rows(fields, asmrep %>% group_by(`Sequence-Role`) %>% summarise(values = as.character(sum(`Sequence-Length`))) %>% rename(names = `Sequence-Role`))
  fields$docid <- attributes(thisdoc)$uid
  xmlinfo[[i]] <- fields
  asmreports[[i]] <- asmrep
}

xmlinfo <- bind_rows(xmlinfo)
xmlinfo <- xmlinfo %>% pivot_wider(names_from = "names", values_from="values")
xmlinfo <- xmlinfo %>% 
  mutate(across(`assembled-molecule`:`unplaced-scaffold`, as.integer)) %>% 
  rowwise() %>% 
  mutate(totallen = sum(`assembled-molecule`, `unlocalized-scaffold`, `unplaced-scaffold`, na.rm = T), unplacedpid = (totallen - `assembled-molecule`) * 100 / totallen) %>%
  group_by(Taxid) %>%
  mutate(asmcounts = n())

write_delim(xmlinfo, "Downloads/birdgenomes.xmlinfo.txt", delim = "\t")

asmreports <- bind_rows(asmreports)
write_delim(asmreports, "Downloads/assemblyreports.txt", delim = "\t")
