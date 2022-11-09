library(dada2)
packageVersion("dada2")

# Load sequence tables
load("lp2017watercolumn16s_seqtab.rda")
seqtab2017 <- seqtab; rm(seqtab)

load("lp2018watercolumn16s_seqtab.rda")
seqtab2018 <- seqtab; rm(seqtab)

load("lp2019watercolumn16sreseq_concat_seqtab.rda")
seqtab2019 <- seqtab; rm(seqtab)

# Merge sequence tables
seqtab_all <- mergeSequenceTables(seqtab2017, seqtab2018, seqtab2019)
save(seqtab_all, file = "lp2017-2019watercolumn16sreseq2_seqtab_all.rda")

# Remove bimeras
seqtab_nochim <- removeBimeraDenovo(seqtab_all, method = "consensus", multithread = TRUE)
save(seqtab_nochim, file = "lp2017-2019watercolumn16sreseq2_seqtab_nochim.rda")
