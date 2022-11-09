# Curate samples

# Load libraries
library(tidyverse)
library(seqinr)


#### Import ASV data ####
# Import melted sequence table
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn16sreseq2_melt_rmpos.tsv", col_names = TRUE)


# #### Remove negative control, mock community, and ELA lake samples ####
asv_melt <- asv_melt %>%
  filter(grepl("\\d\\d\\-\\d\\d\\d", sample_id))
length(unique(asv_melt$sample_id))


#### Create bacteria-only dataset ####
# Create taxonomy table
taxonomy <- asv_melt %>%
  distinct(asv_code, .keep_all = TRUE) %>%
  select(asv_code, sequence, kingdom, phylum, class, order, lineage, clade, tribe) %>%
  arrange(asv_code)

# Identify ASVs unassigned at the domain (kingdom) rank
kingdom_na <- taxonomy %>%
  filter(is.na(kingdom)) %>%
  pull(asv_code)
length(kingdom_na)

# Identify ASVs assigned to archaea
archaea <- taxonomy %>%
  filter(kingdom == "Archaea") %>%
  pull(asv_code)
length(archaea)

# Identify ASVs assigned to chloroplasts
chloroplasts <- taxonomy %>%
  filter(order == "Chloroplast") %>%
  pull(asv_code)
length(chloroplasts)

# Identify ASVs assigned to eukaryotes
eukaryota <- taxonomy %>%
  filter(kingdom == "Eukaryota") %>%
  pull(asv_code)
length(eukaryota)

# Calculate number of ASVs to remove
length(c(kingdom_na, archaea, chloroplasts, eukaryota))

# Remove kingdom-level NAs, archaea, chloroplasts, and eukaryotes from melted ASV table
asv_melt <- asv_melt %>%
  filter(!asv_code %in% kingdom_na & !asv_code %in% archaea & !asv_code %in% chloroplasts & !asv_code %in% eukaryota)
length(unique(asv_melt$asv_code))


#### Remove samples with low sequence counts ####
# Plot histogram of sequence counts by sample
asv_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ggplot() +
  geom_histogram(aes(x = nseqs, y = ..count..), binwidth = 1000) +
  scale_x_continuous() +
  theme_bw()

# Identify samples with fewer than 10,000 sequences
lakes_rm <- asv_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs = sum(nseqs)) %>%
  filter(nseqs < 10000) %>%
  pull(sample_id)

# Filter out low-sequence samples
asv_melt <- asv_melt %>%
  filter(!sample_id %in% lakes_rm)
length(unique(asv_melt$sample_id))  # Number of lakes

# Calculate number of retained samples
(nsamples <- length(unique(asv_melt$sample_id)))


#### Calculate relative sequence abundance ####
samples_nseqs <- asv_melt %>%
  group_by(sample_id) %>%
  summarize(sample_nseqs = sum(nseqs))

asv_melt <- asv_melt %>%
  left_join(samples_nseqs, by = "sample_id") %>%
  mutate(relseqs = nseqs/sample_nseqs) %>%
  select(-sample_nseqs)
sum(asv_melt$relseqs) == nsamples  # Should evaluate to TRUE


#### Save curated sample data ####
# Create and write ASV table based on sequence counts
asvtable_nseqs <- asv_melt %>%
  dplyr::select(sample_id, asv_code, nseqs) %>%
  pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0)

asvtable_nseqs <- asvtable_nseqs %>%
  arrange(sample_id) %>%
  dplyr::select(sample_id, order(colnames(asvtable_nseqs)))
sum(asvtable_nseqs[,-1]) == sum(asv_melt$nseqs)  # Should evaluate to TRUE
# asvtable_nseqs %>%
#   write_tsv(paste0("output/tables/lp2017-2019watercolumn16sreseq2_", nsamples, "lakes_asvtable_nseqs.tsv"), col_names = TRUE)

# Create and write new taxonomy table
taxonomy <- asv_melt %>%
  distinct(asv_code, .keep_all = TRUE) %>%
  select(asv_code, sequence, kingdom, phylum, class, order, lineage, clade, tribe) %>%
  arrange(asv_code)
# taxonomy %>%
#   write_tsv(paste0("output/tables/lp2017-2019watercolumn16sreseq2_", nsamples, "lakes_taxonomy.tsv"), col_names = TRUE)

# Create and write combined ASV/taxonomy table
asv_by_site <- asv_melt %>%
  dplyr::select(sample_id, asv_code, nseqs) %>%
  pivot_wider(names_from = sample_id, values_from = nseqs, values_fill = 0)

asv_by_site <- asv_by_site %>%
  arrange(asv_code) %>%
  dplyr::select(asv_code, order(colnames(asv_by_site)))
sum(asv_by_site[,-1]) == sum(asv_melt$nseqs)  # Should evaluate to TRUE

asv_by_site_taxonomy <- taxonomy %>%
  left_join(asv_by_site, by = "asv_code")
asv_by_site_taxonomy %>%
  write_tsv(paste0("output/tables/lp2017-2019watercolumn16sreseq2_", nsamples, "lakes_asvtable_taxonomy.tsv"), col_names = TRUE)

# Write new melted ASV table
# asv_melt %>%
#   write_tsv(paste0("output/melted/lp2017-2019watercolumn16sreseq2_", nsamples, "lakes_melt.tsv"), col_names = TRUE)

# Write new fasta file
# write.fasta(sequences = as.list(taxonomy$sequence),
#             names = taxonomy$asv_code,
#             file.out = paste0("output/fasta/lp2017-2019watercolumn16sreseq2_", nsamples, "lakes.fasta"))
