# Filter ASVs

# Load libraries
library(tidyverse)
library(seqinr)


#### Import ASV data ####
# Import melted sequence table
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn16sreseq2_melt_all.tsv", col_names = TRUE)


#### Examine ASV distributions ####
asv_nlakes <- asv_melt %>%
  group_by(asv_code, sequence) %>%
  dplyr::count(name = "nlakes") %>%
  ungroup()


#### Calculate ASV sequence counts ####
asv_nseqs <- asv_melt %>%
  group_by(asv_code, sequence) %>%
  summarize(nseqs = sum(nseqs))

asv_nlakes %>%
  left_join(asv_nseqs, by = c("asv_code", "sequence")) %>%
  mutate(nseqs_category = case_when(nseqs == 1 ~ "1 seq",
                                    nseqs == 2 ~ "2 seqs",
                                    nseqs > 2 & nseqs <= 10 ~ "3 - 10 seqs",
                                    nseqs > 10 ~ ">10 seqs")) %>%
  ggplot() +
  geom_histogram(aes(x = nlakes, y = ..count.., fill = nseqs_category),
                 binwidth = 1) +
  theme_classic()


#### Import alignment ####
alignment <- read.fasta("output/fasta/lp2017-2019watercolumn16sreseq2_all_sina.fasta", "DNA", as.string = TRUE, forceDNAtolower = FALSE)
alignment <- tibble(asv_code = names(alignment),
                    sequence = unlist(getSequence(alignment, as.string = TRUE)))


#### Trim poorly aligned positions ####
asv_start_position <- 182L
asv_end_position <- 1400L

asv_trim <- alignment %>%
  mutate(sequence = str_replace_all(sequence, "U", "T")) %>%
  mutate(trimmed = str_sub(sequence, start = asv_start_position, end = asv_end_position)) %>%
  mutate(trimmed = str_replace_all(trimmed, "-", "")) %>%
  mutate(seq_length = nchar(trimmed)) %>%
  mutate(nt_outside_pos = case_when(grepl("A|T|C|G", str_sub(sequence, 1, asv_start_position - 1)) |
                                      grepl("A|T|C|G", str_sub(sequence, asv_end_position + 1, nchar(sequence))) ~ TRUE,
                                    TRUE ~ FALSE))

# Assess alignment of nucleotides outside positions
asv_trim %>%
  group_by(nt_outside_pos) %>%
  dplyr::count(name = "nasvs")

# Assess sequence space occupied by ASVs with nucleotides outside positions
asv_trim %>%
  filter(nt_outside_pos) %>%
  left_join(asv_nseqs, by = "asv_code") %>%
  arrange(-nseqs)

# Assess fragment length lower and upper limits
asv_trim %>%
  select(asv_code, seq_length, nt_outside_pos) %>%
  ggplot() +
  geom_histogram(aes(x = seq_length, y = ..count.., fill = nt_outside_pos),
                 binwidth = 1) +
  theme_classic()

asv_trim %>%
  filter(!nt_outside_pos) %>%
  group_by(seq_length) %>%
  dplyr::count()

# Assign minimum and maximum sequence lengths
seq_length_min <- 251L
seq_length_max <- 257L

# Assess sequence loss associated with removing ASVs
asv_trim %>%
  filter(nt_outside_pos | seq_length < seq_length_min | seq_length > seq_length_max) %>%
  select(asv_code, seq_length) %>%
  left_join(asv_nseqs, by = "asv_code") %>%
  summarize(nseqs = sum(nseqs)) %>%
  mutate(pct_seqs = nseqs/sum(asv_melt$nseqs))

# Remove ASVs with sequence lengths outside lower and upper limits
asv_rmpos <- asv_trim %>%
  filter(!nt_outside_pos) %>%
  filter(seq_length >= seq_length_min & seq_length <= seq_length_max) %>%
  select(asv_code, trimmed) %>%
  arrange(asv_code)


#### Write new files for retained ASVs ####
# Write new fasta file for retained ASVs
# write.fasta(sequences = as.list(asv_rmpos$trimmed),
#             names = asv_rmpos$asv_code,
#             file.out = "output/fasta/lp2017-2019watercolumn16sreseq2_rmpos.fasta")

# Write new melted ASV table for retained ASVs
asv_melt_rmpos <- asv_melt %>%
  filter(asv_code %in% asv_rmpos$asv_code)
# asv_melt_rmpos %>%
#   write_tsv("output/melted/lp2017-2019watercolumn16sreseq2_melt_rmpos.tsv")
