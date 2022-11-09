library(tidyverse)

#### Import and format data ####
sampling_years <- read_tsv("data/lakepulse_sampling_years.tsv")

filepath_nreads <- "output/nreads/"

(files_nreads <- dir(path = filepath_nreads, pattern = "*_nreads.txt"))

nreads <- tibble(filename = files_nreads) %>%
  filter(grepl("^\\d\\d\\-\\d\\d\\d", filename)) %>%
  mutate(sample_id = str_remove_all(filename, "_.*"),
         read_direction = case_when(grepl("_R1_", filename) ~ "fwd",
                                    grepl("_R2_", filename) ~ "rev"),
         processing_step = case_when(grepl("trimmed", filename) ~ "trimmed",
                                     grepl("filtered", filename) ~ "filtered",
                                     TRUE ~ "raw")) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(filepath_nreads, .),
                                                  col_names = c("nreads")))) %>%
  unnest(c(file_contents)) %>%
  dplyr::select(-filename)

nreads <- nreads %>%
  unite("processing_step", c(processing_step, read_direction), sep = "_") %>%
  pivot_wider(names_from = processing_step, values_from = nreads) %>%
  dplyr::select(sample_id, raw_fwd, raw_rev, trimmed_fwd, trimmed_rev, filtered_fwd, filtered_rev)

readtracking2017 <- read_tsv("output/read_tracking/lp2017watercolumn16s_dada2_readtracking.tsv") %>%
  rename(sample_id = ...1)

readtracking2018 <- read_tsv("output/read_tracking/lp2018watercolumn16s_dada2_readtracking.tsv") %>%
  rename(sample_id = ...1)

readtracking2019 <- read_tsv("output/read_tracking/lp2019watercolumn16sreseq_concat_dada2_readtracking.tsv") %>%
  rename(sample_id = ...1)

readtracking_all <- bind_rows(readtracking2017,
                              readtracking2018,
                              readtracking2019)

load("output/sequence_tables/lp2017-2019watercolumn16sreseq2_seqtab_nochim.rda")
seqtab_nochim <- rowSums(seqtab_nochim) %>%
  enframe(name = "sample_id", value = "nochim")

nreads_all <- nreads %>%
  left_join(readtracking_all, by = "sample_id") %>%
  left_join(seqtab_nochim, by = "sample_id") %>%
  left_join(sampling_years, by = c("sample_id" = "lakepulse_id")) %>%
  relocate(sampling_year, .after = sample_id)

# Write n reads summary to file
# nreads_all %>%
#   write_tsv("output/read_tracking/lp2017-2019watercolumn16sreseq2_readtracking.tsv")

# Import n reads summary
nreads_all <- read_tsv("output/read_tracking/lp2017-2019watercolumn16sreseq2_readtracking.tsv")


#### Visualize read loss ####
nreads_long <- nreads_all %>%
  dplyr::select(sample_id, sampling_year, ends_with("_fwd"), merged, nochim) %>%
  rename_with(~str_remove(.x, "_.*"), ends_with("_fwd")) %>%
  pivot_longer(!c(sample_id, sampling_year), names_to = "processing_step", values_to = "nreads")

(nreads_lineplot <- nreads_long %>%
    ggplot(aes(x = factor(processing_step, levels = c("raw", "trimmed", "filtered", "dada", "merged", "nochim")),
               y = nreads, group = sample_id)) +
    facet_wrap(~sampling_year) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(labels = scales::comma, breaks = seq(0, max(nreads_long$nreads), by = 5e4)) +
    labs(y = "Number of sequences") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))
#ggsave("figures/lp2017-2019watercolumn16sreseq2_readtracking_lineplot.pdf", nreads_lineplot, width = 12, height = 6)

(nreads2019_lineplot <- nreads_long %>%
    filter(sampling_year == 2019) %>%
    ggplot(aes(x = factor(processing_step, levels = c("raw", "trimmed", "filtered", "dada", "merged", "nochim")),
               y = nreads, group = sample_id)) +
    facet_wrap(~sampling_year) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(labels = scales::comma, breaks = seq(0, max(nreads_long$nreads), by = 5e3)) +
    labs(y = "Number of sequences") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))
#ggsave("figures/lp2019watercolumn16sreseq_concat_readtracking_lineplot.pdf", nreads2019_lineplot, width = 5, height = 6)

(nreads2017_lineplot <- nreads_long %>%
    filter(sampling_year == 2017) %>%
    ggplot(aes(x = factor(processing_step, levels = c("raw", "trimmed", "filtered", "dada", "merged", "nochim")),
               y = nreads, group = sample_id)) +
    facet_wrap(~sampling_year) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(labels = scales::comma, breaks = seq(0, max(nreads_long$nreads), by = 1e5)) +
    labs(y = "Number of sequences") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))
#ggsave("figures/lp2017watercolumn16sreseq_concat_readtracking_lineplot.pdf", nreads2017_lineplot, width = 5, height = 6)

(nreads2018_lineplot <- nreads_long %>%
    filter(sampling_year == 2018) %>%
    ggplot(aes(x = factor(processing_step, levels = c("raw", "trimmed", "filtered", "dada", "merged", "nochim")),
               y = nreads, group = sample_id)) +
    facet_wrap(~sampling_year) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(labels = scales::comma, breaks = seq(0, max(nreads_long$nreads), by = 1e5)) +
    labs(y = "Number of sequences") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))
#ggsave("figures/lp2018watercolumn16sreseq_concat_readtracking_lineplot.pdf", nreads2018_lineplot, width = 5, height = 6)

(nreads_boxplot <- nreads_long %>%
    ggplot(aes(x = factor(processing_step, levels = c("raw", "trimmed", "filtered", "dada", "merged", "nochim")),
               y = nreads,
               colour = factor(processing_step, levels = c("raw", "trimmed", "filtered", "dada", "merged", "nochim")))) +
    facet_wrap(~sampling_year) +
    geom_jitter(alpha = 0.5, height = 0) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::comma, breaks = seq(0, max(nreads_long$nreads), by = 5e4)) +
    labs(y = "Number of sequences") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          legend.position = "none"))
#ggsave("figures/lp2017-2019watercolumn16sreseq2_readtracking_boxplot.pdf", nreads_boxplot, width = 12, height = 6)


#### Count ASVs ####
load("output/sequence_tables/lp2017-2019watercolumn16sreseq2_seqtab_nochim.rda")

seqtab_long <- seqtab_nochim %>%
  as_tibble(rownames = "sample_id") %>%
  filter(grepl("^\\d\\d\\-\\d\\d\\d", sample_id)) %>%
  pivot_longer(!sample_id, names_to = "sequence", values_to = "ncontigs") %>%
  filter(ncontigs > 0)

nasvs <- seqtab_long %>%
  group_by(sample_id) %>%
  count(name = "nasvs") %>%
  ungroup() %>%
  left_join(sampling_years, by = c("sample_id" = "lakepulse_id")) %>%
  relocate(sampling_year, .after = sample_id)

(nasvs_histogram <- nasvs %>%
    ggplot() +
    geom_histogram(aes(x = nasvs, y = ..count.., fill = as.character(sampling_year)), binwidth = 50) +
    scale_y_continuous(breaks = seq(0, nrow(nasvs), by = 10)) +
    labs(x = "Number of ASVs",
         y = "Number of samples",
         fill = "Sampling year") +
    theme_bw())
#ggsave("figures/lp2017-2019watercolumn16sreseq2_nasvs_histogram.pdf", nasvs_histogram, width = 5, height = 4)
