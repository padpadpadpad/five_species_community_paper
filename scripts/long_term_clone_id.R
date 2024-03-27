# ---------------------------
# Purpose of script: to check the identification of clones of the five species community after >1 year in culture
#
# What this script does:
# 1. Reads in a processed phyloseq object
# 2. Visualises the dataset
# 3. Filters the dataset
# 4. Look at how often the clone was identified correctly
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-03-27
#
# Copyright (c) Daniel Padfield, 2024
#
# ---------------------------
#
# Notes: None
#   
#
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian")
}
# load packages
librarian::shelf(phyloseq, mikemc/speedyseq, david-barnett/microViz, tidyverse)

## ---------------------------

# load in the data ####
ps <- readRDS("data/long_term_assignment_16s.rds")

# Visualise dataset ####

# look at number of reads per sample
hist(sample_sums(ps))
summary(sample_sums(ps))
# Minumum read number is 135 (probably failed)!

# look at number of ASVs per sample
d_diversity <- psmelt(ps) %>%
  janitor::clean_names() %>%
  filter(abundance > 0) %>%
  group_by(sample) %>%
  tally() %>%
  ungroup()
# every sample has lots of ASVs in it! This might be concerning!

# look at prevalence of individual ASVS
d_prevalence <- psmelt(ps) %>%
  janitor::clean_names() %>%
  filter(abundance > 0) %>%
  group_by(otu) %>%
  tally() %>%
  ungroup()
# 425 unique ASVs (we should just have 5)!

# wrangle data ####

# turn into dataframe
d_ps <- psmelt(ps) %>%
  janitor::clean_names()

# get metadata
d_meta <- sample_data(ps)

# look at unique communities used
d_meta$community %>% unique()

# check the number of each isolate sequenced
group_by(d_meta, isolate) %>%
  tally()

# find the unique number of samples
unique(d_ps$sample) %>% length()
# ok there are 169 samples

# calculate the relative abundance of each genus, as each isolate is known to be a different genus
d_ps <- d_ps %>%
  group_by_at(vars(-abundance, -otu, -species)) %>%
  summarise(abundance = sum(abundance)) %>%
  group_by(sample) %>%
  mutate(prop = abundance / sum(abundance)) %>%
  ungroup()

# plot distribution of relative abundances
ggplot(d_ps, aes(prop)) +
  geom_histogram(fill = 'grey', col = 'black') +
  theme_bw()
# most things are very rare!
# a few things are very abundant

# keep the most abundant ASV in each sample
d_ps2 <- d_ps %>%
  group_by(sample) %>%
  slice_max(prop, n = 1) %>%
  ungroup()

ggplot(d_ps2, aes(prop)) +
  geom_boxplot(fill = 'grey', col = 'black') +
  theme_bw()
# couple of samples are not good

# filter any that do not have a relative abundance of the most common thing of >0.95
d_not_sure <- filter(d_ps2, prop < 0.95)
d_ps2 <- filter(d_ps2, prop >= 0.95)

# check what sample 127 looks like
filter(d_ps, sample == 'sample127') %>%
  filter(abundance > 0) %>%
  View()
# sample127 and sample99 we have no clue about - sample 99 failed and sample 127 is 60% Achromobacter and 40% Ochrobactrum

# rename isolate so A = Achromobacter, O = Ochrobactrum, S = Stenotrophomonas, P = Pseudomonas, V = Variovorax
d_ps2 <- d_ps2 %>%
  mutate(isolate = case_when(
    isolate == 'A' ~ 'Achromobacter',
    isolate == 'O' ~ 'Ochrobactrum',
    isolate == 'S' ~ 'Stenotrophomonas',
    isolate == 'P' ~ 'Pseudomonas',
    isolate == 'V' ~ 'Variovorax'
  ))

# compare isolate identification to assigned genus from 16S sequencing ####

# add column if they are the same
d_ps2 <- mutate(d_ps2, same = ifelse(isolate == genus, 1, 0))  
sum(d_ps2$same)

the_same <- filter(d_ps2, same == 1)
not_the_same <- filter(d_ps2, same == 0)

nrow(the_same)/nrow(d_ps2) * 100
# Meaghan ID'ed the correct clone 95.8% of the time

# look for any patterns in incorrectly ID'ed things
select(not_the_same, isolate, genus)

# 2 Achromobacter were ID'ed as Ochrobactrum
# 2 Pseudomonas were ID'ed as Variovorax
# 2 Varioborax were ID'ed as Pseudomonas
# 1 Achromobacter was ID'ed as Pseudomonas

# create a table to export
d_table <- select(d_ps2, isolate, genus) %>%
  group_by(isolate, genus) %>%
  summarise(n = n(),
            .groups = 'drop') %>%
  group_by(isolate) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  arrange(isolate, genus)
