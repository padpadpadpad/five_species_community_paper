# ---------------------------
# Purpose of script: Visualise the abundance of different ASVs in the disassembly experiment
#
# What this script does:
# 1. Analyses the prevalence filtered 16S object from the disassembly experiment
# 2. Looks at how dominant the 5 abundant ASVs
# 3.
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-07-11
#
# Copyright (c) Daniel Padfield, 2024
#
# ---------------------------
#
# Notes:
#
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)){
  install.packages("librarian")
}
# load packages
librarian::shelf(phyloseq, david-barnett/microViz, dada2, tidyverse, Biostrings, DECIPHER, patchwork)

## ---------------------------

# load in phyloseq
ps <- readRDS("data/disassembly_16s.rds")

# do prevalence filtering 
# have to be present in at least 5% of samples and have at least 200 reads in total
ps <- tax_filter(ps, min_prevalence = 0.05, min_total_abundance = 200)

sample_sums(ps) %>%
  sort() %>% sum()
ps
# only 12 taxa

# turn into a dataset
d <- psmelt(ps) %>%
  janitor::clean_names()

head(d)

# find out which OTU links to which of our 5 species
seqs <- d$otu %>% unique()

assignment <- assignSpecies(seqs, refFasta = 'data/five_spp_16S.fasta', tryRC = TRUE) %>%
  data.frame() %>%
  filter(!is.na(Genus)) %>%
  rownames_to_column(var = "otu") %>%
  select(otu, five_spp = Genus)

d <- left_join(d, assignment) %>%
  mutate(in_comm = ifelse(is.na(five_spp), "no", 'yes'))

d <- group_by(d, replicate) %>%
  mutate(prop = abundance/sum(abundance)) %>%
  ungroup()

d2 <- group_by(d, replicate, in_comm) %>%
  summarise(mean_prop = sum(prop)) %>%
  ungroup()

# look at the % of reads mapped to those in the community or not in the community
ggplot(d2, aes(in_comm, mean_prop)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point()

d3 <- filter(d, abundance > 0) %>%
  filter(prop > 0.001) %>%
  group_by(genus, five_spp, otu) %>%
  summarise(presence = n()/48,
    mean_prop = mean(prop, na.rm = TRUE),
    .groups = 'drop')

# look at relative abundance curve for each community
d_filt <- d %>%
  group_by(replicate) %>%
  arrange(-abundance) %>%
  mutate(rank = 1:n()) %>% 
  ungroup() %>% 
  filter(prop > 0.001) %>%
  group_by(genus, otu) %>%
  mutate(id = cur_group_id()) %>%
  ungroup() %>%
  mutate(genus2 = paste(genus, id, sep = ' '))

# plot of rank abundance of each species

# 

# set up colours
color_scheme <- c(a = '#e41a1c', o = '#377eb8', p = '#4daf4a', s = '#984ea3', v = '#ff7f00', zother = 'light grey')

# set labels to be used in plots
labs <- c(expression(italic("Achromobacter")~sp.), expression(italic("Ochrobactrum")~sp.), expression(italic("Pseudomonas")~sp.), expression(italic("Stenotrophomonas")~sp.), expression(italic("Variovorax")~sp.), 'other')

# create column for colours
d_filt <- mutate(d_filt, species = case_when(genus == 'Achromobacter' & in_comm == 'yes' ~ 'a',
                                             genus == 'Ochrobactrum' & in_comm == 'yes' ~ 'o',
                                             genus == 'Pseudomonas' & in_comm == 'yes' ~ 'p',
                                             genus == 'Stenotrophomonas' & in_comm == 'yes' ~ 's',
                                             genus == 'Variovorax' & in_comm == 'yes' ~ 'v',
                                             TRUE ~ 'zother'))
dodge=0.2

p1 <- mutate(d_filt, rank2 = ifelse(in_comm=='yes', rank + dodge, rank - dodge),
         rank2 = rank2 + runif(n(), -0.1, 0.1)) %>%
  ggplot(aes(rank2, prop)) +
  geom_line(aes(group = replicate), col = 'light grey', alpha = 0.25) +
  geom_point(aes(col = species), size = 3, shape = 21, fill = 'white') +
  theme_bw() +
  labs(y = 'Relative abundance',
       x = 'Species rank') +
  scale_color_manual('Species', values = color_scheme, labels = labs) +
  scale_x_continuous(breaks = c(1:11)) +
  ylim(c(0,1)) +
  theme(panel.grid.minor = element_blank()) +
  NULL

p1

# 
d_incomm <- filter(d_filt, in_comm == 'yes')
d_nocomm <- filter(d_filt, in_comm == 'no')

# look at relative abundance of each species, and its species rank
p2 <- ggplot(d_incomm, aes(fct_reorder(species, -prop, mean), prop)) +
  geom_point(aes(col = species), shape = 21, fill = 'white', position = position_jitter(width = 0.1), show.legend = FALSE, size = 3) + 
  scale_color_manual('Species', values = color_scheme[1:5], labels = labs[1:5]) +
  ylim(c(0,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = 'OTU',
       y = 'Relative abundance',
       title = '(c)') +
  scale_x_discrete(labels = c(labs[3], labs[1], labs[5], labs[4], labs[2]), guide = guide_axis(n.dodge = 2)) +
  NULL

p2

d_nocomm_names <- group_by(d_nocomm, otu, genus) %>%
  summarise(prop = mean(prop), .groups = 'drop') %>%
  group_by(genus) %>%
  arrange(-prop) %>%
  mutate(order = 1:n(),
         species = paste(genus, order, sep = ' ')) %>%
  ungroup() 

d_nocomm <- left_join(select(d_nocomm, -species), select(d_nocomm_names, otu, species))

p3 <- d_nocomm %>%
  ggplot(aes(fct_reorder(species, -prop, mean), prop)) +
  geom_point(shape = 21, fill = 'white', col =  'light grey', position = position_jitter(width = 0.1), show.legend = FALSE, size = 3) +
  ylim(c(0,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = 'OTU',
       y = 'Relative abundance',
       title = '(c)') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  NULL

p3

# visualise all of these ASVs for Pseudomonas
p4 <- bind_rows(d_nocomm, d_incomm) %>%
  filter(genus == 'Pseudomonas') %>%
  arrange(replicate, species) %>%
  ggplot(aes(species, prop)) +
  geom_line(aes(group = replicate), col = 'light grey', alpha = 0.2, position = position_jitter(width = 0.1, seed = 42)) +
  geom_point(aes(col = species, group = replicate), shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1, seed = 42), show.legend = FALSE) +
  theme_bw() +
  labs(title = '(d)',
       y = 'Relative abundance',
       x = 'OTU') +
  scale_x_discrete(labels = c(expression(italic("Pseudomonas")~sp.), 'Pseudomonas 1', 'Pseudomonas 2')) +
  scale_color_manual(values = c('#4daf4a', 'light grey', 'light grey')) +
  ylim(c(0,1)) +
  NULL  

p1 + p2 + p3 + p4 + plot_layout(guides = 'collect')

ggsave('plots/plot_16s.png', last_plot(), height = 8, width = 13)
