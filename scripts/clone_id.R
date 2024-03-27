# look at the 16S sequencing of the clones after 1 year in culture

# 1. load in necessary packages
librarian::shelf(phyloseq, mikemc/speedyseq, david-barnett/microViz, tidyverse)

# 2. load in the data
ps <- readRDS("data/run_merge/output/ps_no_tree.rds")

# 3. Visualise dataset

# look at number of reads per sample
hist(sample_sums(ps))
summary(sample_sums(ps))

# look at number of ASVs per sample
d_diversity <- psmelt(ps) %>%
  janitor::clean_names() %>%
  filter(abundance > 0) %>%
  group_by(sample) %>%
  tally() %>%
  ungroup()

# look at prevalence of individual ASVS
d_prevalence <- psmelt(ps) %>%
  janitor::clean_names() %>%
  filter(abundance > 0) %>%
  group_by(otu) %>%
  tally() %>%
  ungroup()

# 4. melt data
d_ps <- psmelt(ps) %>%
  janitor::clean_names()

d_meta <- sample_data(ps)
d_meta$community %>% unique()

group_by(d_meta, isolate) %>%
  tally()

# 5. calculate the relative abundance of each genus
d_ps <- d_ps %>%
  group_by_at(vars(-abundance, -otu, -species)) %>%
  summarise(abundance = sum(abundance)) %>%
  group_by(sample) %>%
  mutate(prop = abundance / sum(abundance)) %>%
  ungroup()

# find the unique number of samples
unique(d_ps$sample) %>% length()
# ok there are 169 samples

# plot distribution of relative abundances
ggplot(d_ps, aes(prop)) +
  geom_histogram(fill = 'grey', col = 'black') +
  theme_bw()
# most things are very rare!
# a few things are very abundant

# 6. keep the most abundant ASVs that make up 95% of abundance in each sample
d_ps2 <- d_ps %>%
  group_by(sample) %>%
  slice_max(prop, n = 1) %>%
  ungroup()

ggplot(d_ps2, aes(prop)) +
  geom_boxplot(fill = 'grey', col = 'black') +
  theme_bw()
# couple of samples are not good

# filter any that do not have a relative abundance of the most common thing of >0.88
d_not_sure <- filter(d_ps2, prop < 0.95)
d_ps2 <- filter(d_ps2, prop >= 0.95)

filter(d_ps, sample == 'sample127') %>%
  filter(abundance > 0) %>%
  View()

# sample127 and sample99 we have no clue about - sample 99 failed and sample 127 is 50% Achromobacter and 40% Ochrobactrum

# rename isolate so A = Achromobacter, O = Ochrobactrum, S = Stenotrophomonas, P = Pseudomonas, V = Variovorax
d_ps2 <- d_ps2 %>%
  mutate(isolate = case_when(
    isolate == 'A' ~ 'Achromobacter',
    isolate == 'O' ~ 'Ochrobactrum',
    isolate == 'S' ~ 'Stenotrophomonas',
    isolate == 'P' ~ 'Pseudomonas',
    isolate == 'V' ~ 'Variovorax'
  ))

d_ps2 <- mutate(d_ps2, same = ifelse(isolate == genus, 1, 0))  

sum(d_ps2$same)

the_same <- filter(d_ps2, same == 1)
not_the_same <- filter(d_ps2, same == 0)

160/167 * 100

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

write.csv(d_table, 'data/clone_id.csv', row.names = FALSE)
