# ---------------------------
# Purpose of script: To create a summary table of the five species' genomes
#
# What this script does:
# 1. Create summary table of the five species' genomes
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
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)){
  install.packages("librarian")
}
# load packages
librarian::shelf(flextable, officer, tidyverse)

## ---------------------------

# load in data from assigning taxonomy from the 16S sequences present in the genome
d_16s <- read.csv('data/assignment_16s.csv') %>%
  group_by(file)%>%
  # keep only the best hit
  slice_min(., dist_from_ref_16s) %>%
  mutate(species = replace_na(species, ''),
         assignment_16s = trimws(paste(genus, species, sep = ' '))) %>%
  select(sample = file, assignment_16s, num_16s) 

# read in summary data from GTDBtk
d_gtdbtk <- read.csv('data/gtdb_short.csv')

# read in summary data from checkm2 and checkm
d_checkm2 <- read.csv('data/checkm2.csv')
d_checkm <- read.csv('data/checkm.csv')

# combine these datasets together
d_taxa <- left_join(d_gtdbtk, d_16s) %>%
  left_join(., d_checkm2) %>%
  left_join(., select(d_checkm, sample, n_contigs)) %>%
  arrange(sample) %>%
  # add in column for ncbi strain names
  mutate(community_member = 1:n(),
         species2 = c('Achromobacter veterisilvae', 'Ochrobactrum teleogrylli', 'Pseudomonas fluorescens', 'Stenotrophomonas', 'Variovorax'),
         species3 = c('AB1', 'AB1', 'AB1', 'sp. AB1 (2024)', 'sp. AB1 (2024)'))

# select only columns we want to display
d_table <- select(d_taxa, community_member, sample, species, species2, species3, num_16s, genome_size, n_contigs, GC, total_coding_sequences, completeness, contamination) %>%
  mutate(across(genome_size:contamination, \(x) round(x, 2)))

# make table in flextable
table1 <- flextable(d_table, col_keys = c('community_member', 'dummy', 'sample', 'species', 'num_16s', 'genome_size', 'n_contigs', 'GC', 'total_coding_sequences', 'completeness', 'contamination')) %>%
  set_header_labels(community_member = 'Isolate',
                    sample = 'Sanger sequence\nassignment',
                    species = 'GTDBtk\nassignment',
                    dummy = 'NCBI submission\nname',
                    num_16s = '16s copy\nnumber',
                    genome_size = 'Genome\nsize (bp)',
                    n_contigs = 'Number\nof contigs',
                    GC = 'GC\ncontent',
                    total_coding_sequences = 'Total\ncoding\nsequences',
                    completeness = 'Completeness',
                    contamination = 'Contamination') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  fix_border_issues() %>%
  italic(j = c('species')) %>%
  mk_par(j = "dummy", value = as_paragraph(as_i(species2), ' ', (species3))) %>%
  autofit(add_w = 0.0) %>%
  bg(bg = 'white', part = 'all') 

table1

# save out table
save_as_image(table1, 'plots/table_1.png', webshot = 'webshot2')
save_as_docx(table1, path = 'plots/table_1.docx')
