# ---------------------------
# Purpose of script: Run QC on the sanger sequencing of the 5 species community morphs
#
# What this script does:
# 1. installs an old version of sangeranalyseR
# 2. plot chromatogram for each file
# 3. trim files and resave in separate folder
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-03-28
#
# Copyright (c) Daniel Padfield, 2024
#
# ---------------------------
#
# Notes:
#
# ---------------------------

# install an old version of sangeranalyseR to allow my old clunky code to work
remotes::install_github('roblanf/sangeranalyseR@0e658db5d0cac2430c79b83b1e2f3f0da435e65f')

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)){
  install.packages("librarian")
}

# load packages
librarian::shelf(DECIPHER, sangerseqR, ShortRead, Biostrings, stringr, sangeranalyseR, tidyverse)

## ---------------------------

#---------------------#
# custom functions ####
#---------------------#

# function for all chromatograms
get_chromatogram <- function(seq_file, output_folder){
  temp <- sangerseqR::read.abif(seq_file)
  seq_sang <- sangerseqR::sangerseq(temp)
  sangeranalyseR::secondary.peaks(seq_sang, output.folder = output_folder, file.prefix = tools::file_path_sans_ext(basename(seq_file)))
}

# trim and resave trimmed sequence
trim_and_save <- function(seq_file, output_path, trim_cutoff = 1e-04){
  temp <- sangerseqR::read.abif(seq_file)
  trims <- sangeranalyseR::trim.mott(temp, cutoff = trim_cutoff)
  seq <- substring(temp@data$PBAS.2, trims$start, trims$finish)
  write(seq, paste(output_path, '/', tools::file_path_sans_ext(basename(seq_file)), '.txt', sep = ''))
}

# looking at the .ab1 files - trim and clean some of these files, then save the sequence as a .seq file

# output folder figs ####
fig_path <- 'plots/chromatogram'

# output data folder
trimmed_path <- 'data/sanger/trimmed'

# list files in the data folder - .ab1 files
files <- list.files('data/sanger', pattern = '.ab1', full.names = TRUE)

# look at how many of each morph we have
morphs <- tibble(file = basename(files)) %>%
  separate(., file, c('blah1', 'blah2', 'morph', 'blah3'), sep = '_') %>%
  select(morph) %>%
  group_by(morph) %>%
  tally()
morphs

# read in a single sequence
seq_abif = read.abif(files[3])
seq_sang <- sangerseq(seq_abif)

# filter out files for morphs used in the 5 species community. This was community 18.
# B - gamma, g
# C - delta, d
# D - epsilon, e
# E - zeta, z
# F - kappa, k

# get all chromatograms ####
map(files, get_chromatogram, output_folder = fig_path)

# remove all .csv files from the fig_path folder
file.remove(list.files(fig_path, pattern = '.csv', full.names = TRUE))

# trim each file, save sequence out as the raw sequence ####

# get summary data for each file
file_sum <- summarise.abi.folder('data/sanger')

write.csv(file_sum$summaries, 'data/sanger_seq_qualcheck.csv')

# do trimming and re-save files
walk(files, trim_and_save, output_path = trimmed_path, trim_cutoff = 1e-04)
