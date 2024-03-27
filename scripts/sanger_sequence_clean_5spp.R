#----------------------------------------------------------#
# process, clean and visualise raw Sanger sequencing files #
#----------------------------------------------------------#

# looking at the .ab1 files - trim and clean some of these files, then save the sequence as a .seq file

# load packages
# install an old version of sangeranalyseR to allow my old clunky code to work
remotes::install_github('roblanf/sangeranalyseR@0e658db5d0cac2430c79b83b1e2f3f0da435e65f')

librarian::shelf(DECIPHER, sangerseqR, ShortRead, msa, annotate, ggtree, Biostrings, stringr, sangeranalyseR, tidyverse)

# Primers used ####
# 515F - GTGYCAGCMGCCGCGGTAA
# 806R - GGACTACNVGGGTWTCTAAT
fwd_515F <- 'GTGYCAGCMGCCGCGGTAA'
rev_806R <- 'GGACTACNVGGGTWTCTAAT'

#---------------------#
# Custom functions ####
#---------------------#

# function for all chromatograms
get_chromatogram <- function(seq_file, output_folder){
  temp <- sangerseqR::read.abif(seq_file)
  seq_sang <- sangerseqR::sangerseq(temp)
  sangeranalyseR::secondary.peaks(seq_sang, output.folder = output_folder, file.prefix = tools::file_path_sans_ext(basename(seq_file)))
}

# check whether primer is present in each sequence
primer_presence <- function(seq_file, fwd_primer_seq, rev_primer_seq){
  
  # create reverse complement sequence
  rev_primer_seq <- chartr('ATGC', 'TACG', rev_primer_seq)
  
  # load in .ab1 file and get sequence
  temp <- sangerseqR::read.scf(seq_file)
  seq_sang <- sangerseqR::sangerseq(temp)
  seq <- temp@data$PBAS.2
  
  # create dataframe
  d <- data.frame(name = tools::file_path_sans_ext(basename(seq_file)), 
                  fwd_primer_present = stringr::str_detect(fwd_primer_seq, seq),
                  rev_primer_present = stringr::str_detect(rev_primer_seq, seq), 
                  stringsAsFactors = FALSE)
  return(d)
  
}

# trim and resave trimmed sequence
trim_and_save <- function(seq_file, output_path, trim_cutoff = 1e-04){
  temp <- sangerseqR::read.abif(seq_file)
  trims <- sangeranalyseR::trim.mott(temp, cutoff = trim_cutoff)
  seq <- substring(temp@data$PBAS.2, trims$start, trims$finish)
  write(seq, paste(output_path, '/', tools::file_path_sans_ext(basename(seq_file)), '.txt', sep = ''))
}

# output folder figs ####
fig_path <- 'sequencing/sanger_final/chromatogram'

# output data folder
trimmed_path <- 'sequencing/sanger_final/trimmed'

# list files in the data folder - .ab1 files
files <- list.files('sequencing/sanger_final/raw', pattern = '.ab1', full.names = TRUE)

morphs <- tibble(file = basename(files)) %>%
  separate(., file, c('blah1', 'blah2', 'morph', 'blah3'), sep = '_') %>%
  select(morph) %>%
  group_by(morph) %>%
  tally()

# read in a single sequence
seq_abif = read.abif(files[3])
seq_sang <- sangerseq(seq_abif)

# filter out files for morphs used in the 5 species community. This was community 18.
# B - gamma
# C - delta
# D - epsilon
# E - zeta
# F - kappa

# to keep
files_to_keep <- c('_g_', '_d_', '_e_', '_z_', '_k_')

files2 <- stringr::str_subset(files, paste(files_to_keep, collapse = '|'))

morphs <- tibble(file = basename(files2)) %>%
  separate(., file, c('blah1', 'blah2', 'morph', 'blah3'), sep = '_') %>%
  select(morph) %>%
  group_by(morph) %>%
  tally()

# get all chromatograms - already run
# map(files, get_chromatogram, output_folder = NA)

# trim each file, save sequence out as the raw sequence ####

# get summary data for each file
file_sum <- summarise.abi.folder('sequencing/sanger_final/raw')
write.csv(file_sum$summaries, 'sequencing/sanger_final/sanger_seq_qualcheck.csv')
# the trimming parameters will improve the quality of the files

# do trimming and re-save files
walk(files2, trim_and_save, output_path = trimmed_path, trim_cutoff = 1e-04)
