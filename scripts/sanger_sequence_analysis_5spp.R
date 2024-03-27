# sanger sequence analysis of 5 morphs

# load packages
librarian::shelf(sangeranalyseR, Biostrings, ShortRead, msa, tidyverse, stringr, dada2, ggtree, patchwork, extrafont, flextable, officer)

loadfonts()

#---------------------#
# Custom functions ####
#---------------------#

# load in trimmed files and make into a StringSet from which we can find consensus sequences
read_and_bind <- function(trimmed_file){
  temp <- NULL
  try(temp <- data.frame(file = basename(trimmed_file), seq = read.table(trimmed_file, stringsAsFactors = FALSE)$V1, stringsAsFactors = FALSE), silent = TRUE)
  if(!is.null(temp)) return(temp)
  
}

# a function to return the longest sequence that only contains CATG
clean_base_pairs <- function(sequence, to_contain = "CATG"){
  temp <- stringr::str_split(sequence, paste("[^", to_contain,"]", sep = ''))
  temp <- data.frame(seq = unlist(temp), stringsAsFactors = FALSE)
  temp$length = nchar(temp$seq)
  return(temp[temp$length == max(temp$length),]$seq)
}

#------------------------------------------#
# load in trimmed sanger sequence files ####
#------------------------------------------#

# load in cleaned sanger sequence data
trimmed_files <- list.files('sequencing/sanger_final/trimmed/', full.names = TRUE, pattern = '.txt')

# filter out files for morphs used in the 5 species community. This was community 18.
# B - gamma
# C - delta
# D - epsilon
# E - zeta
# F - kappa

# to keep
files_to_keep <- c('_g_', '_d_', '_e_', '_z_', '_k_')

trimmed_files <- stringr::str_subset(trimmed_files, paste(files_to_keep, collapse = '|'))

# bind all files together
d_trim <- map_df(trimmed_files, read_and_bind) %>%
  mutate(., seq_len = nchar(seq))

#-----------------------------------------------#
# quality control on these trimmed sequences ####
#-----------------------------------------------#

# filter out files that have < 100 base pairs
d_trim <- filter(d_trim, seq_len >= 100)
# 4 have dropped out

# filter out files that have > 5 secondary peaks or an average quality below 30
to_trim <- read_csv('sequencing/sanger_final/sanger_seq_qualcheck.csv') %>% 
  filter(trimmed.secondary.peaks > 5 | trimmed.mean.quality < 30) %>%
  pull(file.name) %>%
  tools::file_path_sans_ext()

d_trim <- filter(d_trim,! tools::file_path_sans_ext(file) %in% to_trim)
# to 43 from 44

# split up file name into morphs to create alignments for each morph
d_trim <- mutate(d_trim, file = gsub('-', '_', file)) %>%
  separate(., file, c('week', 'microcosm', 'morph', 'blah1'), sep = '_', remove = FALSE) %>%
  select(-contains('blah'))

# align sequences for each morph ####

# how many morphs
unique(d_trim$morph)
# "k" "g" "d" "e" "z"

# consensus sequence for morph b - gamma 
d_trim_b <- filter(d_trim, morph == 'g') 
b_SS <- DNAStringSet(d_trim_b$seq)
names(b_SS) <- d_trim_b$file
b_align <- msa(b_SS, method = 'ClustalW')

b_consensus <- msaConsensusSequence(b_align)

b_align_conv <- as.DNAbin(unmasked(b_align))

b_dist <- ape::dist.dna(b_align_conv, model = 'JC69') %>% ape::njs()

d_labs <- tibble(label = b_dist$tip.label) %>%
  separate(label, c('week', 'sample', 'morph', 'fwd'), sep = '_', remove = FALSE) %>%
  mutate(., label2 = paste('microcosm', sample))

p_b <- ggtree(b_dist) %<+% d_labs + 
  geom_tiplab(aes(label=label2), size = MicrobioUoE::pts(9), hjust = -0.1) +
  geom_tippoint(size = 3) +
  labs(title = '(a) morph B')

# consensus sequence for c - delta
d_trim_c <- filter(d_trim, morph == 'd') 
c_SS <- DNAStringSet(d_trim_c$seq)
names(c_SS) <- d_trim_c$file
c_align <- msa(c_SS, method = 'ClustalW')
c_consensus <- msaConsensusSequence(c_align) 

c_align_conv <- as.DNAbin(unmasked(c_align))

c_dist <- ape::dist.dna(c_align_conv, model = 'JC69') %>% njs()

d_labs <- tibble(label = c_dist$tip.label) %>%
  separate(label, c('week', 'sample', 'morph', 'fwd'), sep = '_', remove = FALSE) %>%
  mutate(., label2 = paste('microcosm', sample))

p_c <- ggtree(c_dist) %<+% d_labs + 
  geom_tiplab(aes(label=label2), size = MicrobioUoE::pts(9), hjust = -0.1) +
  geom_tippoint(size = 3) +
  labs(title = '(c) morph C')

# consensus sequence for d - epsilon
d_trim_d <- filter(d_trim, morph == 'e') 
d_SS <- DNAStringSet(d_trim_d$seq)
names(d_SS) <- d_trim_d$file
d_align <- msa(d_SS, method = 'ClustalW')

d_consensus <- msaConsensusSequence(d_align) 

d_align_conv <- as.DNAbin(unmasked(d_align))

d_dist <- ape::dist.dna(d_align_conv, model = 'JC69') %>% njs()

d_labs <- tibble(label = d_dist$tip.label) %>%
  separate(label, c('week', 'sample', 'morph', 'fwd'), sep = '_', remove = FALSE) %>%
  mutate(., label2 = paste('microcosm', sample))

p_d <- ggtree(d_dist) %<+% d_labs + 
  geom_tiplab(aes(label=label2), size = MicrobioUoE::pts(9), hjust = -0.1) +
  geom_tippoint(size = 3) +
  labs(title = '(d) morph D') +
  scale_x_continuous(expand = c(.1, .1))

# need 3 DNA sequences for d
d_trim_d1 <- filter(d_trim_d, file %in% c('3_15_e_515F.txt', '3_39_e_515F.txt', '3_27_e_515F.txt', '6_18_e_515F.txt')) 
d_SS <- DNAStringSet(d_trim_d1$seq)
names(d_SS) <- d_trim_d1$file
d_align <- msa(d_SS, method = 'ClustalW')
d_consensus1 <- msaConsensusSequence(d_align) 
d_trim_d2 <- filter(d_trim_d, file %in% c('3_17_e_515F.txt', '6_9_e_515F.txt', '6_27_e_515F.txt')) 
d_SS <- DNAStringSet(d_trim_d2$seq)
names(d_SS) <- d_trim_d2$file
d_align <- msa(d_SS, method = 'ClustalW')
d_consensus2 <- msaConsensusSequence(d_align) 
d_trim_d3 <- filter(d_trim_d, file %in% c('6_39_e_515F.txt')) 
d_SS <- DNAStringSet(c(d_trim_d3$seq, d_trim_d3$seq))
d_align <- msa(d_SS, method = 'ClustalW')
d_consensus3 <- msaConsensusSequence(d_align)

d_consensus <- c(d_consensus1, d_consensus2, d_consensus3)

# consensus sequence for e - zeta
d_trim_e <- filter(d_trim, morph == 'z') 
e_SS <- DNAStringSet(d_trim_e$seq)
names(e_SS) <- d_trim_e$file
e_align <- msa(e_SS, method = 'ClustalW')

e_consensus <- msaConsensusSequence(e_align) 

e_align_conv <- as.DNAbin(unmasked(e_align))

e_dist <- ape::dist.dna(e_align_conv, model = 'JC69') %>% ape::njs()

d_labs <- tibble(label = e_dist$tip.label) %>%
  separate(label, c('week', 'sample', 'morph', 'fwd'), sep = '_', remove = FALSE) %>%
  mutate(., label2 = paste('microcosm', sample))

p_e <- ggtree(e_dist) %<+% d_labs + 
  geom_tiplab(aes(label=label2), size = MicrobioUoE::pts(9), hjust = -0.1) +
  geom_tippoint(size = 3) +
  labs(title = '(e) morph E') +
  scale_x_continuous(expand = c(.1, .1))

# consensus sequence for f - kappa
d_trim_f <- filter(d_trim, morph == 'k') 
f_SS <- DNAStringSet(d_trim_f$seq)
names(f_SS) <- d_trim_f$file
f_align <- msa(f_SS, method = 'ClustalW')

f_consensus <- msaConsensusSequence(f_align) 

f_align_conv <- as.DNAbin(unmasked(f_align))

f_dist <- ape::dist.dna(f_align_conv, model = 'JC69') %>% ape::njs()

d_labs <- tibble(label = f_dist$tip.label) %>%
  separate(label, c('week', 'sample', 'morph', 'fwd'), sep = '_', remove = FALSE) %>%
  mutate(., label2 = paste('microcosm', sample))

p_f <- ggtree(f_dist) %<+% d_labs + 
  geom_tiplab(aes(label=label2), size = MicrobioUoE::pts(9), hjust = -0.1) +
  geom_tippoint(size = 3) +
  labs(title = '(f) morph F') +
  scale_x_continuous(expand = c(.1, .1))

# need 2 DNA sequences for f
d_trim_f1 <- filter(d_trim_f, file != '6_9_k_515F.txt') 
f_SS <- DNAStringSet(d_trim_f1$seq)
f_align <- msa(f_SS, method = 'ClustalW')
f_consensus1 <- msaConsensusSequence(f_align) 
d_trim_f2 <- filter(d_trim_f, file == '6_9_k_515F.txt') 
f_SS <- DNAStringSet(c(d_trim_f2$seq, d_trim_f2$seq))
f_align <- msa(f_SS, method = 'ClustalW')
f_consensus2 <- msaConsensusSequence(f_align) 

f_consensus <- c(f_consensus1, f_consensus2)

# combine all the consensus sequences together
names <- paste(c('b', 'c', 'd', 'e', 'f'), 'consensus', sep = '_')
all_consensus <- mget(names) %>% unlist(., use.names = FALSE)
morph <- c('b', 'c', 'd', 'd', 'd', 'e', 'f', 'f')

all_consensus <- tibble(morph = morph, seq = all_consensus) %>%
  mutate(., seq = gsub('\\?', 'N', seq),
         seq = gsub('-', 'N', seq)) %>%
  group_by(morph) %>%
  mutate(., id = 1:n()) %>%
  ungroup() %>%
  mutate(., morph = paste(morph, id, sep = ''))

# do a phylogenetic tree of all consensus sequences
all_SS <- DNAStringSet(all_consensus$seq)
names(all_SS) <- all_consensus$morph

tree <- msa(all_SS, method = 'ClustalW') %>%
  as.DNAbin(unmasked(.)) %>%
  ape::dist.dna(., model = 'JC69') %>%
  ape::njs()

d_labs <- tibble(label = tree$tip.label) %>%
  mutate(., label2 = paste('morph ', toupper(substr(label,1,1)),'\nconsensus seq ', substr(label,2,2), sep = ''))

p_all <- ggtree(tree) %<+% d_labs + 
  geom_tiplab(aes(label=label2), size = MicrobioUoE::pts(12), hjust = -0.1) +
  geom_tippoint(size = 3) +
  scale_x_continuous(expand = c(.1, .1))

ggsave('manu_figs/consensus_seq_tree.pdf', p_all, height = 9, width =7)
ggsave('manu_figs/consensus_seq_tree.png', p_all, height = 9, width =7)

# save these files out into fasta files
list_consensus <- mutate(all_consensus, file = paste('sequencing/sanger_final/consensus_seqs/consensus_', morph, '.fasta', sep = ''),
                         row = 1:n()) %>%
  nest(-row) %>%
  pull(data)

walk(.x = list_consensus, .f = ~seqinr::write.fasta(.x$seq, names = .x$morph, nbchar = 1000, file.out = .x$file, as.string=TRUE))

# assign taxonomy of consensus sequences using dada2
all_consensus_taxa <- all_consensus %>%
  nest(seq) %>%
  mutate(seq_no_n = map(data, clean_base_pairs)) %>%
  unnest_legacy(data, seq_no_n) %>%
  mutate(., length_no_ns = nchar(seq_no_n),
         length_seq = nchar(seq))

taxa <- dada2::assignSpecies(all_consensus_taxa$seq_no_n, 'sequencing/dada2_databases/rdp_species_assignment_16.fa.gz') %>%
  data.frame(row.names = NULL, stringsAsFactors = FALSE) %>%
  bind_cols(., all_consensus)

taxa2 <- dada2::assignTaxonomy(all_consensus_taxa$seq, 'sequencing/dada2_databases/rdp_train_set_16.fa.gz', minBoot = 40) %>%
  data.frame(row.names = NULL, stringsAsFactors = FALSE) %>%
  bind_cols(., all_consensus)

d_labs <- tibble(label = tree$tip.label) %>%
  mutate(., morph2 = case_when(substr(label,1,1) == 'b' ~ 'c',
                              substr(label,1,1) == 'd' ~ 'a',
                              substr(label,1,1) == 'c' ~ 'd',
                              substr(label,1,1) == 'e' ~ 'e',
                              substr(label,1,1) == 'f' ~ 'b')) %>%
  mutate(num_seqs = c('9/9 sequences', '1/8 sequences', '6/6 sequences', '3/8 sequences', '10/11 sequences', '1/11 sequences', '4/8 sequences', '9/9 sequences'),
         prop = c(9/9, 1/8, 6/6, 3/8, 10/11,1/11,4/8,9/9)) %>%
  mutate(., label2 = paste('morph ', morph2,'\n', num_seqs, sep = '')) %>%
  full_join(., select(taxa2, genus = Genus, label = morph)) %>%
  full_join(., select(taxa, species = Species, label = morph)) %>%
  mutate(species = ifelse(is.na(species), 'spp', species),
         label3 = paste(genus, 'spp', sep = ' '))

head(d_labs)

p_all <- ggtree(tree) %<+% d_labs + 
  geom_tiplab(aes(label=label2), size = MicrobioUoE::pts(12), offset = 0.01, family = 'Helvetica') +
  geom_tiplab(aes(label=paste0('italic(', genus,')~sp.')), 
              parse=TRUE, size = MicrobioUoE::pts(12), offset = 0.09, family = 'Helvetica') +
  geom_tippoint(aes(size =  3*prop)) +
  scale_x_continuous(expand = c(.10, .1)) +
  theme(legend.position = 'none')

ggsave('manu_figs/invasion_vs_density/consensus_seq_tree.pdf', p_all, height = 6, width =7)
ggsave('manu_figs/invasion_vs_density/consensus_seq_tree.png', p_all, height = 6, width =7)

# make summary table, not consensus sequence tree
head(d_labs)

d_table <- select(d_labs, morph2, num_seqs, prop, genus) %>%
  mutate(., isolate = case_when(morph2 == 'a' ~ 1,
                                morph2 == 'b' ~ 2,
                                morph2 == 'c' ~ 3,
                                morph2 == 'd' ~ 4,
                                morph2 == 'e' ~ 5),
         num_seqs = word(num_seqs, 1)) %>%
  select(-morph2) %>%
  arrange(isolate, -prop)

table <- select(d_table, isolate, genus, num_seqs) %>%
  flextable() %>%
  set_header_labels(isolate = 'Isolate',
                    genus = 'Sanger sequence assignment',
                    num_seqs = 'Number of clones assigned\n from disassembled communities') %>%
  align(align = 'center', part = 'all') %>% # align column names centrally
  align(align = 'left', part = 'body', j =2) %>% 
  font(fontname = 'Times', part = 'all') %>% # set font name for the table
  fontsize(size = 12, part = 'all') %>% # set font size for the table
  hline(i = c(3,5,6,7), border = fp_border_default()) %>%
  compose(j = "genus",
          value = as_paragraph(as_i(genus), ' sp.')
  ) %>% 
  bold(part = 'body', j = 2, i = c(1,4,6,7,8)) %>% # bold the column names
  autofit() # fix any random size issues

# make alternative table
d_table <- group_by(d_table, isolate) %>%
  mutate(genus2 = ifelse(prop < max(prop), 'Other', genus),
         num_seqs = parse_number(num_seqs),
         tot_seqs = sum(num_seqs)) %>%
  ungroup() %>%
  group_by(isolate, genus2) %>%
  summarise(num_seqs = sum(num_seqs),
            tot_seqs = unique(tot_seqs)) %>%
  ungroup() %>%
  mutate(num_seqs = paste(num_seqs, tot_seqs, sep = '/'))

# load in long-term stability summary data
# change isolate column to 1-5 based on alphabet of genus
d_long <- read.csv('data/processed/clone_id.csv') %>%
  mutate(isolate = case_when(isolate == 'Achromobacter' ~ 1,
                             isolate == 'Ochrobactrum' ~ 2,
                             isolate == 'Pseudomonas' ~ 3,
                             isolate == 'Stenotrophomonas' ~ 4,
                             isolate == 'Variovorax' ~ 5)) %>%
  group_by(isolate) %>%
  mutate(genus2 = ifelse(n < max(n), 'Other', genus)) %>%
  group_by(isolate, genus2) %>%
  summarise(num_seqs = sum(n),
            tot_seqs = unique(total)) %>%
  ungroup() %>%
  mutate(num_seqs = paste(num_seqs, tot_seqs, sep = '/')) %>%
  select(isolate, genus2, long_term = num_seqs)

table <- full_join(d_table, d_long) %>%
  select(-tot_seqs) %>%
  arrange(isolate, -parse_number(long_term)) %>%
  flextable() %>%
  set_header_labels(isolate = 'Isolate',
                    genus2 = 'Sequence assignment',
                    num_seqs = 'Clones assigned\nfrom disassembled\ncommunities',
                    long_term = 'Clones assigned\nfrom long-term\ncommunities') %>%
  align(align = 'center', part = 'all') %>% # align column names centrally
  align(align = 'left', part = 'body', j =2) %>% 
  font(fontname = 'Times', part = 'all') %>% # set font name for the table
  fontsize(size = 12, part = 'all') %>% # set font size for the table
  hline(i = c(2,4,6,7), border = fp_border_default()) %>%
  compose(j = "genus2",
          i = c(1,3,5,7,8),
          value = as_paragraph(as_i(genus2), ' sp.')
  ) %>% 
  bold(part = 'body', j = 2, i = c(1,3,5,7,8)) %>% # bold the column names
  autofit() # fix any random size issues

# save this table
save_as_image(table, 'figures/new/colony_table.png')
 