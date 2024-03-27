#-----------------------------#
# genome table information ####
#-----------------------------#

# make genome summary table of the five species community
librarian::shelf(flextable, officer, ftExtra, tidyverse)

# set path to reference genome folder
path <- 'five_spp_genomes'

d_16s <- read.csv('summaries/assignment_16s.csv') %>%
  group_by(file)%>%
  slice_min(., dist_from_ref_16s) %>%
  mutate(species = replace_na(species, ''),
         assignment_16s = trimws(paste(genus, species, sep = ' '))) %>%
  select(sample = file, assignment_16s, num_16s) 

d_gtdbtk <- read.csv('summaries/gtdb_short.csv')

d_checkm2 <- read.csv('summaries/checkm2.csv')
d_checkm <- read.csv('summaries/checkm.csv')

d_taxa <- left_join(d_gtdbtk, d_16s) %>%
  left_join(., d_checkm2) %>%
  left_join(., select(d_checkm, sample, n_contigs)) %>%
  mutate(community_member = 1:n(),
         species2 = c('Achromobacter veterisilvae', 'Ochrobactrum teleogrylli', 'Pseudomonas fluorescens', 'Stenotrophomonas', 'Variovorax'),
         species3 = c('AB1', 'AB1', 'AB1', 'sp. AB1 (2024)', 'sp. AB1 (2024)'))

d_table <- select(d_taxa, community_member, sample, species, species2, species3, num_16s, genome_size, n_contigs, GC, total_coding_sequences, completeness, contamination) %>%
  mutate(across(genome_size:contamination, \(x) round(x, 2)))

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

save_as_image(table1, 'tables/paper_table.png', webshot = 'webshot2')

d_table <- select(d_taxa, community_member, sample, species, assignment_16s, ani, align_frac, num_16s, genome_size, n_contigs, N50_contig, GC, total_coding_sequences, completeness, contamination) %>%
  mutate(across(genome_size:contamination, \(x) round(x, 2)))

table1 <- flextable(select(d_table, community_member:num_16s)) %>%
  set_header_labels(community_member = 'isolate',
                    sample = 'Sanger sequence assignment',
                    species = 'GTDBtk assignment',
                    ani = 'ANI',
                    align_frac = 'Aligned fraction',
                    num_16s = '16s copy\nnumber') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  fix_border_issues() %>%
  italic(j = c('species', 'assignment_16s')) %>%
  autofit() %>%
  bg(bg = 'white', part = 'all')
table1

save_as_image(table1, 'five_spp_genomes/tables/taxa_assignment.png', webshot = 'webshot2')

table1 <- flextable(select(d_table, community_member, sample, genome_size:contamination)) %>%
  set_header_labels(community_member = 'isolate',
                    sample = 'Sanger sequence assignment',
                    genome_size = 'Genome size',
                    n_contigs = 'Number of contigs',
                    N50_contig = 'N50 contig',
                    GC = 'GC content',
                    total_coding_sequences = 'total coding sequences',
                    completeness = 'Completeness',
                    contamination = 'Contamination') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  fix_border_issues() %>%
  autofit() %>%
  bg(bg = 'white', part = 'all')

table1

save_as_image(table1, 'five_spp_genomes/tables/genomes_stats.png', webshot = 'webshot2')

# make a table of the amr genes
d_amr <- read.csv('five_spp_genomes/summaries/amrfinderplus.csv')

table_amr <- select(d_amr, n, sample, contig_id, start, stop, subclass, sequence_name, percent_identity_to_reference_sequence) %>%
  flextable() %>%
  set_header_labels(n = 'isolate',
                    sample = 'Sanger sequence assignment',
                    contig_id = 'contig ID',
                    start = 'start',
                    stop = 'stop',
                    subclass = 'Subclass',
                    sequence_name = 'Sequence name',
                    percent_identity_to_reference_sequence = '% identity to ref') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  fix_border_issues() %>%
  autofit() %>%
  bg(bg = 'white', part = 'all')

save_as_image(table_amr, 'five_spp_genomes/tables/amr_genes.png', webshot = 'webshot2')

# have a look at padloc output
d_padloc <- read.csv('five_spp_genomes/summaries/padloc.csv')
colnames(d_padloc)

d_table <- group_by(d_padloc, sample, system) %>%
  tally() %>%
  group_by(sample) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  select(sample, total, system, n)

table <- flextable(d_table) %>%
  set_header_labels(sample = 'Sanger sequence assignment',
                    total = "Total systems ID'ed",
                    system = 'System name',
                    n = 'Number identified') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  merge_v(j = c('sample', 'total')) %>%
  valign(j = c(1, 2), valign = 'top') %>%
  fix_border_issues() %>%
  autofit() %>%
  bg(bg = 'white', part = 'all')

table

save_as_image(table, 'five_spp_genomes/tables/padloc.png', webshot = 'webshot2')

# genomad
d_genomad <- read.csv('five_spp_genomes/summaries/virus_info.csv')

table <- select(d_genomad, n, sample, topology, length, coordinates, virus_score, n_hallmarks, marker_enrichment) %>%
  flextable() %>%
  set_header_labels(n = 'Isolate',
                    sample = 'Sanger sequence\nassignment',
                    topology = "Topology",
                    length = 'Sequence\nlength',
                    coordinates = 'Provirus\nposition',
                    virus_score = 'geNomad confidence\nscore',
                    n_hallmarks = 'Number of hallmark\ngeNomad markers',
                    marker_enrichment = 'Enrichment of\nviral markers') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  fix_border_issues() %>%
  autofit() %>%
  bg(bg = 'white', part = 'all')

save_as_image(table, 'five_spp_genomes/tables/virus.png', webshot = 'webshot2')

# genomad
d_genomad <- read.csv('five_spp_genomes/summaries/plasmid_info.csv')

table <- select(d_genomad, n, sample, topology, n_genes, length, plasmid_score, n_hallmarks, marker_enrichment, conjugation_genes) %>%
  mutate(conjugation_genes = str_count(conjugation_genes, ';') + 1) %>%
  flextable() %>%
  set_header_labels(n = 'Isolate',
                    sample = 'Sanger sequence\nassignment',
                    topology = "Topology",
                    length = 'Sequence\nlength',
                    coordinates = 'Provirus\nposition',
                    plasmid_score = 'geNomad confidence\nscore',
                    n_genes = 'Number of\ngenes',
                    n_hallmarks = 'Number of hallmark\ngeNomad markers',
                    marker_enrichment = 'Enrichment of\nplasmid markers',
                    conjugation_genes = 'NUmber of\nconjugation genes') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  fix_border_issues() %>%
  autofit() %>%
  bg(bg = 'white', part = 'all')

save_as_image(table, 'five_spp_genomes/tables/plasmid.png', webshot = 'webshot2')

fegenie <- read.csv('five_spp_genomes/fegenie/FeGenie-geneSummary.csv', skip = 1) %>%
  group_by(genome.assembly, category) %>%
  tally() %>%
  mutate(genome.assembly = gsub('.fasta', '', genome.assembly)) %>%
  ungroup()

table <- 
  flextable(fegenie) %>%
  set_header_labels(genome.assembly = 'Sanger sequence\nassignment',
                    category = "FeGenie gene\ncategory",
                    n = 'Number of\ngenes') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  align(align = 'center', part = 'all') %>%
  merge_v(j = c('genome.assembly')) %>%
  valign(j = c(1), valign = 'top') %>%
  fix_border_issues() %>%
  autofit() %>%
  bg(bg = 'white', part = 'all')

save_as_image(table, 'five_spp_genomes/tables/fegenie.png', webshot = 'webshot2')
