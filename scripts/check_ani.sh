# look at how similar genomes are

# look at all the genomes for P fluorescens on ncbi

# use ncbi-genome-download to get the genomes

# create a new environment and install ncbi-genome-download in it
mamba create -n ncbi_genome_download -c bioconda ncbi-genome-download

# activate the environment
micromamba activate ncbi_genome_download

# move to the right directory
cd ~/google_drive/work_googledrive/five_species_community/coexisting_community_paper/github/data/genome_comparison

ncbi-genome-download --genera "Pseudomonas fluorescens" bacteria --dry-run

# download the genomes
ncbi-genome-download --genera "Pseudomonas fluorescens SBW25" bacteria --assembly-levels "all" --dry-run --formats fasta --section "genbank"

# download the genomes
ncbi-genome-download --genera "Pseudomonas fluorescens SBW25" bacteria --assembly-levels "all" --formats "fasta,assembly-report" --section "genbank" --parallel 4

# close the environment
micromamba deactivate

# create a new environment called fastani and install bioconda/fastani in it
mamba create -n fastani

# activate the environment
micromamba activate fastani
conda install bioconda::fastani

fastani -h
fastani -q ~/google_drive/work_googledrive/five_species_community/five_spp_genomes/assemblies/Pseudomonas.fasta -r ~/google_drive/work_googledrive/five_species_community/coexisting_community_paper/github/data/genome_comparison/genbank/bacteria/GCA_931907645.1/GCA_931907645.1_MPBAS00001_genomic.fna.gz -o ~/google_drive/work_googledrive/five_species_community/coexisting_community_paper/github/data/genome_comparison/genbank/bacteria/comparison.txt

/GCA_000012445.1_genomic.fna.gz

# run fastani on two genomes

mamba create -n skani
micromamba activate skani
conda install bioconda::skani

skani dist ~/google_drive/work_googledrive/five_species_community/five_spp_genomes/assemblies/Pseudomonas.fasta ~/google_drive/work_googledrive/five_species_community/coexisting_community_paper/github/data/genome_comparison/genbank/bacteria/GCA_931907645.1/GCA_931907645.1_MPBAS00001_genomic.fna.gz

skani dist ~/google_drive/work_googledrive/five_species_community/five_spp_genomes/assemblies/Pseudomonas.fasta ~/Downloads/GCF_902498135.1_PS900_genomic.fna

skani dist ~/google_drive/work_googledrive/five_species_community/coexisting_community_paper/github/data/genome_comparison/genbank/bacteria/GCA_931907645.1/GCA_931907645.1_MPBAS00001_genomic.fna.gz ~/Downloads/GCF_902498135.1_PS900_genomic.fna
