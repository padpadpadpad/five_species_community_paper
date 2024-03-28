# Repository for characterising the 5 species community

## Outline

This is the repository for processed data and code of the paper: "*Characterising a stable five-species microbial community for use in experimental evolution and ecology*". It mostly contains the processed datasets and code used to recreate the plots in the manuscript, as well as the raw Sanger sequencing data.

## Using this repository

The GitHub repo should be easy to use on everyone’s machine using RStudio projects. RStudio projects automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. Simply double-click on the .Rproj folder that is in the base folder of this repository and it should open the project in RStudio. All the packages used in the script should be at the top of each script, and RStudio usually requests you install them if they are not installed already.

## Data

A more detailed explanation of each data file is present in the script that uses it.

-   **assignment_16s.csv** - cleaned output of 16s taxonomic assignment from the genome assemblies.

-   **checkm.csv -** cleaned output of CheckM from the genome assemblies.

-   **checkm2.csv -** cleaned output of CheckM from the genome assemblies.

-   **gtdb_short.csv -** cleaned output of GTDBtk from the genome assemblies.

-   **interaction_data.csv -** processed interaction data.

-   **invasion_from_rare_data.csv -** processed invasion from rare data.

-   **long_term_abundance_data.csv** - processed long term abundance data.

-   **supernatant_data.csv -** processed supernatant data.

-   **invasion_from_rare_raw.csv** **-** raw invasion from rare data.

-   **long_term_assignment_16s.rds** - processed phyloseq object of the clone identification from long-term culturing of the 5 species community.

-   **sanger/**\* - the raw files from the sanger sequencing of clones.

-   **sanger/trimmed/\*** - the trimmed files from the sanger sequencing of clones.

## Scripts

-   **sanger_sequence_clean_5spp.R -** sets trimming parameters for the raw sanger sequencing files in **sanger/**.
-   **sanger_sequence_analysis_5spp.R -** trims the sanger sequencing files, assigns taxonomy using **dada2**, and recreates the table in Figure 1c.
-   **invasion_from_rare.R** - analysis of invasion from rare data and recreates Figure 2.
-   **phenotypic_assays.R** - analyses the phenotypic data (supernatant assays, invasion-from-rare assays, medium-term persistence assays) and recreates Figures 3, and 4.
-   **long_term_clone_id.R** - analyses **long_term_assignment_16s.rds** to look at how well the 5 species can be identified after \>1 year in co-culture.
-   **make_genome_summary.R** - creates the genome summary table (Table 1) from **assignment_16s.csv**, **checkm.csv**, **checkm2.csv**, and **gtdb_short.csv**.
