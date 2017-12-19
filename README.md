TISSUES 2.0: An integrative web resource on mammalian tissue expression
==============

Code to reproduce the fold enrichment analyses and figures from the article
--------------
**Project structure**
- README.md --> *Markdown file*
- makefile  -->	*main script*
- generate\_files.pl --> *script to generate the necessary files. This script uses the original datasets files (with no filter for unconfident gene-tissue associations)*
- add\_orthology\_info.pl --> *script to generate files containing scored gene-tissue pairs and 1:1 orthologs with the other three organisms, if these exist, for each gene.*
- analyses.R --> *script orchestrating the generation of the figures*
- R/ --> *All the scripts necessary to reproduce the figures*
    - summary\_figure\_drawing.R --> *Initial figure with the tissues covered by each dataset*
    - fold\_enrichment\_analysis\_per\_dataset.R --> *Generates the fold-enrichment plots for each dataset*
    - fold\_enrichment\_score\_calibration\_analysis.R --> *Score calibration figure*
    - correlation\_analysis.R --> *Heatmap figures based on Pearson's correlation coefficients between final star scores for gene tissue pairs across datasets*
    - functions.R --> *functions used across the other four scripts* 
- data/ 
 - datasests/ --> *Datasets files need to be stored here*
 - dictionary/ --> *Contains the files necessary to perform the tissue backtracking bassed on BRENDA Ontology*
  - labels.tsv --> *The bto terms corresponding to the 21 tissues of interest: tissues\_code  tissue\_name  BTO*
  - bto\_entities.tsv --> *mapping of bto terms to internal identifiers: internal\_code  tissues\_code  BTO*
  - bto\_groups .tsv --> *the parent-children relationships used to do the backtracking: internal\_code  parent\_internal\_code*
 - orthology/ --> *Contains the eggNOG orthology groups files necessary to extract 1:1 orthologs across organisms*
  - roNOG.tsv --> *ortholgy groups for rodents, used for extracting 1:1 orthologs between mouse and rat*
  - maNOG.tsv --> *ortholgy groups for mammals, used for extracting orthologs between the other organism pairs*
- figures/ --> *Folder where all the figures generated are stored*

**Run the analyses**

1. Download the project
2. Make sure you have a default CRAN repository set, for example by
   putting the following into your `~/.Rprofile`:

     ```r
     local({r <- getOption("repos")
            r["CRAN"] <- "http://cran.us.r-project.org"
            options(repos=r)})
     ```

3. Download the datasets from https://doi.org/10.6084/m9.figshare.5715880 and place them in the data/datasets folder 
4. Execute the makefile script from the command line:
  `> make`
5. All the files will be generated in the data folder
6. All the figures from the analyses will be created in the figures/ folder

**Generated files**

./data/
- Fold enrichment analyses result files: *DATASET*\_uniprot\_fold\_enrichment\_analysis.tsv
- Scored gene-tissue pairs for the 21 major tissues: pairs\_major\_tissues.tsv 
./figures/
- Summary of covered tissues: Summary\_figure\_tissues.png
- Fold enrichment plots: datasets\_fold\_enrichment.png
- Score calibration plot: datasets_score_calibration.png
- Heatmaps with Pearson's correlation scores between pairs of datasets 

**Requirements**

- Perl
- R
