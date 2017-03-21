.PHONY : all
.DEFAULT: all
PACKAGE      = package
VERSION      = ` date "+%Y.%m%d%" `
RELEASE_DIR  = ..
RELEASE_FILE = $(PACKAGE)-$(VERSION)
SILENT = yes
DATA_DIR     = data/datasets

# download the datasets from https://figshare.com/articles/datasets/4640734 and place them in the data/datasets folder

all: generate_files generate_figures

#Generates all the files necessary for the analyses
generate_files: 
	perl generate_files.pl;
	perl update_consistency_with_orthologs.pl
	#Generates all the figures from the analyses described in the article: TISSUES: An integrative web resource on mammalian tissue expression
generate_figures:
	Rscript analyses.R




